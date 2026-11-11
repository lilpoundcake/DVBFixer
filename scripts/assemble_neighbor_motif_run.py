#!/usr/bin/env python3
"""Independently assemble native proposals from one call per neighbor component.

No component refinement, no sampling. Replaces every generated atom in a common
observed frame; any retained patchwork is an output record template only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main():
    import numpy as np

    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.gap_neighbors import source_gap_components
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.joint_gap_research import (
        MOTIF_COMPONENT_REVISION,
        plan_motif_component,
    )
    from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((args.run / "declaration.json").read_text())
    first = args.run / "sampling-component-0"
    parent_bytes = (first / "parent-request.json").read_bytes()
    source_bytes = (first / "source.pdb").read_bytes()
    baseline = (first / "baseline-output-only.pdb").read_bytes()
    parent = DiffusionRequest.from_json(parent_bytes.decode())
    components, _ = source_gap_components(parent, source_bytes.decode())
    assert [list(g) for g in components] == manifest["components"]
    source_xyz = _coordinates(source_bytes.decode())[0]
    proposed = {}
    receipts = []
    for i, group in enumerate(components):
        root = args.run / f"sampling-component-{i}"
        for name, expected in (("parent-request.json", parent_bytes), ("source.pdb", source_bytes),
                               ("baseline-output-only.pdb", baseline)):
            if (root / name).read_bytes() != expected:
                raise ValueError("component custody differs")
        trace = json.loads((root / "trace.json").read_text())
        if (trace["kind"] != MOTIF_COMPONENT_REVISION or trace["status"] != "native-proposal-retained"
                or trace["sample_invocations"] != 1 or trace["ode_updates"] != 500
                or trace["gap_indices"] != list(group) or trace["seed"] != manifest["seed"]):
            raise ValueError("component is not one completed shared diffusion trajectory")
        view = plan_motif_component(parent, source_bytes, group, seed=trace["seed"])
        if (root / "research-group.json").read_bytes() != view.request.to_json().encode():
            raise ValueError("component context/mask differs")
        native_bytes = (root / "native.npz").read_bytes()
        if hashlib.sha256(native_bytes).hexdigest() != trace["native_sha256"]:
            raise ValueError("component native digest differs")
        with np.load(root / "native.npz", allow_pickle=False) as archive:
            assert set(archive.files) == {"atom37"}
            native = archive["atom37"]
        if native.shape != (len(view.positions), 37, 3) or not np.isfinite(native).all():
            raise ValueError("component native shape/values differ")
        axis = {}
        cursor = 0
        for target, placement in zip(view.engine_view.request.target_sequences,
                                     view.engine_view.request.sequence_placements, strict=True):
            mapped = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
            for gap in view.engine_view.request.gaps:
                if gap.chain == target.chain:
                    mapped.update(zip(range(gap.target_interval.start, gap.target_interval.stop),
                                      gap.generated_residues, strict=True))
            assert set(mapped) == set(range(len(target.sequence)))
            axis.update({(r.chain, r.residue_number, r.insertion_code): cursor + index for index, r in mapped.items()})
            cursor += len(target.sequence)

        def point(atom, *, native=native, axis=axis):
            return native[axis[atom.chain, atom.residue_number, atom.insertion_code], ATOM37_ORDER.index(atom.atom_name)]

        mapping = view.engine_view.atom_to_source
        fixed = [(e, s) for e, s in mapping if s in source_xyz and e.atom_name != "OXT"]
        fit = weighted_kabsch(np.asarray([point(e) for e, _ in fixed]), np.asarray([source_xyz[s] for _, s in fixed]))
        generated = set(view.request.generated_atoms)
        coordinates = {s: fit.apply(np.asarray([point(e)]))[0] for e, s in mapping if s in generated}
        if set(coordinates) != generated or set(proposed) & generated:
            raise ValueError("component ownership overlaps or omits generated atoms")
        with np.load(root / "generated.npz", allow_pickle=False) as archive:
            assert set(archive.files) == {"coordinates"}
            if not np.allclose(archive["coordinates"], np.asarray([coordinates[a] for a in view.request.generated_atoms]),
                               rtol=0, atol=1e-8):
                raise ValueError("independent component native materialization differs")
        proposed.update(coordinates)
        receipts.append({"gap_indices": list(group), "tokens": len(view.positions),
            "generated_atoms": len(generated), "sample_invocations": 1, "ode_updates": 500,
            "independent_native_materialization_matches": True})
    if set(proposed) != set(parent.generated_atoms):
        raise ValueError("assembly does not own every generated atom exactly once")
    raw = rewrite_generated_coordinates(baseline.decode(), proposed)
    assessment = assessment_adapter(parent, source_bytes.decode())(raw)
    if assessment.admission_failures:
        raise ValueError("assembly violates immutable source contract")
    for name, data in (("raw-assembled.pdb", raw.encode()), ("request.json", parent_bytes),
                       ("source.pdb", source_bytes), ("run-script.py", Path(__file__).read_bytes())):
        with (args.output / name).open("xb") as stream:
            stream.write(data)
    report = {"components": receipts, "new_sampling": False, "no_component_refinement": True,
        "all_generated_atoms_replaced": len(proposed), "observations_unchanged": True,
        "raw_sha256": hashlib.sha256(raw.encode()).hexdigest(), "scientific_qualification": False}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
