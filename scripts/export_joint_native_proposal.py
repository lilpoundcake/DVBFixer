#!/usr/bin/env python3
"""Independently reconstruct a historical local-union joint proposal; no refinement."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path


def main():
    import numpy as np

    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.joint_gap_research import (
        STRUCTURE_REVISION,
        build_joint_structure_view,
    )
    from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.sample
    trace = json.loads((root / "trace.json").read_text())
    if (trace["kind"] != STRUCTURE_REVISION or trace["status"] != "native-proposal-retained"
            or trace["sample_invocations"] != 1 or trace["ode_updates"] != 500):
        raise ValueError("not a completed whole-gap local-union trajectory")
    parent = DiffusionRequest.from_json((root / "parent-request.json").read_text())
    source = (root / "source.pdb").read_bytes()
    baseline = (root / "baseline-output-only.pdb").read_bytes()
    view = build_joint_structure_view(parent, source, seed=trace["seed"])
    if (trace["gap_indices"] != list(range(len(parent.gaps)))
            or (root / "research-group.json").read_bytes() != view.request.to_json().encode()
            or hashlib.sha256(parent.to_json().encode()).hexdigest() != trace["parent_request_sha256"]
            or hashlib.sha256(source).hexdigest() != trace["source_sha256"]
            or hashlib.sha256(baseline).hexdigest() != trace["baseline_output_only_sha256"]
            or hashlib.sha256((root / "native.npz").read_bytes()).hexdigest() != trace["native_sha256"]):
        raise ValueError("joint input/native custody differs")
    with np.load(root / "native.npz", allow_pickle=False) as archive:
        if set(archive.files) != {"atom37"}:
            raise ValueError("unexpected native arrays")
        native = archive["atom37"]
    if native.shape != (len(view.positions), 37, 3) or not np.isfinite(native).all():
        raise ValueError("native shape/coordinates differ")
    axis = {}
    cursor = 0
    for target, placement in zip(view.engine_view.request.target_sequences,
                                 view.engine_view.request.sequence_placements, strict=True):
        mapped = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for gap in view.engine_view.request.gaps:
            if gap.chain == target.chain:
                mapped.update(zip(range(gap.target_interval.start, gap.target_interval.stop),
                                  gap.generated_residues, strict=True))
        if set(mapped) != set(range(len(target.sequence))):
            raise ValueError("native axis has missing tokens")
        axis.update({(r.chain, r.residue_number, r.insertion_code): cursor + i for i, r in mapped.items()})
        cursor += len(target.sequence)

    def point(atom):
        return native[axis[atom.chain, atom.residue_number, atom.insertion_code], ATOM37_ORDER.index(atom.atom_name)]

    source_xyz = _coordinates(source.decode())[0]
    mapping = view.engine_view.atom_to_source
    pairs = [(e, s) for e, s in mapping if s in source_xyz and e.atom_name != "OXT"]
    fit = weighted_kabsch(np.asarray([point(e) for e, _ in pairs]), np.asarray([source_xyz[s] for _, s in pairs]))
    generated = set(view.request.generated_atoms)
    coordinates = {s: fit.apply(np.asarray([point(e)]))[0] for e, s in mapping if s in generated}
    if set(coordinates) != set(parent.generated_atoms):
        raise ValueError("proposal does not replace every generated atom")
    with np.load(root / "generated.npz", allow_pickle=False) as archive:
        if set(archive.files) != {"coordinates"} or not np.allclose(archive["coordinates"],
                np.asarray([coordinates[a] for a in view.request.generated_atoms]), rtol=0, atol=1e-8):
            raise ValueError("independent native proposal differs")
    raw = rewrite_generated_coordinates(baseline.decode(), coordinates)
    assessment = assessment_adapter(parent, source.decode())(raw)
    if assessment.admission_failures:
        raise ValueError("raw proposal changed immutable source/chemistry")
    filename = f'1HE8_all-gaps_seed-{trace["seed"]}_raw.pdb'
    with (args.output / filename).open("xb") as stream:
        stream.write(raw.encode())
    report = {"seed": trace["seed"], "file": filename,
        "sha256": hashlib.sha256(raw.encode()).hexdigest(), "tokens": len(view.positions),
        "generated_gaps": len(parent.gaps), "generated_atoms": len(coordinates),
        "independent_native_materialization_matches": True, "observations_unchanged": True,
        "baseline_generated_coordinates_all_replaced": True, "new_sampling_during_export": False,
        "refined": False, "scientific_qualification": False, "published": False}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "raw-assessment.json").write_text(json.dumps(
        asdict(assessment), ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
