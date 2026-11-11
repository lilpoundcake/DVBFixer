#!/usr/bin/env python3
"""Map local RFD3 outputs into a retained assembly without inventing a sampler trace.

Native exports remain intact. Observed-only fitting precedes coordinate copying;
only declared newly sampled atoms replace the corresponding retained donors.
Other generated donors are post-sampling assembly context, never conditioning.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_group(inputs, run, *, mapping_recovery=None):
    manifest = json.loads((inputs / "manifest.json").read_text())
    predictions = list((run / "predictions").glob("*.cif.gz"))
    if len(predictions) != 1:
        raise ValueError("missing/ambiguous local native export")
    prediction = predictions[0]
    metadata = json.loads(prediction.with_suffix("").with_suffix(".json").read_text())
    atom_map = metadata["diffused_index_map"]
    recovery_receipt = None
    if mapping_recovery is not None:
        recovery = json.loads(mapping_recovery.read_text())
        if (recovery["input_manifest_sha256"] != digest(inputs / "manifest.json") or
                recovery["original_metadata_sha256"] != digest(prediction.with_suffix("").with_suffix(".json")) or
                not recovery["same_coordinates_sequence_mask_atom_axis"]):
            raise ValueError("explicit mapping recovery differs from its original input/output")
        restored_map = recovery["diffused_index_map"]
        if any(restored_map[k] != value for k, value in atom_map.items() if k in restored_map):
            raise ValueError("mapping recovery contradicts an intact upstream identity")
        atom_map = restored_map
        recovery_receipt = {"path": str(mapping_recovery), "sha256": digest(mapping_recovery)}
    native = inventory(prediction)
    source = inventory(inputs / "source-projection.pdb")
    selected = {}
    expected_native = set()
    for item in manifest["axis"]:
        atom = AtomIdentity(**item["atom"])
        label = f"{atom.chain}{atom.residue_number}"
        if label not in atom_map:
            raise ValueError(f"RFD3 source identity {label} is absent/truncated; explicit verified mapping recovery required")
        match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", atom_map[label])
        if match is None:
            raise ValueError("unsupported RFD3 output identity mapping")
        chain, number, code = match.groups()
        mapped = AtomIdentity(chain, number, code, atom.atom_name)
        if mapped in expected_native or mapped not in native or native[mapped][0] != item["residue_name"]:
            raise ValueError("local sequence/atom-axis mapping is incomplete or ambiguous")
        expected_native.add(mapped)
        selected[atom] = native[mapped][2]
    if expected_native != native.keys():
        raise ValueError("local native export has unrequested/missing atoms")
    fixed = tuple(sorted(source))
    transform = weighted_kabsch(np.array([selected[a] for a in fixed]), np.array([source[a][2] for a in fixed]))
    aligned = {a: transform.apply(np.asarray([xyz]))[0] for a, xyz in selected.items()}
    drift = np.array([np.linalg.norm(aligned[a] - source[a][2]) for a in fixed])
    generated = {a: xyz for a, xyz in aligned.items() if a not in source}
    if len(generated) != manifest["generated_atoms"]:
        raise ValueError("generated mask differs after local mapping")
    return generated, {"native_export": str(prediction), "native_sha256": digest(prediction),
        "input_manifest_sha256": digest(inputs / "manifest.json"), "metadata_sha256": digest(prediction.with_suffix("").with_suffix(".json")),
        "native_observed_drift_rms_angstrom": float(np.sqrt(np.mean(drift**2))),
        "native_observed_drift_max_angstrom": float(drift.max()), "sampled_atoms": len(generated),
        "atomic_gap_indices": manifest["selected_gap_indices"], "rotation": transform.rotation.tolist(),
        "translation": transform.translation.tolist(), "source_sha256": manifest["source_sha256"],
        "explicit_mapping_recovery": recovery_receipt, "native_observations_within_pdb_precision": bool(drift.max() <= 1e-3)}


def run(root, output, *, context_only=False, mapping_recovery=None):
    base = root / "1he8-visual-compatible-selection"
    request = DiffusionRequest.from_json((base / "evaluation-mask.json").read_text())
    source = (base / "source.pdb").read_text()
    baseline = (base / "selected.pdb").read_text()
    if hashlib.sha256(source.encode()).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("full source differs from original request")
    if context_only:
        all_groups = {"context": native_group(root / "1he8-rfd3-pair-context-input", root / "1he8-rfd3-pair-context",
                                             mapping_recovery=mapping_recovery)}
        assemblies = (("joint-context", ("context",)),)
    else:
        all_groups = {name: native_group(root / f"1he8-rfd3-{name}-input", root / f"1he8-rfd3-local-controls/{name}")
                      for name in ("pair", "gap2", "gap3")}
        assemblies = (("joint-pair", ("pair",)), ("independent-pair", ("gap2", "gap3")))
    output.mkdir(exist_ok=False)
    for name, arms in assemblies:
        work = output / name
        work.mkdir()
        coordinates = {}
        lineage = []
        for arm in arms:
            xyz, receipt = all_groups[arm]
            if receipt["source_sha256"] != request.normalized_pdb.sha256 or coordinates.keys() & xyz.keys():
                raise ValueError("source/group ownership differs or overlaps")
            coordinates.update(xyz)
            lineage.append(receipt)
        wanted = {a for a in request.generated_atoms if any(
            (a.chain, a.residue_number, a.insertion_code) == (r.chain, r.residue_number, r.insertion_code)
            for i in (2, 3) for r in request.gaps[i].generated_residues)}
        if coordinates.keys() != wanted:
            raise ValueError("local proposals do not cover exactly both original gaps")
        text = rewrite_generated_coordinates(baseline, coordinates)
        assessment = assessment_adapter(request, source)(text)
        try:
            structure = _structure(text, "local-pair.pdb")
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
            chirality = {"all_l": True}
        except ChiralityError as exc:
            chirality = {"all_l": False, "error": str(exc)}
        (work / "raw.pdb").write_text(text)
        (work / "source-restored.pdb").write_text(text)
        (work / "source.pdb").write_text(source)
        (work / request.normalized_pdb.path).write_text(source)
        (work / "evaluation-mask.json").write_text(request.to_json())
        report = {"arm": name, "new_backend_groups": lineage, "retained_other_eight_donors_sha256": digest(base / "selected.pdb"),
            "retained_selection_recipe_sha256": digest(base / "result.json"), "new_sampled_atoms": len(coordinates),
            "whole_assembly_generated_atoms": len(request.generated_atoms), "other_donors_used_only_after_sampling": True,
            "assessment": asdict(assessment), "geometry_passed": assessment.geometry_passed, "chirality": chirality,
            "public_sampler_trace_not_fabricated": True, "scientific_qualification": False}
        (work / "lineage-and-assessment.json").write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps({"arm": name, "actionable_findings": sum(i.actionable for i in assessment.issues), "all_l": chirality["all_l"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--context-only", action="store_true")
    parser.add_argument("--mapping-recovery", type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute(), context_only=args.context_only, mapping_recovery=args.mapping_recovery)
