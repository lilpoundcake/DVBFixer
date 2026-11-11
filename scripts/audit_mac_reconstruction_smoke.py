#!/usr/bin/env python3
"""Independent post-sampling audit of explicitly declared reference projections.

Reference coordinates enter only held-out scoring after fitting to observations.
Native drift is retained before any source restoration. No sampling trace is
fabricated and no refinement or acceptance threshold is changed here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import gemmi
import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request


def inventory(path):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("expected exactly one model")
    result = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                identity = AtomIdentity(chain.name, str(residue.seqid.num), residue.seqid.icode.strip(), atom.name)
                if identity in result or atom.altloc not in ("\x00", " ") or atom.occ != 1:
                    raise ValueError("duplicate/alternate/partial atom identity")
                result[identity] = (residue.name, atom.element.name, np.array(list(atom.pos)))
    return result


def audit(input_dir, prediction, output, *, seed=7):
    manifest = json.loads((input_dir / "manifest.json").read_text())
    for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb"):
        if hashlib.sha256((input_dir / name).read_bytes()).hexdigest() != manifest.get("files", {}).get(name):
            raise ValueError("declared source/reference projection bytes changed")
    source_path = input_dir / "source.pdb"
    source = inventory(source_path)
    reference = inventory(input_dir / "reference-projection-NOT-CONDITIONING.pdb")
    native = inventory(prediction)
    if set(native) != set(reference) or any(native[k][:2] != reference[k][:2] for k in reference):
        raise ValueError("native sequence/atom inventory differs from declared projection")
    if not all(np.isfinite(value[2]).all() for value in native.values()):
        raise ValueError("nonfinite native coordinates")
    request = build_cli_diffusion_request(source_path, {manifest.get("chain", "A"): manifest["sequence"]}, seeds=(seed,),
                                         profile="protpardelle-1c-cc94-mps-v1")
    # This is an evaluation mask only, not a claim that another model used the
    # Protpardelle sampler/profile or generated its trace.
    axis = tuple(sorted(native))
    fixed = tuple(sorted(source))
    transform = weighted_kabsch(np.array([native[k][2] for k in fixed]), np.array([source[k][2] for k in fixed]))
    aligned = transform.apply(np.array([native[k][2] for k in axis]))
    positions = dict(zip(axis, aligned))
    drift = np.array([np.linalg.norm(positions[k] - source[k][2]) for k in fixed])
    generated = set(axis) - set(source)
    if generated != set(request.generated_atoms):
        raise ValueError("evaluation generated mask differs from source inventory")
    restored = np.array([source[k][2] if k in source else positions[k] for k in axis])
    source_text = source_path.read_bytes().decode()
    names = tuple(native[k][0] for k in axis)
    elements = tuple(native[k][1] for k in axis)
    aligned_text = materialize_candidate_pdb(source_text, axis, aligned, names, request, elements=elements)
    restored_text = materialize_candidate_pdb(source_text, axis, restored, names, request, elements=elements)
    assess = assessment_adapter(request, source_text)
    findings = {}
    for label, text in (("aligned-native", aligned_text), ("source-restored", restored_text)):
        result = assess(text)
        findings[label] = asdict(result)
        findings[label]["geometry_passed"] = result.geometry_passed
    error = np.array([np.linalg.norm(positions[k] - reference[k][2]) for k in sorted(generated)])
    junctions = {}
    boundaries = []
    for gap in request.gaps:
        if gap.left_anchor is not None:
            boundaries.append((gap.left_anchor, gap.generated_residues[0]))
        if gap.right_anchor is not None:
            boundaries.append((gap.generated_residues[-1], gap.right_anchor))
    for label, coords in (("native", positions), ("source-restored", dict(zip(axis, restored)))):
        junctions[label] = {f"{left.chain}{left.residue_number}{left.insertion_code}-{right.chain}{right.residue_number}{right.insertion_code}": float(np.linalg.norm(
            coords[AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")] -
            coords[AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")]))
            for left, right in boundaries}
    report = {"prediction": str(prediction), "prediction_sha256": hashlib.sha256(prediction.read_bytes()).hexdigest(),
        "input_manifest_sha256": hashlib.sha256((input_dir / "manifest.json").read_bytes()).hexdigest(),
        "evaluation_mask_only": True, "sampled_profile_not_inferred": True,
        "exact_sequence_and_atom_inventory": True, "observed_atoms": len(source), "generated_atoms": len(generated),
        "native_observed_drift_after_observed_only_fit_angstrom": {"rms": float(np.sqrt(np.mean(drift**2))), "max": float(drift.max())},
        "heldout_generated_rmsd_angstrom": float(np.sqrt(np.mean(error**2))),
        "junction_cn_angstrom": junctions, "assessments": findings,
        "scientific_qualification": False, "new_refinement": False}
    output.mkdir(parents=True, exist_ok=False)
    (output / request.normalized_pdb.path).write_text(source_text)
    (output / "evaluation-mask.json").write_text(request.to_json())
    (output / "aligned-native.pdb").write_text(aligned_text)
    (output / "source-restored.pdb").write_text(restored_text)
    np.savez_compressed(output / "coordinates.npz", native=np.array([native[k][2] for k in axis]), aligned=aligned,
                        restored=restored, rotation=transform.rotation, translation=transform.translation)
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    print(json.dumps(audit(args.input, args.prediction, args.output, seed=args.seed), indent=2))
