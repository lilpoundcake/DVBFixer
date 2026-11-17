#!/usr/bin/env python3
"""Export independently checked private heavy-coordinate candidates for review."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest


def export(analysis, readback, output):
    receipt_path = readback / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    for name, expected in receipt["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"independent readback payload changed: {name}")
    completed_path = analysis / "completed.json"
    if receipt["bindings"].get(str(completed_path)) != digest(completed_path):
        raise ValueError("independent readback is for another analysis")
    completed = json.loads(completed_path.read_text())
    analysis_plan = analysis / "plan.json"
    if receipt["bindings"].get(str(analysis_plan)) != digest(analysis_plan):
        raise ValueError("independent readback lacks the analysis declaration")
    for name, expected in json.loads(analysis_plan.read_text())["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"analysis lineage payload changed: {name}")
    selected = [row for row in receipt["rows"] if row["accepted"]]
    if not selected or len({row["entry_id"] for row in selected}) != len(selected):
        raise ValueError("require unique independently accepted entries; no donor ranking")
    output.mkdir(parents=True, exist_ok=False)
    exports = []
    for evidence in selected:
        if (evidence["changed_native_observations"] or evidence["new_geometry_findings"] or not evidence["global_all_l"]
                or not evidence["canonical_dat_ownership_exact"] or not evidence["original_full_context_record_subsequence_exact"]):
            raise ValueError("incomplete independent heavy-coordinate acceptance")
        result = next(row for row in completed["rows"] if row["entry_id"] == evidence["entry_id"] and row["stage"] == evidence["stage"])
        directory, case = Path(result["directory"]), Path(result["case"])
        lineage = json.loads((case / "lineage.json").read_text())
        original_cif = Path(lineage["original_source"])
        if digest(original_cif) != lineage["original_source_sha256"]:
            raise ValueError("original deposited source changed")
        destination = output / evidence["entry_id"]
        destination.mkdir(exist_ok=False)
        inputs = {
            f"{evidence['entry_id']}-restored-heavy.pdb": directory / "full-context.pdb",
            f"{evidence['entry_id']}-restored-heavy.dat": directory / "generated.dat",
            "original-deposited-source.cif.gz": original_cif,
            "original-full-context.pdb": case / "original-full-context.pdb",
            "prepared-full-context.pdb": case / "prepared-full-context.pdb",
            "prepared-protein.pdb": case / "source.pdb",
            "request.json": case / "request.json",
            "lineage.json": case / "lineage.json",
            "native-model.pdb": Path(result["native"]),
            "native-geometry-analysis.json": directory / "result.json",
            "independent-readback.json": receipt_path,
        }
        original_request = Path(lineage["original_request"])
        if digest(original_request) != lineage["original_request_sha256"]:
            raise ValueError("original preparation request changed")
        preparation = original_request.parent
        inputs["source-preparation-result.json"] = preparation / "result.json"
        inputs["construct-target.fasta"] = preparation / "target.fasta"
        inputs["original-preparation-request.json"] = original_request
        if lineage["prepared_source_derivative"]:
            result_metadata = json.loads((preparation / "result.json").read_text())
            completion_root = Path(result_metadata["source_input"]).parent
            for path in sorted(completion_root.iterdir()):
                if path.is_file():
                    inputs[f"prior-input-completion/{path.name}"] = path
        candidate = inputs[f"{evidence['entry_id']}-restored-heavy.pdb"]
        stream = iter(candidate.read_text().splitlines(keepends=True))
        if not all(any(line == original for line in stream) for original in inputs["original-full-context.pdb"].read_text().splitlines(keepends=True)):
            raise ValueError("original deposited context records were not preserved")
        request = DiffusionRequest.from_json(inputs["request.json"].read_text())
        if digest(inputs["prepared-protein.pdb"]) != request.normalized_pdb.sha256:
            raise ValueError("frozen request source alias changed")
        dat = DatRecord.load(inputs[f"{evidence['entry_id']}-restored-heavy.dat"])
        owned = {(row["chain"], row["resid"], row["icode"], row["atom"]) for row in dat.added_atoms}
        generated = {(atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name) for atom in request.generated_atoms}
        if owned != generated or len(dat.added_atoms) != len(generated):
            raise ValueError("exported ownership differs from the exact completion mask")
        for label, row in evidence["CLI"].items():
            inputs[f"diagnose-{label}.json"] = Path(row["report"])
        for name, source in inputs.items():
            (destination / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination / name)
            if digest(destination / name) != digest(source):
                raise ValueError("export copy digest mismatch")
        manifest = {
            "kind": "private-experimental-restoration-review", "entry_id": evidence["entry_id"],
            "backend": "MODELLER 10.8 CPU", "native_seed": result["seed"], "stage": result["stage"],
            "heavy_coordinate_contract_passed": True, "fixed_full_context_atoms": evidence["fixed_atoms"],
            "fixed_reference_is_prepared_derivative": lineage["prepared_source_derivative"],
            "prior_completed_fixed_conditioning_atoms": sum(line.startswith(("ATOM  ", "HETATM")) for line in inputs["prepared-full-context.pdb"].read_text().splitlines())
                - sum(line.startswith(("ATOM  ", "HETATM")) for line in inputs["original-full-context.pdb"].read_text().splitlines()),
            "completed_heavy_atoms": len(generated), "new_geometry_findings": 0, "global_all_l": True,
            "source_diagnostics_by_severity": evidence["CLI"]["prepared"]["by_severity"],
            "candidate_diagnostics": evidence["CLI"]["candidate"]["by_severity"],
            "terminal_preparation_findings": evidence["CLI"]["candidate"]["terminal_preparation_findings"],
            "no_atoms_outside_frozen_mask_added_to_export": True,
            "observations_and_all_retained_context_unchanged": True,
            "request_source_alias": {"frozen_path": request.normalized_pdb.path, "artifact": "prepared-protein.pdb",
                                     "sha256": request.normalized_pdb.sha256},
            "coordinate_accuracy_MD_readiness_or_broad_qualification_established": False,
            "manual_donor_approval": False, "public_diffusion_bundle_or_capability": False,
            "modeled_atoms_are_not_experimental_observations": True,
            "files": {name: digest(destination / name) for name in inputs},
        }
        # A review manifest is intentionally not a public bundle completion index.
        save(destination / "review-manifest.json", manifest)
        exports.append({"entry_id": evidence["entry_id"], "directory": str(destination),
                        "manifest_sha256": digest(destination / "review-manifest.json"), "completed_heavy_atoms": len(generated)})
    save(output / "index.json", {"exports": exports, "readback_sha256": digest(receipt_path),
                                 "model_calls": 0, "optimizer_calls": 0, "public_capability_changed": False})
    print(json.dumps({"exported": exports}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.analysis.resolve(), args.readback.resolve(), args.output.resolve())
