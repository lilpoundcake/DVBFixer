#!/usr/bin/env python3
"""Unapproved private review copies, never public bundles or donor selection."""
from __future__ import annotations

import argparse
import json
import shutil
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates


def generated_record(request, text):
    _, records = _coordinates(text)
    added = []
    for atom in request.generated_atoms:
        line = records[atom]
        added.append({"chain": atom.chain, "resid": atom.residue_number, "icode": atom.insertion_code,
            "resname": line[17:20].strip(), "atom": atom.atom_name,
            "element": line[76:78].strip() or atom.atom_name[0]})
    if {entry["chain"] for entry in added} - {atom.chain for atom in request.generated_atoms}:
        raise ValueError("generated sidecar changed chain identity")
    return DatRecord(description="PRIVATE proposal review: generated ownership only; not donor approval", added_atoms=added)


def native_provenance(row, request, optimizer_seed):
    backend = row.get("native_backend", "Protenix")
    if backend not in {"Protenix", "MODELLER"}:
        raise ValueError("unknown private review backend")
    seed = row.get("native_seed", request.seeds[0])
    if seed != request.seeds[0]:
        raise ValueError("review native seed differs from its bound request")
    return {"native_backend": backend, "native_seed": seed, "optimizer_seed": optimizer_seed,
        "native_full_inference_status": ("resource-refused; retained sampler return only" if backend == "Protenix"
                                         else "native MODELLER proposal completed; source restoration assessed separately")}


def export(repair_audit, repair_root, output):
    evidence = json.loads(repair_audit.read_text())
    for name, expected in evidence["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("independent repair evidence changed")
    plan_path = repair_root / "plan.json"
    if evidence["bindings"].get(str(plan_path)) != digest(plan_path):
        raise ValueError("review root differs from the independently audited repair plan")
    output.mkdir(parents=True, exist_ok=False)
    with (output / "exporter-source.py").open("xb") as stream:
        stream.write(Path(__file__).read_bytes())
    repair_plan = json.loads((repair_root / "plan.json").read_text())
    jobs = {job["case_id"]: job for job in repair_plan["jobs"]}
    rows = []
    for row in evidence["rows"]:
        case = row["case"]
        if case not in {"1NON", "2BE1"}:
            raise ValueError("undeclared private review case")
        optimizer_seed = jobs[case]["seed"]
        work = repair_root / "execution" / f"{case}-seed-{optimizer_seed}"
        directory = output / case
        directory.mkdir(exist_ok=False)
        request = DiffusionRequest.from_json((work / "request.json").read_text())
        provenance = native_provenance(row, request, optimizer_seed)
        provenance["native_protocol"] = repair_plan.get("native_protocol", "baseline")
        files, copies = {}, [(work / "source.pdb", "01-observed-source.pdb", "prepared-protein-only-source", None),
            (work / "raw.pdb", "02-raw-geometry-failed.pdb", "source-restored-proposal", False)]
        if row.get("repair_execution_completed"):
            for index, stage in enumerate(row["stages"]):
                copies.append((Path(stage["final_path"]), f"{index+3:02d}-{stage['stage']}.pdb", stage["stage"], stage["passed"]))
        metadata = []
        commands = []
        for source, name, stage, passed in copies:
            expected = evidence["bindings"].get(str(source))
            if expected is None or not source.resolve().is_relative_to(repair_root):
                raise ValueError("review source is not owned/bound by this repair audit")
            destination = directory / name
            shutil.copyfile(source, destination)
            if expected != digest(source) or expected != digest(destination):
                raise ValueError("review copy changed")
            files[name] = digest(destination)
            if passed is not None:
                sidecar = destination.with_suffix(".dat")
                generated_record(request, destination.read_text()).save(sidecar)
                checked = DatRecord.load(sidecar)
                identities = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in checked.added_atoms}
                expected = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
                if identities != expected or len(checked.added_atoms) != len(expected):
                    raise ValueError("review sidecar included observations or lost generated atoms")
                files[sidecar.name] = digest(sidecar)
            object_name = f"{provenance['native_backend'].lower()}_s{provenance['native_seed']}_{stage.replace('-', '_')}"
            commands.append(f"load {json.dumps(str(destination.resolve()), ensure_ascii=False)}, {object_name}")
            metadata.append({"file": name, "stage": stage, "geometry_passed": passed,
                "sha256": files[name], "source_path": str(source), "manual_approval": None})
        (directory / "overlay.pml").write_text("\n".join(commands)+"\n", encoding="utf-8")
        files["overlay.pml"] = digest(directory / "overlay.pml")
        shutil.copyfile(work / "request.json", directory / "request.json")
        if digest(directory / "request.json") != evidence["bindings"].get(str(work / "request.json")):
            raise ValueError("review request copy differs from audit")
        files["request.json"] = digest(directory / "request.json")
        recipe = {"entry_id": case, **provenance, "joint_group": "all-requested-gaps-one-native-proposal",
            "gaps": [asdict(gap) for gap in request.gaps], "stages": metadata,
            "selected_donor": None, "allowed_donors": [], "manual_approval": None,
            "geometry_acceptance_is_not_manual_approval": True,
            "unknown_reference_coordinates_not_available_or_used": True,
            "all_structures_already_in_source_frame_no_alignment": True,
            "scope": "prepared protein-only projection; not completion of omitted chemistry",
            "source_request_sha256": digest(work / "request.json"), "files": files}
        save(directory / "review.json", recipe)
        rows.append({"entry_id": case, **provenance, "review_path": f"{case}/review.json",
                     "review_sha256": digest(directory / "review.json"), "manual_approval": None})
    (output / "README.md").write_text(
        "# Private natural-gap review — not approved donors\n\n"
        "Each folder contains the actual prepared protein-only source, failed raw\n"
        "source-restored proposal, available repaired stages and explicit pass/fail metadata.\n"
        "There is no complete hidden reference. Native backend/seed/status and optimizer\n"
        "seed are distinct fields; Protenix full-inference resource refusals remain visible.\n"
        "Open overlay.pml manually in PyMOL; no fitting is applied. The .dat sidecars\n"
        "contain only whole-gap generated atoms, never fixed or released observations.\n"
        "No candidate is selected, approved, published or registered in a public GUI.\n",
        encoding="utf-8")
    save(output / "review-index.json", {"kind": "private-unapproved-review-v1", "cases": rows,
        "repair_audit_path": str(repair_audit), "repair_audit_sha256": digest(repair_audit),
        "exporter_sha256": digest(Path(__file__)), "exporter_source_snapshot": "exporter-source.py",
        "public_bundle": False, "all_manual_approvals_unset": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repair-audit", type=Path, required=True)
    parser.add_argument("--repair-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.repair_audit.resolve(), args.repair_root.resolve(), args.output.resolve())
