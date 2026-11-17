#!/usr/bin/env python3
"""Reassess saved full-context native stages; zero model/optimizer invocations."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_modeller_full_context import (
    compare_findings,
    findings,
    fixed_audit,
    merge_projection,
    native_transport,
    ownership,
    read_atoms,
)

from dvbfixer.model.diffusion.contract import DiffusionRequest


def declare(previous, output):
    plan_path = previous / "plan.json"
    old = json.loads(plan_path.read_text())
    snapshots = previous / "executed-code/index.json"
    index = json.loads(snapshots.read_text())
    if digest(plan_path) != index["plan_sha256"]:
        raise ValueError("historical plan snapshot changed")
    aliases = {row["original_path"]: Path(row["snapshot_path"]) for row in index["snapshots"]}
    files = [plan_path, snapshots, Path(__file__), Path(__file__).with_name("retained_modeller_full_context.py")]
    for name, expected in old["source_code_and_input_sha256"].items():
        path = Path(name)
        if digest(path) != expected:
            if path.suffix != ".py" or name not in aliases or digest(aliases[name]) != expected:
                raise ValueError("historical noncode input or code snapshot changed")
            path = aliases[name]
        files.append(path)
    jobs = []
    for job in old["jobs"]:
        work = previous / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        files += [work / "execution.json", work / "native/native-interface-calls.json", work / "native/native-outputs.json"]
        for label, filename in (("raw", "native.pdb"), ("generated-only-CG", "generated-only-CG-native.pdb")):
            native = work / "native" / filename
            if native.is_file():
                jobs.append({"entry_id": job["case_id"], "seed": job["seed"], "stage": label,
                             "case": job["case"], "native": str(native), "original_execution": str(work / "execution.json")})
                files.append(native)
    if [(j["entry_id"], j["stage"]) for j in jobs] != [("1CC6", "raw"), ("1CC6", "generated-only-CG"), ("9M3K", "raw")]:
        raise ValueError("recovery binds only the three retained native stages")
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", {"previous": str(previous), "jobs": jobs,
        "bindings": {str(p): digest(p) for p in files}, "model_calls": 0, "optimizer_calls": 0,
        "repair": "serial-only identity alias resolution; apply record-custody check after native-coordinate audit",
        "source_coordinates_or_frozen_masks_changed": False, "historical_results_overwritten": False})


def run(output):
    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"saved-analysis input changed: {name}")
    reports = []
    for job in plan["jobs"]:
        case, native = Path(job["case"]), Path(job["native"])
        directory = output / job["entry_id"] / job["stage"]
        directory.mkdir(parents=True, exist_ok=False)
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source = (case / "source.pdb").read_text()
        context = (case / "prepared-full-context.pdb").read_text()
        projection, transported, mapping = native_transport(request, context, native.read_text())
        audit = fixed_audit(context, transported)
        save(directory / "native-fixation.json", audit)
        save(directory / "mapping.json", mapping)
        atoms = read_atoms(projection)
        identities = tuple(request.generated_atoms)
        restored = materialize_candidate_pdb(source, identities, [atoms[key].xyz for key in identities],
            tuple(atoms[key].name for key in identities), request, elements=tuple(atoms[key].element for key in identities))
        merged = merge_projection(restored, context)
        (directory / "projection.pdb").write_text(restored)
        (directory / "full-context.pdb").write_text(merged)
        ownership(request, projection, directory / "generated.dat")
        assessment = assessment_adapter(request, source)(restored)
        structure = _structure(projection, "saved-native.pdb")
        all_l = True
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except ChiralityError:
            all_l = False
        baseline = findings(context, directory / "baseline.pdb")
        full = compare_findings(baseline, findings(merged, directory / "diagnostic-input.pdb"))
        raw_job = next(row for row in plan["jobs"] if row["entry_id"] == job["entry_id"] and row["stage"] == "raw")
        raw_projection, _, _ = native_transport(request, context, Path(raw_job["native"]).read_text())
        raw_atoms = read_atoms(raw_projection)
        max_motion = max(math.dist(atoms[key].xyz, raw_atoms[key].xyz) for key in identities)
        accepted = (audit["exact_at_PDB_precision"] and all_l and assessment.geometry_passed
                    and not assessment.admission_failures and not full["new_or_changed_findings"] and max_motion <= 3.)
        result = {**job, "directory": str(directory), "native_fixed": audit, "global_all_l": all_l,
                  "assessment": asdict(assessment), "full_diagnostics": full, "accepted": accepted,
                  "max_generated_motion_from_raw_angstrom": max_motion, "movement_limit_angstrom": 3.,
                  "source_findings": [asdict(row) for row in baseline],
                  "model_calls": 0, "optimizer_calls": 0, "public_output_published": False}
        save(directory / "result.json", result)
        reports.append(result)
    save(output / "completed.json", {"rows": reports, "native_stages": len(reports), "model_calls": 0, "optimizer_calls": 0,
                                     "historical_outcomes_reclassified": False})
    print(json.dumps({"rows": [{"entry": r["entry_id"], "stage": r["stage"], "fixed_changes": len(r["native_fixed"]["changed_atoms"]),
                               "geometry_issues": len(r["assessment"]["issues"]), "new_findings": len(r["full_diagnostics"]["new_or_changed_findings"]),
                               "accepted": r["accepted"]} for r in reports]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.previous is None:
            parser.error("declare requires --previous")
        declare(args.previous.resolve(), args.output.resolve())
    else:
        run(args.output.resolve())
