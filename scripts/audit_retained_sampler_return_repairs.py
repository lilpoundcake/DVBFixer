#!/usr/bin/env python3
"""Independent PDB-precision readback; raw and full-inference failures stay failures."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from scripts.audit_retained_loop_context_sampling import fixed_records
from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_1he8_traced_mobility import anchor_residues

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def check_observed_records(source, candidate, fixed, allowed):
    before, after = fixed_records(source, fixed), fixed_records(candidate, fixed)
    for atom, line in before.items():
        current = after[atom]
        if atom not in allowed and current != line:
            raise ValueError("nonmobile source observation changed")
        if atom in allowed and (current[:30] != line[:30] or current[54:] != line[54:]):
            raise ValueError("released observation metadata changed")


def audit(root, output):
    plan = json.loads((root / "plan.json").read_text())
    source_snapshots = []
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            snapshot = root / "refiner-entry-executed-snapshot.py"
            if Path(name).name != "refine_retained_sampler_returns.py" or not snapshot.is_file() or digest(snapshot) != expected:
                raise ValueError("repair declaration input/code changed")
            source_snapshots.append({"declared_path": name, "executed_snapshot": str(snapshot), "sha256": expected})
    rows, bindings = [], {str(root / "plan.json"): digest(root / "plan.json")}
    for job in plan["jobs"]:
        work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        execution_path = work / "execution.json"
        receipt = json.loads(execution_path.read_text())
        if not receipt["completed"] or receipt["returncode"] or receipt["guard_stop_reason"]:
            rows.append({"case": job["case_id"], "repair_execution_completed": False,
                         "resource_stop": receipt["guard_stop_reason"], "returncode": receipt["returncode"]})
            bindings[str(execution_path)] = digest(execution_path)
            continue
        request_path, source_path, raw_path = work / "request.json", work / "source.pdb", work / "raw.pdb"
        request = DiffusionRequest.from_json(request_path.read_text())
        if digest(source_path) != request.normalized_pdb.sha256:
            raise ValueError("repair source changed")
        source, raw = source_path.read_text(), raw_path.read_text()
        check_observed_records(source, raw, set(request.fixed_atoms), frozenset())
        raw_xyz, _ = _coordinates(raw)
        assess = assessment_adapter(request, source)
        summary = json.loads((work / "result.json").read_text())
        if summary != receipt["result"]:
            raise ValueError("repair result/receipt mismatch")
        verified = []
        prior_pass = False
        for index, declared in enumerate(summary["stages"]):
            stage = declared["stage"]
            expected_stage = ("generated-only", "generated-plus-sequence2")[index]
            if stage != expected_stage or prior_pass:
                raise ValueError("unnecessary/undeclared anchor-release stage")
            path = work / stage / "refined.pdb"
            if digest(path) != declared["final_sha256"]:
                raise ValueError("repair endpoint changed")
            text = path.read_text()
            selected = () if index == 0 else anchor_residues(request, source, raw, "sequence2")
            restraints = build_anchor_restraints(source, selected, observed_atoms=request.fixed_atoms,
                base_k_kj_mol_nm2=10000., bfactor_source="disabled") if selected else ()
            allowed = frozenset(restraint.atom for restraint in restraints)
            check_observed_records(source, text, set(request.fixed_atoms), allowed)
            coordinates, _ = _coordinates(text)
            if set(coordinates) != set(raw_xyz):
                raise ValueError("repair atom ownership changed")
            assessment = assess(text, allowed)
            movements = check_anchor_movements(restraints, coordinates, AnchorMovementLimits())
            structure = _structure(text, "repair.pdb")
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
            passed = assessment.geometry_passed and not assessment.admission_failures and all(m.within_limits for m in movements)
            if passed != declared["passed"] or asdict(assessment) != declared["assessment"]:
                # JSON serializes tuple containers as lists.
                if passed != declared["passed"] or json.loads(json.dumps(asdict(assessment))) != declared["assessment"]:
                    raise ValueError("independent repair assessment disagrees")
            verified.append({"stage": stage, "passed": passed, "assessment": asdict(assessment),
                "movements": [asdict(m) for m in movements], "global_all_l": True,
                "nonmobile_observed_records_exact": True, "observed_metadata_exact": True,
                "generated_atoms": len(request.generated_atoms), "observed_atoms_released": len(allowed),
                "final_path": str(path), "final_sha256": digest(path)})
            prior_pass = passed
            bindings[str(path)] = digest(path)
        accepted = bool(verified and verified[-1]["passed"])
        if accepted != summary["accepted"]:
            raise ValueError("repair acceptance mismatch")
        rows.append({"case": job["case_id"], "repair_execution_completed": True, "accepted": accepted,
            "native_backend": summary.get("native_backend", "Protenix"),
            "native_seed": summary.get("native_seed", request.seeds[0]),
            "stages": verified, "raw_assessment": asdict(assess(raw)), "raw_failure_not_reclassified": True,
            "native_execution_status_not_reclassified": True, "new_sampling_calls": 0,
            "solver_precision_coordinates_separately_saved": False, "sparse_trajectory_or_topology_qualified": False,
            "not_manual_donor_approval_or_publication": True})
        for p in (execution_path, request_path, source_path, raw_path, work / "result.json"):
            bindings[str(p)] = digest(p)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "independent-readback.json", {"rows": rows, "bindings": bindings,
        "declaration_source_snapshots": source_snapshots,
        "auditor_sha256": digest(Path(__file__)), "readback_precision": "saved output PDB, not unsaved solver state"})
    print(json.dumps({"cases": len(rows), "repair_accepted": sum(row.get("accepted", False) for row in rows)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
