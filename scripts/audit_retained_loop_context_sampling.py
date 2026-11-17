#!/usr/bin/env python3
"""Independent complete-endpoint readback for the two private development calls."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def fixed_records(text, fixed):
    rows = {}
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in fixed:
            if atom in rows:
                raise ValueError("duplicate source observation in endpoint")
            rows[atom] = line
    if set(rows) != fixed:
        raise ValueError("endpoint omitted original fixed observations")
    return rows


def audit(root, output):
    plan = json.loads((root / "plan.json").read_text())
    completed = json.loads((root / "completed.json").read_text())
    if completed["native_calls_completed"] != len(plan["jobs"]):
        raise ValueError("native technical cohort incomplete")
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"sampling declaration input changed: {name}")
    rows, bindings = [], {str(root / "plan.json"): digest(root / "plan.json")}
    for job in plan["jobs"]:
        work = Path(job["work"])
        execution = json.loads((work / "execution.json").read_text())
        if execution["returncode"] or execution["guard_stop_reason"]:
            raise ValueError("technical/resource failure is not an executed scientific outcome")
        device_path = work / "device-final.json"
        device = json.loads(device_path.read_text())
        if not device["model_mps_calls"] or device["automatic_cpu_fallback"]:
            raise ValueError("actual native model execution was not demonstrated")
        policy_path = work / "context-wiring-policy.json"
        policy = json.loads(policy_path.read_text())
        command = job["command"]
        case = Path(command[command.index("--case")+1])
        request_path = case / "context-request.json"
        request = DiffusionRequest.from_json(request_path.read_text())
        source_path = case / request.normalized_pdb.path
        if digest(request_path) != policy["request_sha256"] or digest(source_path) != request.normalized_pdb.sha256:
            raise ValueError("executed source/request differs from private policy")
        prediction = work / "predictions"
        candidate_path = prediction / "candidate.pdb"
        source, candidate = source_path.read_text(), candidate_path.read_text()
        if fixed_records(source, set(request.fixed_atoms)) != fixed_records(candidate, set(request.fixed_atoms)):
            raise ValueError("complete-assembly source observations changed")
        coordinates, _ = _coordinates(candidate)
        if set(coordinates) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("final atom ownership differs from declared complete endpoint")
        steps_path = prediction / "reinjection-steps.json"
        steps = json.loads(steps_path.read_text())
        if (len(steps) != 200 or [s["step_index"] for s in steps] != list(range(200))
                or any(s["step_count"] != 200 or s["post_overwrite_max_angstrom"] != 0 for s in steps)):
            raise ValueError("per-update fixed-coordinate trace incomplete or nonzero")
        assessment = assessment_adapter(request, source)(candidate)
        global_l, chirality_failure = True, None
        structure = _structure(candidate, "raw-private-development.pdb")
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except ChiralityError as exc:
            global_l, chirality_failure = False, str(exc)
        rows.append({"entry_id": job["case_id"], "raw_geometry_passed": assessment.geometry_passed,
            "raw_admission_failures": list(assessment.admission_failures), "assessment": asdict(assessment),
            "global_all_l": global_l, "chirality_failure": chirality_failure,
            "complete_source_fixed_records_exact": True, "generated_atoms": len(request.generated_atoms),
            "represented_fixed_atoms": policy["represented_fixed_atoms"], "preserved_only_atoms": policy["preserved_only_atoms"],
            "native_updates": len(steps), "actual_MPS_model_calls": device["model_mps_calls"],
            "raw_endpoint_passed": assessment.geometry_passed and not assessment.admission_failures and global_l,
            "candidate_path": str(candidate_path), "candidate_sha256": digest(candidate_path),
            "crystal_reference_used": False, "repair_calls": 0,
            "historical_native_atom_identity_axis_separately_dumped": False,
            "no_continuous_topology_or_public_qualification": True})
        for path in (work / "execution.json", device_path, policy_path, request_path, source_path, candidate_path, steps_path):
            bindings[str(path)] = digest(path)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "raw-endpoint-readback.json", {"rows": rows, "bindings": bindings, "auditor_sha256": digest(Path(__file__)),
        "native_axis_trace_limit": "Endpoint identities/custody and native input preflight are checked; full historical sampled axis was not separately serialized.",
        "technical_development_controls_not_final_evaluation": True})
    print(json.dumps({"readback_cases": len(rows), "raw_endpoint_passes": sum(r["raw_endpoint_passed"] for r in rows),
                      "all_source_records_exact": all(r["complete_source_fixed_records_exact"] for r in rows)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
