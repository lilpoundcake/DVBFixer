"""One private deterministic additional minimizer call on retained RFD3-local.

Same original source anchors, atom mask, base score and geometry/movement gate.
An explicit junction-control mode adds one declared harmonic C-N term. Neither
mode packs, generates proposals, resets H state or re-anchors observations.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import key, launch_worker
from scripts.local_gap_rosetta_worker import checked_junction_restraint
from scripts.mac_reconstruction_comparison import digest, save
from scripts.one_non_local_research import assess, environment, movements
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def run(args):
    case, root = args.case.resolve(), args.output.resolve()
    junction_control = args.junction_control
    budget_control = getattr(args, "budget_control", False)
    polished_parent = junction_control or budget_control
    parent = case / ("polish-rfd3-local-1000-attempt01" if polished_parent else "local-retry/rfd3-local")
    receipt_path = parent / "independent-readback.json" if polished_parent else case / "independent-readback-v2.json"
    receipt = json.loads(receipt_path.read_text())
    result = json.loads((parent / "result.json").read_text())
    if (not result["native_movement_bounds_passed"] or not result["global_all_l"]
            or (receipt["candidate_sha256"] if polished_parent else receipt["cells"]["rfd3-local"]["candidate_sha256"])
            != digest(parent / "candidate.pdb")):
        raise ValueError("requires the independently audited bounded all-L parent")
    parent_plan = parent / "plan.json" if polished_parent else case / "local-retry/plan.json"
    verify_bindings(json.loads(parent_plan.read_text())["source_code_and_input_sha256"], (parent_plan,))
    root.mkdir(parents=True, exist_ok=False)
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    source = atoms(case / "source.pdb")
    payload = json.loads((parent / "payload.json").read_text())
    if {tuple(r["identity"]): r["xyz"] for r in payload["source"]} != {a: list(row["xyz"]) for a, row in source.items()}:
        raise ValueError("parent restraints are not the original source coordinates")
    payload.update(max_iterations=1000, resume_from=str(parent))
    if budget_control:
        payload.update(max_iterations=3000, budget_control={"passes": 3, "iterations_per_pass": 3000})
    if junction_control:
        payload["junction_length_restraint"] = {"first": ["B", "90", "", "C"], "second": ["B", "91", "", "N"],
                                                 "target_A": 1.33, "sd_A": .02, "weight": 1.}
        checked_junction_restraint(payload)
    save(root / "payload.json", payload)
    (root / "starting.pdb").write_bytes((parent / "starting.pdb").read_bytes())
    worker = Path("scripts/local_gap_rosetta_worker.py").resolve()
    files = [case / "source.pdb", case / "request.json", case / "local-mask.json", receipt_path,
             parent_plan, *[parent / n for n in ("payload.json", "starting.pdb", "native-before.json", "native-after.json", "native.json", "result.json", "candidate.pdb")],
             Path(__file__).resolve(), worker, Path("scripts/one_non_local_research.py").resolve(),
             Path("tests/test_one_non_local.py").resolve(),
             Path("scripts/retained_evidence.py").resolve(), root / "starting.pdb", root / "payload.json"]
    save(root / "plan.json", {"protocol": "private-1non-local-budget-control-v1" if budget_control else "private-1non-local-junction-control-v1" if junction_control else "private-1non-local-same-state-polish-v1", "parent": str(parent),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "case": str(case),
        "extra_minimizer_calls": 3 if budget_control else 1, "extra_max_iterations": 3000 if budget_control else 1000,
        "new_packing_calls": 0, "new_diffusion_calls": 0,
        "original_source_restraints_and_mobility_unchanged": True, "LBFGS_history_reinitialized": True,
        "all_native_coordinates_and_parent_base_score_verified_before_new_restraint": True,
        "junction_length_restraint": payload.get("junction_length_restraint"),
        "objective_change": "one new harmonic C-N atom-pair restraint; parent base score verified before addition" if junction_control else None,
        "gate_unchanged": True, "public_capability": False, "timeout_seconds": 1800 if budget_control else 900,
        "resources": "CPU1thread/pressure-v2/RSS8GiB; no solvent"})
    freeze_code(root / "plan.json")
    command = [str(args.python.absolute()), str(worker), str(root)]
    save(root / "command.json", command)
    execution = launch_worker(command, root, environment(), 1800 if budget_control else 900)
    verify_bindings(json.loads((root / "plan.json").read_text())["source_code_and_input_sha256"])
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("polish did not complete; raw failure retained without retry")
    native = json.loads((root / "native.json").read_text())
    coords = {tuple(row["identity"]): row["xyz"] for row in native["atoms"]}
    allowed = set(map(tuple, payload["observed_mobile"]))
    generated = set(map(key, request.generated_atoms))
    if len(coords) != len(native["atoms"]) or set(coords) != set(source) | generated:
        raise ValueError("polish native atom axis changed")
    error = max(math.dist(row["xyz"], coords[a]) for a, row in source.items() if a not in allowed)
    if error > 1e-7:
        raise ValueError("polish moved an observation outside the original window")
    updates = {a: np.asarray(coords[key(a)]) for a in (*request.fixed_atoms, *request.generated_atoms) if key(a) in allowed | generated}
    candidate = root / "candidate.pdb"
    candidate.write_text(rewrite_generated_coordinates((root / "starting.pdb").read_text(), updates))
    final = assess(request, case / "source.pdb", candidate, allowed)
    motion, outside = movements(source, {a: {"xyz": p} for a, p in coords.items()}, allowed)
    within = outside <= 1e-7 and all(r["within_limits"] for r in motion)
    final.update(execution_status="completed", native_local_observed_movements=motion,
                 native_outside_mask_max_A=error, native_movement_bounds_passed=within,
                 native_score=native["score"], parent_native_score=json.loads((parent / "native.json").read_text())["score"])
    final["native_base_score"] = native.get("base_score", native["score"])
    final["junction_C_N_A_native"] = native.get("junction_C_N_A")
    final["junction_length_restraint"] = payload.get("junction_length_restraint")
    final["coordinate_policy_passed"] = final["coordinate_policy_passed"] and within
    final["validation_passed"] = final["validation_passed"] and within
    final["quality_status"] = "accepted" if final["validation_passed"] else "rejected"
    save(root / "result.json", final)
    print(json.dumps({k: final[k] for k in ("execution_status", "quality_status", "geometry_passed", "coordinate_policy_passed", "full_geometry_error_counts", "global_all_l", "native_score", "parent_native_score")}))


def audit(root):
    """Read native states and final bytes independently; no optimization or rewriting."""
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.ffutils.geometry import assert_all_l

    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    case, parent = Path(plan["case"]), Path(plan["parent"])
    payload, old_payload = [json.loads((p / "payload.json").read_text()) for p in (root, parent)]
    restraint = checked_junction_restraint(payload)
    if ({k: v for k, v in payload.items() if k not in {"resume_from", "max_iterations", "junction_length_restraint", "budget_control"}} != {
            k: v for k, v in old_payload.items() if k not in {"max_iterations", "resume_from"}}
            or old_payload.get("junction_length_restraint") is not None):
        raise ValueError("polish changed the original source references or mobility")
    before = json.loads((root / "native-before.json").read_text())
    after = json.loads((root / "native-after.json").read_text())
    old_before = json.loads((parent / "native-before.json").read_text())
    old_after = json.loads((parent / "native-after.json").read_text())
    if (before["atoms"] != old_before["atoms"] or before["coordinates"] != old_after["coordinates"]
            or before["packing_executed"] or after["max_iterations"] != (3000 if payload.get("budget_control") else 1000)):
        raise ValueError("polish did not begin with the exact parent's complete native state")
    if not restraint and not math.isclose(before["score"], json.loads((parent / "native.json").read_text())["score"],
                                         rel_tol=1e-10, abs_tol=1e-6):
        raise ValueError("unaugmented resume changed the parent's objective")
    if restraint:
        proof = json.loads((root / "junction-length-restraint.json").read_text())
        old_score = json.loads((parent / "native.json").read_text())["score"]
        index = {tuple(row["identity"]): i for i, row in enumerate(before["atoms"])}
        for label, state, base, augmented in (("before", before["coordinates"], before["base_score"], before["score"]),
                                             ("after", after["coordinates"], json.loads((root / "native.json").read_text())["base_score"],
                                              json.loads((root / "native.json").read_text())["score"])):
            distance = math.dist(state[index[tuple(restraint["first"])]], state[index[tuple(restraint["second"])]])
            penalty = restraint["weight"] * ((distance-restraint["target_A"])/restraint["sd_A"])**2
            if not math.isclose(augmented-base, penalty, rel_tol=1e-8, abs_tol=1e-6):
                raise ValueError(f"{label} energy delta is not exactly the declared junction term")
        if not math.isclose(before["base_score"], old_score, rel_tol=1e-10, abs_tol=1e-6) or not proof["both_endpoints_mobile"]:
            raise ValueError("original objective or junction endpoint mobility changed")
    source = atoms(case / "source.pdb")
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    if digest(case / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("original source binding changed")
    fixed, generated = set(source), set(map(key, request.generated_atoms))
    allowed = set(map(tuple, payload["observed_mobile"]))
    native = json.loads((root / "native.json").read_text())
    coords = {tuple(r["identity"]): r["xyz"] for r in native["atoms"]}
    x, y = np.array(before["coordinates"]), np.array(after["coordinates"])
    if x.shape != y.shape or x.shape != (len(before["atoms"]), 3) or not np.isfinite([x, y]).all():
        raise ValueError("polish native particle state is nonfinite or incomplete")
    if len(coords) != len(native["atoms"]) or set(coords) != fixed | generated:
        raise ValueError("polish changed heavy-atom ownership")
    for r, p, q in zip(before["atoms"], x, y, strict=True):
        a = tuple(r["identity"])
        mobile = not r["virtual"] and (a in generated | allowed or r["hydrogen"])
        if r["mobile"] != mobile or (not mobile and np.linalg.norm(p-q) > 1e-7):
            raise ValueError("polish native fixation/mobility failed")
        if a in coords and not np.array_equal(q, coords[a]):
            raise ValueError("requested heavy coordinates do not match the native optimizer return")
    candidate = atoms(root / "candidate.pdb")
    if set(candidate) != fixed | generated:
        raise ValueError("polish serialized atom inventory changed")
    for a in candidate:
        if candidate[a]["xyz"] != tuple(float(f"{v:.3f}") for v in coords[a]):
            raise ValueError("polish serialization changed native coordinates")
        if a in fixed:
            first, second = source[a]["line"], candidate[a]["line"]
            if (a not in allowed and first != second) or first[:30] != second[:30] or first[54:] != second[54:]:
                raise ValueError("polish changed outside-window records or observed metadata")
    result = json.loads((root / "result.json").read_text())
    motion, outside = movements(source, {a: {"xyz": p} for a, p in coords.items()}, allowed)
    if motion != result["native_local_observed_movements"] or outside > 1e-7 or not all(r["within_limits"] for r in motion):
        raise ValueError("polish movement limits do not reproduce")
    assessment = assessment_adapter(request, (case / "source.pdb").read_text())((root / "candidate.pdb").read_text(),
        frozenset(a for a in request.fixed_atoms if key(a) in allowed))
    from dataclasses import asdict
    if json.loads(json.dumps(asdict(assessment))) != result["request_geometry"]:
        raise ValueError("polish request geometry decision does not reproduce")
    structure = PDBFile(StringIO((root / "candidate.pdb").read_text()))
    assert_all_l(structure.topology, structure.positions)
    clashes = clashes_python(structure.topology, structure.positions, include_water=True, clash_warn_a=.6, clash_error_a=.9)
    severe = sum(f.severity == Severity.ERROR for f in clashes)
    if (severe != result["full_geometry_error_counts"].get("clash", 0)
            or (result["validation_passed"] and (severe or not assessment.geometry_passed or assessment.admission_failures))):
        raise ValueError("polish steric result/rejected decision does not reproduce")
    return {"plan_sha256": digest(plan_path), "auditor_sha256": digest(Path(__file__)),
        "all_native_parent_coordinates_preserved_at_resume": True, "native_atoms": len(x),
        "original_source_references_and_mobile_mask_preserved": True, "native_outside_mask_max_A": outside,
        "global_all_l": True, "severe_clashes": severe, "native_movement_bounds_passed": True,
        "request_geometry_passed": assessment.geometry_passed, "quality_status": result["quality_status"],
        "junction_length_restraint": restraint, "new_restraint_energy_delta_verified": bool(restraint),
        "candidate_sha256": digest(root / "candidate.pdb"), "new_minimizer_calls_in_audit": 0}


def audit_budget_passes(root):
    """Validate each retained intermediate against the same original source."""
    payload = json.loads((root / "payload.json").read_text())
    if payload.get("budget_control") != {"passes": 3, "iterations_per_pass": 3000} or payload.get("junction_length_restraint"):
        raise ValueError("requires the declared unaugmented budget control")
    plan = json.loads((root / "plan.json").read_text())
    case = Path(plan["case"])
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    source = atoms(case / "source.pdb")
    allowed = set(map(tuple, payload["observed_mobile"]))
    generated = set(map(key, request.generated_atoms))
    before = json.loads((root / "native-before.json").read_text())
    previous = before["coordinates"]
    previous_score = before["score"]
    reports = []
    for stage in range(1, 4):
        state = json.loads((root / f"pass-{stage}.json").read_text())
        current = state["coordinates"]
        if (len(current) != len(before["atoms"]) or not np.isfinite(current).all()
                or state["score_before"] != previous_score):
            raise ValueError("incomplete/nonfinite budget state or broken score chain")
        coords = {}
        for row, initial, point in zip(before["atoms"], before["coordinates"], current, strict=True):
            if not row["mobile"] and initial != point:
                raise ValueError("budget intermediate moved a fixed native particle")
            identity = tuple(row["identity"])
            if identity in source or identity in generated:
                coords[identity] = point
        step = max(math.dist(a, b) for a, b in zip(previous, current, strict=True))
        distance = math.dist(coords['B', '90', '', 'C'], coords['B', '91', '', 'N'])
        if step != state["maximum_atom_step_A"] or distance != state["junction_C_N_A"]:
            raise ValueError("budget progress measurements do not reproduce")
        updates = {a: np.asarray(coords[key(a)]) for a in (*request.fixed_atoms, *request.generated_atoms)
                   if key(a) in allowed | generated}
        candidate = root / f"pass-{stage}-candidate.pdb"
        with candidate.open("x") as stream:
            stream.write(rewrite_generated_coordinates((root / "starting.pdb").read_text(), updates))
        report = assess(request, case / "source.pdb", candidate, allowed)
        motion, outside = movements(source, {a: {"xyz": p} for a, p in coords.items()}, allowed)
        report.update(native_local_observed_movements=motion, native_outside_mask_max_A=outside,
                      native_movement_bounds_passed=outside <= 1e-7 and all(row["within_limits"] for row in motion))
        save(root / f"pass-{stage}-validation.json", report)
        reports.append({k: v for k, v in state.items() if k != "coordinates"} | {
            "geometry_passed": report["geometry_passed"], "quality_status": report["quality_status"],
            "native_movement_bounds_passed": report["native_movement_bounds_passed"],
            "full_geometry_error_counts": report["full_geometry_error_counts"],
            "global_all_l": report["global_all_l"]})
        previous, previous_score = current, state["score_after"]
    if previous != json.loads((root / "native-after.json").read_text())["coordinates"]:
        raise ValueError("last budget checkpoint differs from optimizer return")
    return reports


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path("/Users/lilpoundcake/micromamba/envs/pyrosetta/bin/python"))
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--junction-control", action="store_true",
                        help="Resume the polished parent once with the declared B90:C–B91:N length restraint")
    parser.add_argument("--budget-control", action="store_true",
                        help="Three unaugmented 3000-iteration passes from the pre-junction parent")
    args = parser.parse_args()
    if args.junction_control and args.budget_control:
        parser.error("budget and junction controls are mutually exclusive")
    if args.audit_only:
        result = audit(args.output.resolve())
        save(args.output / "independent-readback.json", result)
        print(json.dumps(result))
    else:
        run(args)
