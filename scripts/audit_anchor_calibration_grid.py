#!/usr/bin/env python3
"""Read-only reanalysis: protein-request view AND independently retained context.

The original grid's protein-only assessment received full-context files. Keep
those reports unchanged, audit all context separately, and construct a coordinate-
identical request-scoped evaluation view. No refinement or sampling occurs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_preflight import independent_windows
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import verify_bindings

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def key(a):
    return a.chain, a.residue_number, a.insertion_code, a.atom_name


def evaluation_view(request, source, full_path, anchors):
    """Every evaluated coordinate comes from final bytes; never restore motion."""
    full = atoms(full_path)
    axis = tuple(request.generated_atoms)
    view = materialize_candidate_pdb(source, axis, [full[key(a)]["xyz"] for a in axis],
        tuple(full[key(a)]["resname"] for a in axis), request,
        elements=tuple(full[key(a)]["element"] for a in axis))
    view = rewrite_generated_coordinates(view, {a: np.asarray(full[key(a)]["xyz"]) for a in anchors})
    return view


def record_audit(before, after, mobile):
    old, new = before.splitlines(keepends=True), after.splitlines(keepends=True)
    if len(old) != len(new):
        raise ValueError("full-context physical record count changed")
    seen = set()
    for first, second in zip(old, new, strict=True):
        atom = None
        if first.startswith(("ATOM  ", "HETATM")):
            atom = (first[21], first[22:26].strip(), first[26].strip(), first[12:16].strip())
        if atom in mobile:
            if first[:30] != second[:30] or first[54:] != second[54:]:
                raise ValueError("mobile metadata changed")
            seen.add(atom)
        elif first != second:
            raise ValueError("outside-window record/ANISOU/annotation changed")
    if seen != mobile:
        raise ValueError("mobile inventory differs")


def movement(rows, reference, points, envelope):
    grouped = {}
    for atom in rows:
        grouped.setdefault(atom[:3], {})[atom[3]] = math.dist(reference[atom]["xyz"], points[atom])
    measured = []
    for residue, values in sorted(grouped.items()):
        rmsd = math.sqrt(sum(values[n]**2 for n in ("N", "CA", "C", "O"))/4)
        maximum, ca = max(values.values()), values["CA"]
        measured.append({"residue": residue, "backbone_rmsd": rmsd, "ca": ca, "heavy_max": maximum,
                         "within_limits": rmsd <= envelope[0] and ca <= envelope[1] and maximum <= envelope[2]})
    return measured


def diagnose(path, output):
    command = [sys.executable, "-m", "dvbfixer.cli", "diagnose", str(path), "--format", "json", "--output", str(output),
               "--include-water", "--clash-cutoff", "0.6,0.9"]
    env = {**os.environ, "PATH": str(Path(sys.executable).parent)+os.pathsep+os.environ.get("PATH", ""),
           "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
    with output.with_suffix(".log").open("xb") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env, timeout=180)
    if result.returncode not in (0, 1) or not output.is_file():
        raise RuntimeError("whole-context CLI diagnostics did not produce a report")
    return json.loads(output.read_text())["findings"]


def finding_difference(before, after):
    signatures = Counter(json.dumps(r, sort_keys=True, ensure_ascii=False) for r in before)
    new = []
    for row in after:
        signature = json.dumps(row, sort_keys=True, ensure_ascii=False)
        if signatures[signature]:
            signatures[signature] -= 1
        else:
            new.append(row)
    return new


def run(root, output):
    from scripts.anchor_calibration_preflight import prepared_context

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure

    output.mkdir(parents=True, exist_ok=False)
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    (output / "executed-audit-code.py").write_bytes(Path(__file__).read_bytes())
    retained = [plan_path]
    for job in plan["jobs"]:
        directory = root / f"execution/{job['case_id']}-seed-{job['seed']}"
        retained += [p for p in directory.rglob("*") if p.is_file() and p.suffix in {".pdb", ".dat", ".npz", ".json", ".xml"}]
    save(output / "declaration.json", {"audit_code_sha256": digest(Path(__file__)),
        "grid_plan_sha256": digest(plan_path), "retained_artifact_bindings": {str(p): digest(p) for p in retained},
        "optimizer_calls": 0, "original_reports_will_not_be_overwritten": True,
        "whole_context_cli": "include-water; unchanged0.6/0.9A steric thresholds",
        "scope": "correct protein-adapter evaluation view plus separate full-context/native/metadata audit"})
    baseline_cache = {}
    results = []
    for job in plan["jobs"]:
        case = Path(job["case"])
        work = root / f"execution/{job['case_id']}-seed-{job['seed']}"
        original = json.loads((work / "result.json").read_text())
        destination = output / job["case_id"]
        destination.mkdir()
        row = {"entry": case.name, "cell": original["cell"], "original_result_sha256": digest(work / "result.json"),
               "minimizer_calls": 0, "native_failures_reclassified_as_passes": False}
        if not original["completed"]:
            from openmm.unit import nanometer
            from scripts.anchor_calibration_grid import read_temporary

            from dvbfixer.ffutils.geometry import find_d_residues

            snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
            states = sorted(work.glob("phase-*/returned-state.npz"))
            if not states:
                raise ValueError("reported failure has no retained returned state")
            topology, _, _ = read_temporary(snapshot)
            with np.load(states[-1], allow_pickle=False) as data:
                failed_points = data["positions_nm"].copy()
            offenders = find_d_residues(topology, failed_points*nanometer)
            if "global_chirality_failure" not in original["outcome_flags"] or not offenders:
                raise ValueError("retained failure does not independently reproduce native D-Cα rejection")
            row.update(final_candidate_available=False, rejection=original["outcome_flags"], reason=original["reason"])
            row.update(independent_D_ca_rejection=True, D_ca_findings=offenders, rejected_state_sha256=digest(states[-1]))
            save(destination / "result.json", row)
            results.append(row)
            continue
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source = (case / "source.pdb").read_text()
        initial = atoms(case / "proposal.pdb")
        final_path = work / "refined.pdb"
        if digest(final_path) != original["final_sha256"]:
            raise ValueError("retained final bytes changed")
        final = atoms(final_path)
        source_atoms = atoms(case / "source.pdb")
        anchor_residues = {(r.chain, r.residue_number, r.insertion_code)
                           for r in independent_windows(request, original["cell"]["window"])}
        anchors = {a for a in request.fixed_atoms if key(a)[:3] in anchor_residues}
        generated = set(request.generated_atoms)
        mobile = {key(a) for a in generated | anchors}
        if set(initial) != set(final):
            raise ValueError("full-context final inventory changed")
        record_audit((case / "proposal.pdb").read_text(), final_path.read_text(), mobile)
        snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
        axis = json.loads((snapshot / "temporary-axis.json").read_text())["atoms"]
        with np.load(work / "phase-4/returned-state.npz", allow_pickle=False) as data:
            points = data["positions_nm"]*10
        native = {tuple(a["identity"][n] for n in ("chain", "residue_number", "insertion_code", "atom_name")): points[i]
                  for i, a in enumerate(axis)}
        for atom, record in initial.items():
            if atom not in mobile and not np.allclose(native[atom], record["xyz"], rtol=0, atol=1.e-4):
                raise ValueError("native full-context fixation failed")
            if atom in mobile and any(f"{v:8.3f}" != final[atom]["line"][start:start+8]
                                      for v, start in zip(native[atom], (30, 38, 46), strict=True)):
                raise ValueError("serialized mobile coordinates do not match native state")
        view = evaluation_view(request, source, final_path, anchors)
        view_path = destination / "protein-evaluation-view.pdb"
        view_path.write_text(view)
        projection = atoms(view_path)
        if set(projection) != {key(a) for a in generated | set(request.fixed_atoms)}:
            raise ValueError("request-scoped evaluation inventory differs")
        if any(projection[a]["xyz"] != final[a]["xyz"] for a in projection):
            raise ValueError("evaluation view altered actual final coordinates")
        parsed = _structure(final_path.read_text(), "independent-full-context.pdb")
        assert_all_l(parsed.pdb.topology, parsed.pdb.positions)
        if case.name == "1KTZ":
            from dvbfixer.model.diffusion.rfd3_custody_records import read_records

            old_pairs = read_records((case / "proposal.pdb").read_text())[2]
            new_pairs = read_records(final_path.read_text())[2]
            if new_pairs != old_pairs or len(new_pairs) != 21:
                raise ValueError("approved21SS annotation target changed")
        dat = DatRecord.load(work / "generated.dat")
        if digest(work / "generated.dat") != original["generated_dat_sha256"]:
            raise ValueError("retained generated dat bytes changed")
        owned = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in dat.added_atoms}
        if owned != {key(a) for a in generated} or len(owned) != len(dat.added_atoms):
            raise ValueError("generated-only dat ownership changed")
        native_moves = movement({key(a) for a in anchors}, source_atoms, native, original["cell"]["envelope"])
        serialized_moves = movement({key(a) for a in anchors}, source_atoms,
                                   {a: r["xyz"] for a, r in final.items()}, original["cell"]["envelope"])
        bounded = all(r["within_limits"] for r in [*native_moves, *serialized_moves])
        assessed = assessment_adapter(request, source)(view, anchors)
        issues = []
        internal = {r for gap in request.gaps if gap.left_anchor is not None and gap.right_anchor is not None for r in gap.generated_residues}
        terminal = {r for gap in request.gaps if gap.left_anchor is None or gap.right_anchor is None for r in gap.generated_residues}
        for issue in assessed.issues:
            issues.append({**asdict(issue), "touches_internal_generated": any(r in internal for r in issue.residues),
                           "touches_terminal_generated": any(r in terminal for r in issue.residues)})
        if case.name not in baseline_cache:
            baseline = (case / "context.pdb").read_text()
            if case.name != "1KTZ":
                baseline, _ = prepared_context(source, baseline)
            baseline_path = output / f"{case.name}-prepared-baseline.pdb"
            baseline_path.write_text(baseline)
            baseline_cache[case.name] = diagnose(baseline_path, output / f"{case.name}-baseline-cli.json")
        cli = diagnose(final_path, destination / "full-context-cli.json")
        introduced = finding_difference(baseline_cache[case.name], cli)
        row.update(final_candidate_available=True, outside_window_records_exact=True,
            native_full_context_fixed=True, global_all_l=True, generated_dat_exact=True,
            protein_view_coordinate_identical=True, protein_view_scope="frozen protein request; all context/approved extra links audited separately",
            protein_geometry_passed=assessed.geometry_passed, protein_admission_failures=list(assessed.admission_failures),
            issues=issues, native_anchor_movements=native_moves, serialized_anchor_movements=serialized_moves,
            movement_passed=bounded, restoration_mask_and_motion_passed=bounded and assessed.geometry_passed and not assessed.admission_failures,
            full_cli_counts=dict(Counter(r["severity"] for r in cli)), introduced_cli_findings=introduced,
            terminal_preparation_findings=[r for r in cli if r["category"] == "missing_terminals"],
            full_cli_no_new_findings=not introduced, original_final_sha256=digest(final_path))
        save(destination / "result.json", row)
        results.append(row)
        print(job["case_id"], "independently audited", flush=True)
    summary = []
    for entry in plan["admitted"]:
        selected = [r for r in results if r["entry"] == entry]
        summary.append({"entry": entry, "reported": len(selected),
            "retained_final_candidates": sum(r["final_candidate_available"] for r in selected),
            "mask_and_motion_passes": sum(bool(r.get("restoration_mask_and_motion_passed")) for r in selected),
            "whole_cli_no_new_findings": sum(bool(r.get("full_cli_no_new_findings")) for r in selected),
            "passing_cells": [r["cell"]["id"] for r in selected if r.get("restoration_mask_and_motion_passed")],
            "original_worker_failure_flags": dict(Counter(f for r in selected for f in r.get("rejection", [])))})
    save(output / "receipt.json", {"grid_plan_sha256": digest(plan_path), "audit_code_sha256": digest(Path(__file__)),
        "cases": summary, "optimizer_calls": 0, "original_results_unchanged": True,
        "original_assessment_scope_error": "protein-only adapter received full-context candidates; original gate flags retained, not scientific success rates",
        "policy_selected": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
