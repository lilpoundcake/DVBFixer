#!/usr/bin/env python3
"""Read-only v2 steric revalidation of the frozen 100-cell anchor study.

No optimization, source repair, new H placement or proposal generation. Full CLI
reports retain the actual selected engine; deterministic Python v2 is additionally
run on every source/proposal/final, and on the six rejected native FF states.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_grid import key
from scripts.audit_anchor_calibration_grid import run as audit_grid
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.diagnose import steric
from dvbfixer.diagnose.report import findings_to_dict_list
from dvbfixer.model.diffusion.contract import DiffusionRequest


def finding_key(row):
    """Compare contact identities without engine-version metadata churn."""
    details = row.get("details", row.get("extra", {}))
    identities = details.get("atoms", [])
    if row["category"] == "clash" and len(identities) == 2:
        return ("clash", row["severity"], tuple(sorted(
            (a["chain"], a["resid"], a.get("icode", ""), a["atom"]) for a in identities)))
    return (row["category"], row["severity"], row.get("chain"), row.get("resid"), row.get("atom"), row["message"])


def difference(before, after):
    counts = Counter(finding_key(row) for row in before)
    result = []
    for row in after:
        signature = finding_key(row)
        if counts[signature]:
            counts[signature] -= 1
        else:
            result.append(row)
    return result


def classify(row, generated, anchors):
    details = row.get("details", row.get("extra", {}))
    axis = {(a["chain"], a["resid"], a.get("icode", ""), a["atom"]) for a in details.get("atoms", [])}
    return {"finding": row, "touches_generated": bool(axis & generated),
            "touches_observed_anchor": bool(axis & anchors),
            "only_immutable_atoms": bool(axis) and not bool(axis & (generated | anchors))}


def python_report(topology, positions, output):
    provenance = {"actual_engine": steric.DETERMINISTIC_STERIC_ENGINE,
                  "engine_version": steric.DETERMINISTIC_STERIC_ENGINE_VERSION,
                  "include_water": True, "warning_threshold_angstrom": .6, "error_threshold_angstrom": .9}
    findings = steric.clashes_python(topology, positions, include_water=True,
                                    clash_warn_a=.6, clash_error_a=.9, provenance=provenance)
    rows = findings_to_dict_list(findings)
    result = {"steric": provenance, "findings": rows}
    save(output, result)
    return result


def run(grid, previous, forces, output):
    from openmm.app import PDBFile
    from openmm.unit import nanometer
    from scripts.anchor_calibration_grid import read_temporary

    if steric.DETERMINISTIC_STERIC_ENGINE_VERSION != "2":
        raise ValueError("this declaration requires deterministic steric engine v2")
    output.mkdir(parents=True, exist_ok=False)
    plan = json.loads((grid / "plan.json").read_text())
    bindings = {}
    for root in (Path("src/dvbfixer"),):
        for path in root.rglob("*.py"):
            bindings[str(path.resolve())] = digest(path)
    for name in ("reaudit_anchor_diagnostics.py", "audit_anchor_calibration_grid.py", "retained_evidence.py"):
        path = Path(__file__).with_name(name)
        bindings[str(path.resolve())] = digest(path)
    originals = [p for p in previous.rglob("*") if p.is_file()]
    originals += [forces]
    historical = {str(p): digest(p) for p in originals}
    save(output / "declaration.json", {"revision": "anchor-steric-v2-readback-v1",
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "code_sha256": bindings, "historical_report_bindings": historical,
        "grid_plan_sha256": digest(grid / "plan.json"), "optimizer_calls": 0,
        "hydrogen_placement_calls": 0, "source_repairs": 0, "original_reports_overwritten": False,
        "thresholds_angstrom": [.6, .9], "scope": "100 cells; 94 serialized candidates + 6 rejected native FF states; original/prepared sources and proposals"})
    (output / "executed-analysis.py").write_bytes(Path(__file__).read_bytes())
    # The existing audit independently checks bytes, native fixation, masks,
    # movement, disulfides and rejected D states; it never optimizes anything.
    readback = output / "readback"
    audit_grid(grid, readback)
    force_rows = {row["cell_id"]: row for row in json.loads(forces.read_text())["rows"]}
    source_cache, rows = {}, []
    for job in plan["jobs"]:
        case, cell = Path(job["case"]), job["case_id"]
        work = grid / f"execution/{cell}-seed-{job['seed']}"
        original = json.loads((work / "result.json").read_text())
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        generated = {key(a) for a in request.generated_atoms}
        fresh = json.loads((readback / cell / "result.json").read_text())
        anchors = {tuple(r["residue"]) for r in fresh.get("native_anchor_movements", [])}
        anchor_atoms = {a for a in atoms(case / "proposal.pdb") if a[:3] in anchors}
        if case.name not in source_cache:
            directory = output / case.name
            directory.mkdir()
            reports = {}
            for label, path in (("original-context", case / "context.pdb"),
                                ("prepared-baseline", readback / f"{case.name}-prepared-baseline.pdb"),
                                ("proposal", case / "proposal.pdb")):
                pdb = PDBFile(str(path))
                report = python_report(pdb.topology, pdb.positions, directory / f"{label}-python-v2.json")
                reports[label] = report
                reports[label]["input_sha256"] = digest(path)
            source_cache[case.name] = reports
        destination = readback / cell
        if original["completed"]:
            pdb = PDBFile(str(work / "refined.pdb"))
            report = python_report(pdb.topology, pdb.positions, destination / "python-v2.json")
            old_cli = json.loads((previous / cell / "full-context-cli.json").read_text())
            new_cli = json.loads((destination / "full-context-cli.json").read_text())
            old_mask = json.loads((previous / cell / "result.json").read_text())
            newly_exposed = difference(old_cli["findings"], report["findings"]) if old_cli["steric"]["actual_engine"] == steric.DETERMINISTIC_STERIC_ENGINE else None
            row = {"cell": cell, "entry": case.name, "serialized_final": True,
                "historical_mask_and_motion_passed": old_mask["restoration_mask_and_motion_passed"],
                "v2_mask_and_motion_passed": fresh["restoration_mask_and_motion_passed"],
                "full_cli_engine": new_cli["steric"], "full_cli_counts": fresh["full_cli_counts"],
                "full_cli_introduced": fresh["introduced_cli_findings"],
                "newly_exposed_python_findings": newly_exposed}
        else:
            snapshot = Path(job["command"][job["command"].index("--snapshot")+1])
            topology, _, _ = read_temporary(snapshot)
            with np.load(work / "phase-4/returned-state.npz", allow_pickle=False) as data:
                points = data["positions_nm"].copy()
            report = python_report(topology, points*nanometer, destination / "rejected-native-python-v2.json")
            # Native states include auxiliary H/OXT, unlike the heavy PDBs.
            # Do not compare their finding counts to serialized candidate counts.
            row = {"cell": cell, "entry": case.name, "serialized_final": False,
                "native_scope": "recorded full FF topology with auxiliary H/OXT; not a serialized candidate",
                "v2_mask_and_motion_passed": False, "rejection": fresh["rejection"],
                "D_ca_findings": fresh["D_ca_findings"]}
        baseline = source_cache[case.name]["prepared-baseline"]["findings"]
        row.update(python_connectivity=report["steric"],
                   python_counts=dict(Counter(f["severity"] for f in report["findings"])),
                   classified_findings=[classify(f, generated, anchor_atoms) for f in report["findings"]],
                   v2_not_in_prepared_baseline=difference(baseline, report["findings"]) if original["completed"] else None,
                   saved_force_diagnostic=force_rows[cell]["force_state"],
                   force_diagnostic_above_configured_scale=force_rows[cell]["force_state"]["mobile_gradient_rms_over_all_components"] > 10)
        save(destination / "v2-comparison.json", row)
        rows.append(row)
    summary = []
    for entry in plan["admitted"]:
        selected = [r for r in rows if r["entry"] == entry]
        finals = [r for r in selected if r["serialized_final"]]
        summary.append({"entry": entry, "cells": len(selected), "serialized_candidates": len(finals),
            "historical_mask_and_motion_passes": sum(r["historical_mask_and_motion_passed"] for r in finals),
            "v2_mask_and_motion_passes": sum(r["v2_mask_and_motion_passed"] for r in finals),
            "partial_serialized_connectivity": sum(r["python_connectivity"]["status"] != "complete" for r in finals),
            "finals_with_immutable_severe_clashes": sum(any(f["only_immutable_atoms"] and f["finding"]["category"] == "clash" and f["finding"]["severity"] == "ERROR" for f in r["classified_findings"]) for r in finals),
            "finals_with_generated_severe_clashes": sum(any(f["touches_generated"] and f["finding"]["category"] == "clash" and f["finding"]["severity"] == "ERROR" for f in r["classified_findings"]) for r in finals),
            "finals_above_force_scale": sum(r["force_diagnostic_above_configured_scale"] for r in finals),
            "full_cli_no_introduced_findings": sum(not r["full_cli_introduced"] for r in finals),
            "newly_exposed_finding_range": [min(len(r["newly_exposed_python_findings"] or []) for r in finals), max(len(r["newly_exposed_python_findings"] or []) for r in finals)]})
    for name, expected in {**historical, **bindings}.items():
        if digest(Path(name)) != expected:
            raise ValueError(f"bound evidence/code changed during reanalysis: {name}")
    save(output / "summary.json", {"cases": summary, "source_reports": source_cache,
        "declaration_sha256": digest(output / "declaration.json"), "optimizer_calls": 0,
        "historical_reports_verified_unchanged": True, "new_policy_selected": False,
        "force_scale_note": "saved diagnostic only, not an exposed OpenMM convergence reason",
        "native_scope_note": "six rejected FF states include temporary H/OXT; counts are not pooled with 94 heavy PDB candidates"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--forces", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.grid.resolve(), args.previous.resolve(), args.forces.resolve(), args.output.resolve())
