"""Private 1NON ordinary MODELLER, then matched fixed/local refinement controls.

Execution and scientific acceptance are separate. No new diffusion sampling,
public capabilities or GUI status changes. Every rejected candidate is retained.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.anchor_calibration_preflight import anchor_residues
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import key, launch_worker
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

SOURCE = Path("tmp/1non-backends/attempt01/protein-only.pdb")
REQUEST = Path("tmp/1non-backends/attempt01/context-request.json")
DIFFUSION = Path("tmp/1non-backends/attempt01/packed103-ff/case/proposal.pdb")


def local_mask(request, source, radius=4.):
    """One source-only mask for both proposal methods, never candidate-selected."""
    flanks = {(r.chain, r.residue_number, r.insertion_code) for r in anchor_residues(request, 2)}
    reference = np.array([row["xyz"] for a, row in source.items() if a[:3] in flanks])
    if not len(reference):
        raise ValueError("no observed flanks")
    residues = set(flanks)
    for a, row in source.items():
        if np.any(np.sum((reference - row["xyz"])**2, axis=-1) <= radius**2):
            residues.add(a[:3])
    return tuple(sorted(a for a in source if a[:3] in residues))


def movements(source, candidate, selected):
    grouped = {}
    selected = set(selected)
    outside = 0.
    for a, row in source.items():
        d = math.dist(row["xyz"], candidate[a]["xyz"])
        if a in selected:
            grouped.setdefault(a[:3], {})[a[3]] = d
        else:
            outside = max(outside, d)
    result = []
    for residue, displacements in sorted(grouped.items()):
        bb = math.sqrt(sum(displacements[n]**2 for n in ("N", "CA", "C", "O"))/4)
        ca, maximum = displacements["CA"], max(displacements.values())
        result.append({"residue": list(residue), "backbone_rms_A": bb, "CA_A": ca, "heavy_max_A": maximum,
                       "within_limits": bb <= .5 and ca <= .5 and maximum <= 1.25})
    return result, outside


def render(request, source_text, rows, *, include_observations=False):
    axis = (*request.fixed_atoms, *request.generated_atoms) if include_observations else request.generated_atoms
    return materialize_candidate_pdb(source_text, axis, [rows[key(a)]["xyz"] for a in axis],
                                     tuple(rows[key(a)]["resname"] for a in axis), request)


def scientific_status(geometry_passed, coordinate_passed, admission_failures):
    """Process exit codes deliberately have no role in scientific acceptance."""
    passed = bool(geometry_passed and coordinate_passed and not admission_failures)
    return passed, "accepted" if passed else "rejected"


def assess(request, source_path, candidate, allowed=()):
    """Exactly the same full-context + request-scoped gate for both backends."""
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l

    source, actual = atoms(source_path), atoms(candidate)
    expected = set(map(key, request.fixed_atoms)) | set(map(key, request.generated_atoms))
    if set(actual) != expected:
        raise ValueError("scientific assessment atom inventory changed")
    assessment = assessment_adapter(request, source_path.read_text())(candidate.read_text(),
        frozenset(a for a in request.fixed_atoms if key(a) in set(allowed)))
    findings = diagnose(candidate, candidate.with_suffix(".diagnose.json"))
    geometry_categories = {"clash", "chirality", "ca_planarity", "coincident_atoms", "bond_length", "bond_angle",
                           "non_planar_amide", "ramachandran", "ramachandran_outlier", "rotamer", "chain_break"}
    errors = [f for f in findings if f["severity"] == "ERROR" and f["category"] in geometry_categories]
    cis = [f for f in findings if f["category"] == "cis_peptide" and f.get("resname") != "PRO"]
    all_l = True
    structure = PDBFile(StringIO(candidate.read_text()))
    try:
        assert_all_l(structure.topology, structure.positions)
    except ChiralityError:
        all_l = False
    motion, outside = movements(source, actual, allowed)
    metadata_exact = all(source[a]["line"][:30] == actual[a]["line"][:30]
                         and source[a]["line"][54:] == actual[a]["line"][54:] for a in source)
    geometry_passed = assessment.geometry_passed and not errors and all_l and not cis
    coordinate_passed = outside == 0 and all(m["within_limits"] for m in motion) and metadata_exact
    passed, status = scientific_status(geometry_passed, coordinate_passed, assessment.admission_failures)
    report = {"validation_passed": passed, "quality_status": status,
        "geometry_passed": geometry_passed, "coordinate_policy_passed": coordinate_passed,
        "request_geometry": asdict(assessment), "full_geometry_error_counts": dict(Counter(f["category"] for f in errors)),
        "non_PRO_cis_findings": cis, "global_all_l": all_l, "observed_noncoordinate_metadata_exact": metadata_exact,
        "warning_counts": dict(Counter(f["category"] for f in findings if f["severity"] == "WARNING")),
        "outside_mobile_mask_max_A": outside, "local_observed_movements": motion,
        "all_observed_max_A": max(math.dist(r["xyz"], actual[a]["xyz"]) for a, r in source.items()),
        "readiness_error_counts": dict(Counter(f["category"] for f in findings if f["severity"] == "ERROR"
                                               and f["category"] not in geometry_categories)),
        "candidate_sha256": digest(candidate), "thresholds": {"clash_warning_overlap_A": .6, "clash_error_overlap_A": .9}}
    save(candidate.with_suffix(".validation.json"), report)
    return report


def environment():
    return {**os.environ, "PYTHONPATH": str(Path.cwd()), "PYTHONNOUSERSITE": "1", "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def baseline(args):
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    request = DiffusionRequest.from_json(REQUEST.read_text())
    if digest(SOURCE) != request.normalized_pdb.sha256 or len(request.fixed_atoms) != 5040 or len(request.generated_atoms) != 556:
        raise ValueError("requires the retained protein-only 1NON case")
    source = atoms(SOURCE)
    selected = local_mask(request, source)
    (root / "source.pdb").write_bytes(SOURCE.read_bytes())
    (root / "request.json").write_bytes(REQUEST.read_bytes())
    save(root / "local-mask.json", {"source_only": True, "flank_width": 2, "source_radius_A": 4.,
        "atoms": [list(a) for a in selected], "residues": [list(r) for r in sorted({a[:3] for a in selected})],
        "limits": {"backbone_rms_per_residue_A": .5, "CA_max_A": .5, "heavy_max_A": 1.25}})
    files = [SOURCE, REQUEST, DIFFUSION, Path(__file__), Path("scripts/local_gap_rosetta_worker.py"),
             Path("scripts/one_non_modeller_worker.py"), Path("scripts/run_retained_loop_modeller_controls.py"),
             Path("scripts/audit_one_non_local.py"), Path("scripts/run_mac_modeller_proposal.py"),
             Path("scripts/mac_reconstruction_comparison.py"), Path("scripts/hybrid_gap_reconstruction.py"),
             Path("scripts/retained_evidence.py"), Path("scripts/run_mac_reconstruction_guarded.py"),
             Path("tests/test_one_non_local.py"), root / "source.pdb", root / "request.json", root / "local-mask.json"]
    files += list(Path("src/dvbfixer/model/diffusion").glob("*.py"))
    save(root / "plan.json", {"protocol": "private-1non-common-validation-local-release-v1",
        "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in files},
        "ordinary_MODELLER": "one seed101 initial model + one fast loop, normal initial all-atom optimization; native retained",
        "diffusion": "retained RFD3 seed103 packed proposal; source-bound; not a raw neural endpoint or fresh sampling",
        "new_diffusion_calls": 0, "new_MODELLER_proposals": 1, "phases": ["ordinary-and-validation", "fixed-controls", "local-controls"],
        "paired_refinement": "fixed and local begin from identical source-restored proposal, same pack seed and optimizer budget; not a continuation",
        "cells": [f"{m}-{p}" for p in ("fixed", "local") for m in ("modeller", "rfd3")],
        "scorefunction": "ref2015_cart", "seed": 101, "max_minimizer_iterations": 500,
        "local_observed_atoms": len(selected), "local_observed_residues": len({a[:3] for a in selected}),
        "local_restraints": "coordinate harmonic sd0.2A plus flat-bound walls: backbone0.5A, sidechain1.25A, sd0.05A; final per-residue RMS/CA/heavy audit",
        "amide_restraints": "generated/junction trans omega CircularHarmonic pi, sd0.15rad, same in both cells; no retrospective cis repair",
        "outside_mask": "fixed Cartesian atoms; disposable H mobile; no solvent", "stationarity_claimed": False,
        "acceptance": "same core expanded geometry, full-context geometry errors, global L, no cis non-PRO, source metadata/outside-mask/movement bounds",
        "execution_success_is_not_validation_success": True, "public_capability_or_GUI_changed": False,
        "timeout_seconds_per_worker": args.timeout, "resources": "serial CPU1thread, pressure-v2/RSS8GiB"})
    freeze_code(root / "plan.json")
    execution_dir = root / "ordinary-execution"
    execution_dir.mkdir()
    command = [str(args.python.absolute()), "-m", "scripts.one_non_modeller_worker", str(root / "request.json"),
               str(root / "source.pdb"), str(root / "ordinary-proposal")]
    save(execution_dir / "command.json", command)
    execution = launch_worker(command, execution_dir, environment(), args.timeout)
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("ordinary MODELLER did not complete; no local phase launched")
    ordinary = atoms(root / "ordinary-proposal/mapped-native.pdb")
    # Native coordinates are preserved: only identities/metadata are normalized.
    native = root / "ordinary-native.pdb"
    native.write_text(render(request, SOURCE.read_text(), ordinary, include_observations=True))
    reports = {"modeller-native": assess(request, root / "source.pdb", native)}
    for method, rows in (("modeller", ordinary), ("rfd3", atoms(DIFFUSION))):
        candidate = root / f"{method}-source-restored.pdb"
        candidate.write_text(render(request, SOURCE.read_text(), rows))
        reports[f"{method}-source-restored"] = assess(request, root / "source.pdb", candidate)
    save(root / "ordinary-completed.json", {"execution_status": "completed", "MODELLER_exit_code": 0,
         "geometry_and_coordinate_validation": reports, "local_refinement_executed": False,
         "plan_sha256": digest(root / "plan.json")})
    verify_bindings(json.loads((root / "plan.json").read_text())["source_code_and_input_sha256"])
    print(json.dumps({"execution_status": "completed", "validation": {k: {s: v[s] for s in
          ("quality_status", "geometry_passed", "coordinate_policy_passed", "full_geometry_error_counts", "all_observed_max_A")} for k, v in reports.items()}}))


def refine(args, *, local_retry=False):
    root = args.output.resolve()
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (root / "plan.json",))
    baseline_report = json.loads((root / "ordinary-completed.json").read_text())
    if baseline_report["plan_sha256"] != digest(root / "plan.json"):
        raise ValueError("local phase requires this completed ordinary run")
    registration = root / "registered-baseline"
    ready = json.loads((registration / "ready.json").read_text())
    if ready["parent_plan_sha256"] != digest(root / "plan.json"):
        raise ValueError("registered baseline belongs to a different declaration")
    registration_plan = json.loads((registration / "plan.json").read_text())
    verify_bindings(registration_plan["source_code_and_input_sha256"], (registration / "plan.json",))
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = atoms(root / "source.pdb")
    selected = tuple(map(tuple, json.loads((root / "local-mask.json").read_text())["atoms"]))
    if selected != local_mask(request, source):
        raise ValueError("source-only local mask does not reproduce")
    destination = root / "local-retry" if local_retry else root
    if local_retry:
        destination.mkdir(exist_ok=False)
    declaration = destination if local_retry else root / "refinement-declaration"
    if not local_retry:
        declaration.mkdir(exist_ok=False)
    bindings = [registration / "plan.json", registration / "ready.json", root / "ordinary-completed.json",
                root / "source.pdb", root / "request.json", root / "local-mask.json",
                registration / "modeller-source-restored.pdb", registration / "rfd3-source-restored.pdb",
                Path(__file__).resolve(), Path("scripts/local_gap_rosetta_worker.py").resolve(),
                Path("scripts/audit_one_non_local.py").resolve()]
    if local_retry:
        bindings += [root / name / "native-failure.json" for name in ("modeller-local", "rfd3-local")]
    save(declaration / "plan.json", {"protocol": "private-1non-registered-common-refinement-v1",
         "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in bindings},
         "parent_plan_sha256": digest(root / "plan.json"), "budgets_and_mask_unchanged": True,
         "new_MODELLER_or_diffusion_calls": 0, "cell_order": ["modeller-local", "rfd3-local"] if local_retry else plan["cells"],
         "retry_reason": "correct BoundFunc namespace after pre-packing failure; no parameter/seed/mask change" if local_retry else None,
         "registration": "one proper global observed-CA transform before explicit source restoration; not local refitting"})
    freeze_code(declaration / "plan.json")
    results = {}
    for policy in (("local",) if local_retry else ("fixed", "local")):
        for method in ("modeller", "rfd3"):
            folder = destination / f"{method}-{policy}"
            folder.mkdir(exist_ok=False)
            starting = registration / f"{method}-source-restored.pdb"
            (folder / "starting.pdb").write_bytes(starting.read_bytes())
            allowed = selected if policy == "local" else ()
            payload = {"policy": policy, "generated": [list(key(a)) for a in request.generated_atoms],
                "source": [{"identity": list(a), "xyz": row["xyz"]} for a, row in source.items()],
                "observed_mobile": [list(a) for a in allowed], "seed": 101, "max_iterations": 500}
            save(folder / "payload.json", payload)
            command = [str(args.rosetta_python.absolute()), str(Path("scripts/local_gap_rosetta_worker.py").resolve()), str(folder)]
            save(folder / "command.json", command)
            execution = launch_worker(command, folder, environment(), args.timeout)
            record = {"execution_status": "completed" if execution["returncode"] == 0 and execution["guard_stop_reason"] is None else "failed",
                      "execution_exit_code": execution["returncode"], "quality_status": "not_evaluated"}
            if record["execution_status"] == "completed":
                native = json.loads((folder / "native.json").read_text())
                coords = {tuple(r["identity"]): np.array(r["xyz"]) for r in native["atoms"]}
                expected = set(source) | set(map(key, request.generated_atoms))
                if len(coords) != len(native["atoms"]) or set(coords) != expected:
                    raise ValueError("native refinement atom inventory changed")
                immutable = set(source) - set(allowed)
                error = max(math.dist(source[a]["xyz"], coords[a]) for a in immutable)
                if error > 1e-7:
                    raise ValueError("native outside-window motion before serialization")
                updates = {a: coords[key(a)] for a in (*request.generated_atoms, *request.fixed_atoms)
                           if key(a) in set(allowed) | set(map(key, request.generated_atoms))}
                text = rewrite_generated_coordinates(starting.read_text(), updates)
                candidate = folder / "candidate.pdb"
                candidate.write_text(text)
                record.update(assess(request, root / "source.pdb", candidate, allowed))
                native_motion, native_outside = movements(source, {a: {"xyz": p} for a, p in coords.items()}, allowed)
                native_within = native_outside <= 1e-7 and all(row["within_limits"] for row in native_motion)
                record["native_local_observed_movements"] = native_motion
                record["native_movement_bounds_passed"] = native_within
                record["coordinate_policy_passed"] = record["coordinate_policy_passed"] and native_within
                record["validation_passed"] = record["validation_passed"] and native_within
                record["quality_status"] = "accepted" if record["validation_passed"] else "rejected"
                record["native_outside_mask_max_A"] = error
                record["native_score"] = native["score"]
            results[f"{method}-{policy}"] = record
            save(folder / "result.json", record)
            print(json.dumps({"cell": f"{method}-{policy}", "execution_status": record["execution_status"],
                  "quality_status": record["quality_status"], "errors": record.get("full_geometry_error_counts")}))
    save(destination / "refinement-completed.json", {"execution_status": "completed", "new_diffusion_calls": 0,
         "ordinary_run_preceded_refinement": True, "cells": results, "public_qualification": False})


def register(args):
    """Correct the global coordinate gauge of the saved ordinary output only."""
    from scripts.audit_idpforge_research import register_ca

    root = args.output.resolve()
    parent_plan = json.loads((root / "plan.json").read_text())
    verify_bindings(parent_plan["source_code_and_input_sha256"], (root / "plan.json",))
    if not (root / "ordinary-completed.json").is_file():
        raise ValueError("registration requires completed ordinary validation")
    folder = root / "registered-baseline"
    folder.mkdir(exist_ok=False)
    source = atoms(root / "source.pdb")
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    raw = root / "ordinary-proposal/mapped-native.pdb"
    proposal = atoms(raw)
    axis = list(proposal)
    observed_ca = [a for a in source if a[3] == "CA"]
    registered, transform = register_ca(np.array([proposal[a]["xyz"] for a in axis]),
        np.array([proposal[a]["xyz"] for a in observed_ca]), np.array([source[a]["xyz"] for a in observed_ca]))
    transformed = {a: {**proposal[a], "xyz": p} for a, p in zip(axis, registered, strict=True)}
    save(folder / "plan.json", {"protocol": "private-ordinary-source-frame-registration-v1",
        "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in (
            root / "plan.json", root / "ordinary-completed.json", raw, root / "source.pdb", root / "request.json",
            DIFFUSION, Path(__file__), Path("scripts/audit_idpforge_research.py"))},
        "transform": transform, "correspondence": "all 644 observed CA in the intact assembly; no per-chain fits",
        "new_sampling_or_minimization_calls": 0, "old_unregistered_reports_retained": True})
    freeze_code(folder / "plan.json")
    reports = {}
    candidate = folder / "ordinary-native-registered.pdb"
    candidate.write_text(render(request, (root / "source.pdb").read_text(), transformed, include_observations=True))
    reports["modeller-native-registered"] = assess(request, root / "source.pdb", candidate)
    for method, rows in (("modeller", transformed), ("rfd3", atoms(DIFFUSION))):
        candidate = folder / f"{method}-source-restored.pdb"
        candidate.write_text(render(request, (root / "source.pdb").read_text(), rows))
        reports[f"{method}-source-restored"] = assess(request, root / "source.pdb", candidate)
    save(folder / "ready.json", {"parent_plan_sha256": digest(root / "plan.json"), "transform": transform,
        "new_sampling_or_minimization_calls": 0, "geometry_and_coordinate_validation": reports})
    print(json.dumps({"global_observed_CA_fit": transform, "validation": {name: {k: report[k] for k in
          ("quality_status", "geometry_passed", "coordinate_policy_passed", "full_geometry_error_counts", "all_observed_max_A")}
          for name, report in reports.items()}}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("baseline", "register", "refine", "retry-local"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(os.sys.executable))
    parser.add_argument("--rosetta-python", type=Path, default=Path("/Users/lilpoundcake/micromamba/envs/pyrosetta/bin/python"))
    parser.add_argument("--timeout", type=int, default=900, choices=range(1, 1801))
    args = parser.parse_args()
    {"baseline": baseline, "register": register, "refine": refine,
     "retry-local": lambda a: refine(a, local_retry=True)}[args.stage](args)
