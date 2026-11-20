"""Bounded saved-torsion N-terminal reconstruction versus one MODELLER proposal.

No diffusion execution, observed-coordinate variables, reference coordinates,
closure sweep or production change. Both arms use the same ideal reverse growth,
generated-only packing and Cartesian minimization policy.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_idpforge_research import audit as audit_idpforge
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import (
    audit_native,
    build_backbone_from_right,
    diagnostic_counts,
    dihedral,
    key,
    launch_worker,
)
from scripts.idpforge_research import admit
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def terminal_torsions(order, proposal, length):
    """Only proposal-internal torsions; the first residue's undefined phi is unused."""
    def xyz(i, name):
        r = order[i]
        return proposal[r.chain, r.residue_number, r.insertion_code, name]["xyz"]
    torsions = []
    for i in range(length):
        phi = dihedral(xyz(i-1, "C"), xyz(i, "N"), xyz(i, "CA"), xyz(i, "C")) if i else 0.
        psi = dihedral(xyz(i, "N"), xyz(i, "CA"), xyz(i, "C"), xyz(i+1, "N"))
        torsions.append([phi, psi])
    anchor_phi = dihedral(xyz(length-1, "C"), xyz(length, "N"), xyz(length, "CA"), xyz(length, "C"))
    return torsions, anchor_phi


def prepare_backbone(request, source, proposal):
    transport = admit(request, source)
    order = target_axis(request)[transport["chain"]]
    length = transport["gap_length"]
    if length > 8:
        raise ValueError("initial terminal experiment is bounded to at most eight generated residues")
    fixed, generated = set(map(key, request.fixed_atoms)), set(map(key, request.generated_atoms))
    if set(proposal) != fixed | generated:
        raise ValueError("proposal atom axis differs from the request")
    for a, row in proposal.items():
        expected = next(r["resname"] for r in transport["residues"]
                        if (r["chain"], r["residue_number"], r["insertion_code"]) == a[:3])
        if row["resname"] != expected or not np.isfinite(row["xyz"]).all():
            raise ValueError("proposal sequence/coordinates changed")
    torsions, anchor_phi = terminal_torsions(order, proposal, length)
    anchor = order[length]
    right = [source[anchor.chain, anchor.residue_number, anchor.insertion_code, n]["xyz"] for n in ("N", "CA", "C")]
    rows = build_backbone_from_right(right, torsions, anchor_phi)
    lines = []
    serial = max(int(row["line"][6:11]) for row in source.values())
    for i, residue in enumerate(transport["residues"]):
        rk = residue["chain"], residue["residue_number"], residue["insertion_code"]
        if i < length:
            for name in ("N", "CA", "C", "O"):
                serial += 1
                xyz = rows[i][name]
                lines.append(f'ATOM  {serial:5d} {name:^4s} {residue["resname"]} {rk[0]}{int(rk[1]):4d}{rk[2] or " "}   '
                             f'{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           {name[0]}  \n')
        else:
            lines.extend(row["line"] for a, row in source.items() if a[:3] == rk)
    payload = {"operation": "pack-and-cartmin", "gaps": [], "residues": [asdict(r) for r in order],
               "fixed": [{"identity": list(a), "xyz": r["xyz"]} for a, r in source.items()],
               "generated": [list(key(a)) for a in request.generated_atoms],
               "candidates": 1, "seed": 42, "max_iterations": 500}
    report = {"torsions_radians": torsions, "anchor_phi_radians": anchor_phi,
              "first_generated_phi": "undefined and unused; stored as zero", "ideal_rows": [{k: v.tolist() for k, v in r.items()} for r in rows],
              "boundary_C_N_A": float(np.linalg.norm(rows[-1]["C"] - right[0])),
              "source_variables": 0, "reference_coordinates_read": False, "cis_peptides_sampled": False}
    return "".join(lines) + "TER\nEND\n", payload, report


def materialize(request, source_text, source, native, output):
    coordinates, maximum = audit_native({"fixed": [{"identity": list(a), "xyz": r["xyz"]} for a, r in source.items()],
        "generated": [list(key(a)) for a in request.generated_atoms]}, native)
    transport = admit(request, source)
    names = {(r["chain"], r["residue_number"], r["insertion_code"]): r["resname"] for r in transport["residues"]}
    text = materialize_candidate_pdb(source_text, request.generated_atoms,
        [coordinates[key(a)] for a in request.generated_atoms],
        tuple(names[key(a)[:3]] for a in request.generated_atoms), request)
    output.write_text(text)
    actual = atoms(output)
    if set(actual) != set(source) | set(map(key, request.generated_atoms)) or any(actual[a]["line"] != r["line"] for a, r in source.items()):
        raise ValueError("materialization changed source records or ownership")
    assessment = assessment_adapter(request, source_text)(text)
    save(output.with_suffix(".assessment.json"), asdict(assessment))
    findings = diagnose(output, output.with_suffix(".diagnose.json"))
    # Global chirality is checked independently of the generated-only assessment.
    from io import StringIO

    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l

    structure = PDBFile(StringIO(text))
    all_l = True
    try:
        assert_all_l(structure.topology, structure.positions)
    except ChiralityError:
        all_l = False
    return {"geometry_passed": assessment.geometry_passed and all_l, "admission_failures": list(assessment.admission_failures),
            "diagnose": diagnostic_counts(findings), "native_fixed_max_A": maximum, "source_records_exact": True,
            "global_all_l": all_l, "candidate_sha256": digest(output)}


def run(args):
    root, sample = args.output.resolve(), args.sample.resolve()
    root.mkdir(parents=True, exist_ok=False)
    args.owned_output = True
    parent_audit = audit_idpforge(sample)
    parent_plan = json.loads((sample / "plan.json").read_text())
    if not parent_audit["sampling_completed"] or parent_plan.get("template_ss_policy") != "radian-consistent":
        raise ValueError("require audited completed radian-consistent saved IDPForge sample")
    source, req = sample / "source.pdb", sample / "request.json"
    request = DiffusionRequest.from_json(req.read_text())
    source_rows = atoms(source)
    admit(request, source_rows)
    raw_path = sample / "raw-extraction/native-mapped.pdb"
    from scripts.audit_idpforge_research import audit_extraction
    audit_extraction(raw_path.parent)
    for name, path in (("source.pdb", source), ("request.json", req)):
        (root / name).write_bytes(path.read_bytes())
    rosetta_worker = Path("scripts/hybrid_gap_rosetta.py").resolve()
    modeller_worker = Path("scripts/terminal_modeller_worker.py").resolve()
    files = [source, req, sample / "plan.json", sample / "native-endpoint.npz", raw_path,
             Path(__file__).resolve(), rosetta_worker, modeller_worker,
             Path("scripts/audit_terminal_reconstruction.py").resolve(),
             Path("scripts/hybrid_gap_reconstruction.py").resolve(), Path("scripts/idpforge_research.py").resolve(),
             Path("scripts/audit_idpforge_research.py").resolve(), Path("scripts/retained_evidence.py").resolve(),
             Path("scripts/mac_reconstruction_comparison.py").resolve(),
             Path("scripts/run_mac_reconstruction_guarded.py").resolve(),
             Path("tests/test_terminal_reconstruction.py").resolve(),
             *Path("src/dvbfixer/model/diffusion").glob("*.py")]
    save(root / "plan.json", {"protocol": "private-terminal-reverse-growth-comparison-v1", "new_diffusion_calls": 0,
         "sample": str(sample), "source_code_and_input_sha256": {str(p.resolve()): digest(p) for p in files},
         "arms": ["idpforge", "modeller"], "proposals_per_arm": 1, "MODELLER_initial_and_loop_models": [1, 1],
         "MODELLER_seed": 42, "rosetta_seed": 42, "scorefunction": "ref2015_cart",
         "postprocessing": "same reverse ideal growth, generated-sidechain packing, generated-heavy-only Cartesian LBFGS",
         "max_minimizer_iterations_per_arm": 500, "random_restarts": 0, "closure": False,
         "variable_observed_atoms": 0, "auxiliaries": "disposable Rosetta H are mobile; auxiliary OXT fixed, excluded from output",
         "solvent": False, "reference_coordinates_read": False, "proposal_selection": False, "public_capability": False,
         "raw_MODELLER_baseline_retained_and_assessed": True, "timeout_seconds_per_worker": args.timeout,
         "resources": "serial CPU1thread, pressure-v2/RSS8GiB"})
    freeze_code(root / "plan.json")
    env = {**os.environ, "PYTHONPATH": str(Path.cwd()), "PYTHONNOUSERSITE": "1",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    modeller_dir = root / "modeller-native"
    modeller_dir.mkdir()
    (modeller_dir / "source.pdb").write_bytes(source.read_bytes())
    (modeller_dir / "request.json").write_bytes(req.read_bytes())
    execution = launch_worker([str(args.python.absolute()), str(modeller_worker), str(modeller_dir)], modeller_dir, env, args.timeout)
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("single MODELLER proposal failed; evidence retained, no fallback")
    modeller = atoms(modeller_dir / "mapped-native.pdb")
    native_motion = max(math.dist(r["xyz"], modeller[a]["xyz"]) for a, r in source_rows.items())
    save(modeller_dir / "observed-audit.json", {"native_observed_max_A": native_motion, "source_copied_before_audit": False})
    baseline_native = {"atoms": [{"identity": list(a), "xyz": (source_rows[a] if a in source_rows else r)["xyz"]}
                                  for a, r in modeller.items()]}
    baseline = materialize(request, source.read_text(), source_rows, baseline_native, root / "modeller-source-restored-baseline.pdb")
    baseline["native_MODELLER_observed_max_A"] = native_motion
    baseline["post_restoration_fixed_max_A"] = baseline.pop("native_fixed_max_A")
    baseline["note"] = "explicit source restoration, not a claim of native MODELLER fixation"
    save(root / "modeller-baseline-summary.json", baseline)
    results = {}
    for name, proposal in (("idpforge", atoms(raw_path)), ("modeller", modeller)):
        folder = root / name
        folder.mkdir()
        text, payload, report = prepare_backbone(request, source_rows, proposal)
        (folder / "input-backbone.pdb").write_text(text)
        save(folder / "payload.json", payload)
        save(folder / "construction.json", report)
        execution = launch_worker([str(args.rosetta_python.absolute()), str(rosetta_worker), str(folder)], folder, env, args.timeout)
        if execution["returncode"] or execution["guard_stop_reason"] is not None:
            results[name] = {"completed": False, "execution": {k: v for k, v in execution.items() if k != "telemetry"}}
            continue
        stages = {}
        for stage, filename in (("packed", "packed-native.json"), ("minimized", "native.json")):
            native = json.loads((folder / "candidate-00" / filename).read_text())
            stages[stage] = materialize(request, source.read_text(), source_rows, native, folder / f"{stage}.pdb")
        results[name] = {"completed": True, "stages": stages}
        save(folder / "summary.json", results[name])
    verify_bindings(json.loads((root / "plan.json").read_text())["source_code_and_input_sha256"])
    summary = {"new_diffusion_calls": 0, "MODELLER_proposals": 1, "source_atoms": len(source_rows),
               "generated_atoms": len(request.generated_atoms), "baseline": baseline, "arms": results,
               "public_qualification": False, "broader_superiority_claimed": False}
    save(root / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sample", "output", "python", "rosetta-python"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600, choices=range(1, 1801))
    arguments = parser.parse_args()
    try:
        run(arguments)
    except Exception as exc:
        import traceback
        if getattr(arguments, "owned_output", False):
            save(arguments.output / "failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
