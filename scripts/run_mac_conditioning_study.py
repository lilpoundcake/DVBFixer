#!/usr/bin/env python3
"""Bounded private RFD3 A/C/D study: 30 calls, 24 assembled outcomes.

Full-heavy tiny-gap and joint-pair six-seed arms already exist and are reused.
This is a within-case mechanistic study, not independent qualification.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_conditioning_contracts import materialize, native_coordinates, prepare
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine
from scripts.run_rfd3_context_development_panel import sampling_arguments

PROJECT = Path(__file__).resolve().parents[1]
EVIDENCE = PROJECT / "tmp/mac-reconstruction-comparison-2026-11-11"
ARMS = (
    ("A-flanks8", "tiny", (0,), "flanks8", "all-heavy"),
    ("A-spatial6", "tiny", (0,), "spatial6", "all-heavy"),
    ("C-backbone", "tiny", (0,), "full", "backbone"),
    ("D-independent-gap0", "pair", (0,), "full", "all-heavy"),
    ("D-independent-gap1", "pair", (1,), "full", "all-heavy"),
)


def declare(output: Path):
    output.mkdir(parents=True, exist_ok=False)
    cases = {"tiny": EVIDENCE / "modeller-attempt02/cases/1UBQ-A1-72-gap10-13",
        "pair": EVIDENCE / "expansion-seed7-pilots-attempt01/inputs/1UBQ-A1-72-pair10-13_46-49"}
    descriptions = {}
    for name, kind, gaps, context, information in ARMS:
        descriptions[name] = prepare(cases[kind], output / "inputs" / name, gaps=gaps, context=context, information=information)
    original = json.loads((PROJECT / "tmp/mac-reconstruction-2026-11-11/rfd3-explicit-seed7/command.json").read_text())
    jobs = []
    for seed in DEVELOPMENT_SEEDS:
        for name, kind, *_ in ARMS:
            work = output / "runs" / f"{name}-seed{seed}"
            inputs = output / "inputs" / name
            arguments = sampling_arguments(original[4:], inputs=inputs, work=work, seed=seed)
            command = [sys.executable, "scripts/run_mac_reconstruction_guarded.py", "--python", original[0],
                "--module", original[3], "--output", str(work / "run"), "--resource-policy", "pressure-v2", "--timeout", "1800", "--", *arguments]
            jobs.append({"arm": name, "case": str(cases[kind]), "seed": seed, "work": str(work), "inputs": str(inputs), "command": command,
                "expected_fixed": descriptions[name]["sampler_fixed_atoms"], "expected_unfixed": descriptions[name]["sampler_unfixed_atoms"]})
    files = [PROJECT / "scripts" / name for name in ("run_mac_conditioning_study.py", "mac_conditioning_contracts.py",
        "mac_reconstruction_comparison.py", "run_mac_reconstruction_guarded.py", "rfd3_mps_seeded_research.py",
        "trace_rfd3_seed_transport.py", "refine_mac_reconstruction_smoke.py", "run_rfd3_context_development_panel.py", "preflight_mac_conditioning.py")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += [p for p in (output / "inputs").rglob("*") if p.is_file()]
    files += [case / name for case in cases.values() for name in ("source.pdb", "evaluation-mask.json", "reference-projection-NOT-CONDITIONING.pdb")]
    plan = {"protocol": "private-source-only-rfd3-conditioning-v1", "calls": 30, "assembled_outcomes": 24,
        "jobs": jobs, "arms": {name: {k: v for k, v in value.items() if k != "axis"} for name, value in descriptions.items()},
        "resource_policy": "pressure-v2, serial; 1800s/call; RSS8GiB; MPS6GiB; no new environments/assets",
        "sampler": "same 200 configured timesteps/two recycles/explicit MPS seed; measured calls reported separately",
        "baseline_reuse": {"tiny_full_heavy": "original explicit six-seed 1UBQ gap10-13 arm, no resampling",
            "pair_joint_full_heavy": "completed expansion six-seed joint-pair arm, no resampling"},
        "C_semantics": "observed sidechain coordinates withheld and sampled as nuisance atoms, never new output ownership; restore exact original observations",
        "D_semantics": "two independent calls each omit the other unknown gap entirely but retain every observed residue; explicit /0 entity cuts; coupling/task-size intervention, not weights-only causality",
        "selection": "same-seed independent donors combined without heldout scoring; jointly generated baseline kept atomic",
        "common_repair": "all target-generated atoms together; CPU1thread; seed7; max1000; no restarts/preconditioning/anchors; clean untouched",
        "reference_to_conditioning": False, "reserved_validation_seeds_used": False, "automatic_retries": False,
        "public_capability_enabled": False, "qualification": False, "files": {str(p): digest(p) for p in files}}
    save(output / "plan.json", plan)
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"calls": 30, "assembled_outcomes": 24, "shapes": plan["arms"]}, indent=2))


def run(output: Path):
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["files"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared conditioning input/code changed: {name}")
    for arm in plan["arms"]:
        inputs = output / "inputs" / arm
        receipt = json.loads((inputs / "preflight/preflight.json").read_text())
        if (receipt["manifest_sha256"] != digest(inputs / "manifest.json") or receipt["spec_sha256"] != digest(inputs / "rfd3-input.json")
                or not receipt["exact_source_identity_masks_checked"] or not receipt["unfixed_coordinates_zero"]):
            raise ValueError("pinned-builder preflight missing or differs from declared input")
    (output / "runs").mkdir(exist_ok=False)
    results, proposals = [], {}
    for job in plan["jobs"]:
        work, inputs = Path(job["work"]), Path(job["inputs"])
        env = {**os.environ, "PYTHONPATH": os.pathsep.join(str(PROJECT / p) for p in (".", "src", "scripts")),
            "DVBFIXER_RFD3_MODEL_SEED": str(job["seed"]), "DVBFIXER_RFD3_EXPECTED_FIXED": str(job["expected_fixed"]),
            "DVBFIXER_RFD3_EXPECTED_GENERATED": str(job["expected_unfixed"]), "DVBFIXER_RFD3_TRANSPORT_OUTPUT": str(work / "run/transport")}
        work.mkdir()
        save(work / "command.json", job["command"])
        with (work / "orchestration.log").open("xb") as log:
            code = subprocess.run(job["command"], env=env, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT).returncode
        row = {"arm": job["arm"], "seed": job["seed"], "returncode": code}
        execution_path = work / "run/execution.json"
        row["execution"] = json.loads(execution_path.read_text()) if execution_path.exists() else {}
        try:
            if code != 0 or not row["execution"] or row["execution"].get("guard_stop_reason"):
                raise RuntimeError("technical/resource failure; no implicit retry")
            coordinates, receipt = native_coordinates(inputs, work / "run")
            proposals[job["arm"], job["seed"]] = coordinates
            save(work / "native-custody.json", receipt)
            row["custody"] = receipt
        except Exception:
            row["failure"] = traceback.format_exc()
            save(work / "failure.json", row)
            results.append(row)
            (output / "sampling-results.json").write_text(json.dumps(results, indent=2)+"\n")
            raise
        results.append(row)
        (output / "sampling-results.json").write_text(json.dumps(results, indent=2)+"\n")
        print(json.dumps({"arm": job["arm"], "seed": job["seed"], "returncode": code}), flush=True)
    final_rows = []
    for seed in DEVELOPMENT_SEEDS:
        for arm in ("A-flanks8", "A-spatial6", "C-backbone", "D-independent-pair"):
            names = ("D-independent-gap0", "D-independent-gap1") if arm == "D-independent-pair" else (arm,)
            case = Path(next(j["case"] for j in plan["jobs"] if j["arm"] == names[0]))
            work = output / "assembled" / f"{arm}-seed{seed}"
            work.mkdir(parents=True)
            text = materialize(case, [proposals[name, seed] for name in names])
            request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
            (work / "source-restored.pdb").write_text(text)
            (work / request.normalized_pdb.path).write_bytes((case / "source.pdb").read_bytes())
            (work / "evaluation-mask.json").write_text(request.to_json())
            assess = assessment_adapter(request, (case / "source.pdb").read_text())
            raw = assess(text)
            save(work / "assembly-lineage.json", {"calls": names, "seed": seed, "atomic_joint_group_recombined": False,
                "all_generated_atoms": len(request.generated_atoms), "heldout_based_selection": False})
            refine(work, work / "common-repair", platform="CPU")
            final = work / "common-repair/refined.pdb"
            repaired = final.is_file()
            if not repaired:
                final = work / "source-restored.pdb"
            assessment = assess(final.read_text())
            structure = _structure(final.read_text(), "conditioning-final.pdb")
            all_l = True
            try:
                assert_all_l(structure.pdb.topology, structure.pdb.positions)
            except Exception:
                all_l = False
            # This is the first held-out coordinate read in assembly/refinement.
            reference, coords = inventory(case / "reference-projection-NOT-CONDITIONING.pdb"), inventory(final)
            errors = {a: float(np.sum((reference[a][2]-coords[a][2])**2)) for a in request.generated_atoms}
            final_rows.append({"arm": arm, "seed": seed, "case": str(case), "raw_geometry_passed": raw.geometry_passed,
                "raw_admission_failures": list(raw.admission_failures), "repaired": repaired,
                "final_geometry_passed": assessment.geometry_passed, "final_admission_failures": list(assessment.admission_failures),
                "global_all_l": all_l, "heldout_heavy_rmsd": float(np.sqrt(np.mean(list(errors.values())))),
                "heldout_backbone_rmsd": float(np.sqrt(np.mean([v for a, v in errors.items() if a.atom_name in {"N", "CA", "C", "O"}]))),
                "final_path": str(final), "final_sha256": digest(final), "assessment": asdict(assessment), "qualification": False})
            (output / "final-progress.json").write_text(json.dumps(final_rows, indent=2)+"\n")
    save(output / "final-results.json", final_rows)
    save(output / "completed.json", {"calls": len(results), "assembled_outcomes": len(final_rows), "qualification": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    (declare if args.mode == "declare" else run)(args.output.resolve())
