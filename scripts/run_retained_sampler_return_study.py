#!/usr/bin/env python3
"""Prospective seed11 transport study; preserve raw return, never skip heads."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scripts import run_retained_loop_context_sampling as previous
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_sampler_return_capture import passive_capture

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest

PROJECT, ROOT, RUNTIME, CLONE = previous.PROJECT, previous.ROOT, previous.RUNTIME, previous.CLONE


def declare(output):
    output.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), Path(previous.__file__), PROJECT / "scripts/retained_sampler_return_capture.py",
        PROJECT / "scripts/run_mac_reconstruction_guarded.py", PROJECT / "deploy/protenix-v1/mps_research_smoke.py",
        PROJECT / "deploy/protenix-v1/checkpoint_gap_smoke.py", PROJECT / "deploy/protenix-v1/reinjection.py",
        PROJECT / "deploy/protenix-v1/profile-lock.json"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    for folder in ("configs", "protenix/model", "protenix/data", "runner"):
        files += list((CLONE / folder).rglob("*.py"))
    old_readback = ROOT / "retained-loop-context-sampling-final-readback-attempt01/readback.json"
    if json.loads(old_readback.read_text())["actual_calls_used"] != 2:
        raise ValueError("previous two-call tranche must remain accounted")
    files.append(old_readback)
    jobs = []
    for entry in ("1NON", "2BE1"):
        original = previous.ENCODED / entry
        original_request = original / "context-request.json"
        request = DiffusionRequest.from_json(original_request.read_text())
        if request.seeds != (7,):
            raise ValueError("unexpected original seed declaration")
        derived = replace(request, seeds=(11,))
        before, after = json.loads(request.to_json()), json.loads(derived.to_json())
        before.pop("seeds")
        after.pop("seeds")
        if before != after:
            raise ValueError("seed-only derivative changed the scientific task")
        case = output / "inputs" / entry
        case.mkdir(parents=True, exist_ok=False)
        source = case / derived.normalized_pdb.path
        source.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original / request.normalized_pdb.path, source)
        if digest(source) != request.normalized_pdb.sha256:
            raise ValueError("source bytes changed in derivative")
        (case / "context-request.json").write_text(derived.to_json())
        save(case / "lineage.json", {"original_request_path": str(original_request),
            "original_request_sha256": digest(original_request), "original_seed": 7, "new_seed": 11,
            "only_request_change": "seeds", "encoded_input": str(original / "encoded/input.json"),
            "not_replay_or_retry": True})
        files += [p for p in original.rglob("*") if p.is_file()]
        files += [p for p in case.rglob("*") if p.is_file()]
        work = output / "execution" / entry
        jobs.append({"entry_id": entry, "work": str(work), "case": str(case),
            "command": [sys.executable, str(PROJECT / "scripts/run_mac_reconstruction_guarded.py"),
                "--python", str(RUNTIME / "envs/protenix/bin/python"), "--module", "scripts.run_retained_sampler_return_study",
                "--output", str(work), "--resource-policy", "pressure-v2", "--timeout", "1800", "--",
                "worker", "--case", str(case), "--output", str(work)]})
    save(output / "plan.json", {"jobs": jobs, "bindings": {str(p): digest(p) for p in files},
        "actual_call_budget": 2, "seed": 11, "updates": 200, "cycles": 1, "candidates": 1,
        "purpose": "passive sampler-return transport before downstream heads; not resource-limit calibration",
        "native_tasks_unchanged_except_seed": True, "context_tokens": [354, 365],
        "downstream_heads_skipped": False, "repair_calls": 0, "automatic_sampling_retries": 0,
        "guard": "unchanged serial RSS8GiB/MPS6GiB/1800s pressure-v2; no fallback",
        "original_seed7_failures_preserved": True, "new_sources_or_assets": False,
        "public_or_frozen_runners_changed": False, "not_final_evaluation_or_donor_approval": True})


def worker(case, output):
    sys.path.insert(0, str(CLONE))
    sys.path.insert(0, str(PROJECT / "deploy/protenix-v1"))
    import mps_research_smoke as legacy
    from protenix.model.protenix import Protenix

    path = case / "context-request.json"
    request = DiffusionRequest.from_json(path.read_text())
    source = case / request.normalized_pdb.path
    if request.seeds != (11,) or request.candidate_count != 1 or case.name not in {"1NON", "2BE1"}:
        raise ValueError("undeclared private seed/task")
    previous.validate_context_source(request, source.read_bytes())
    previous.validate_context_capability(request, "protenix")
    lineage = json.loads((case / "lineage.json").read_text())
    if digest(Path(lineage["original_request_path"])) != lineage["original_request_sha256"]:
        raise ValueError("original request changed")
    args = SimpleNamespace(request=path, input=Path(lineage["encoded_input"]), assets=RUNTIME / "assets/protenix",
        kalign=RUNTIME / "tools/bin/kalign", output=output / "predictions", seed=11)
    with previous.context_wiring(legacy, request), passive_capture(legacy, Protenix,
        output=output / "sampler-return", request=request, request_path=path, source=source):
        legacy.run(args, allowed_seeds=(11,))


def run(output):
    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared file changed: {name}")
    (output / "run-reservation").mkdir(exist_ok=False)
    env = {**os.environ, "PYTHONNOUSERSITE": "1", "PYTORCH_ENABLE_MPS_FALLBACK": "0", "LAYERNORM_TYPE": "torch",
        "PROTENIX_ROOT_DIR": str(RUNTIME / "assets/protenix"), "PYTHONPATH": os.pathsep.join((str(PROJECT), str(PROJECT / "src")))}
    rows = []
    for job in plan["jobs"]:
        result = subprocess.run(job["command"], cwd=PROJECT, env=env, check=False)
        work = Path(job["work"])
        receipt = json.loads((work / "execution.json").read_text())
        resource_failure = receipt["guard_stop_reason"] is not None
        row = {"entry_id": job["entry_id"], "returncode": result.returncode,
            "native_full_inference_success": result.returncode == 0,
            "resource_failure": receipt["guard_stop_reason"],
            "sampler_return_committed": (work / "sampler-return/capture.json").is_file()}
        rows.append(row)
        save(output / f"{job['entry_id']}-outcome.json", row)
        if result.returncode and not resource_failure:
            save(output / "blocked.json", {"outcomes": rows, "remaining_calls_not_launched": True,
                                          "automatic_retry": False})
            raise RuntimeError("technical failure retained; no retry")
    save(output / "completed.json", {"actual_calls_accounted": len(rows), "outcomes": rows,
        "scientific_acceptance": False, "independent_raw_readback_required": True})


def audit(root, output):
    from scripts.audit_retained_loop_context_sampling import fixed_records

    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    plan = json.loads((root / "plan.json").read_text())
    completion = json.loads((root / "completed.json").read_text())
    if completion["actual_calls_accounted"] != len(plan["jobs"]):
        raise ValueError("incomplete technical tranche")
    for name, expected in plan["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("execution declaration changed")
    rows, bindings = [], {}
    for job in plan["jobs"]:
        work = Path(job["work"])
        captured = work / "sampler-return"
        manifest_path = captured / "capture.json"
        if not manifest_path.is_file():
            rows.append({"entry_id": job["entry_id"], "raw_endpoint_available": False,
                         "geometry_outcome": "not_assessable"})
            continue
        manifest = json.loads(manifest_path.read_text())
        for name, expected in manifest["files"].items():
            p = captured / name
            if Path(name).name != name or digest(p) != expected:
                raise ValueError("capture payload digest/path mismatch")
            bindings[str(p)] = expected
        request_path, source = Path(manifest["request_path"]), Path(manifest["source_path"])
        request = DiffusionRequest.from_json(request_path.read_text())
        if digest(request_path) != manifest["request_sha256"] or digest(source) != manifest["source_sha256"]:
            raise ValueError("capture request/source changed")
        if request_path != Path(job["case"]) / "context-request.json" or request.seeds != (11,):
            raise ValueError("capture task/seed substitution")
        native_axis = json.loads((captured / "native-axis.json").read_text())
        axis = tuple(AtomIdentity(**atom) for atom in native_axis["atoms"])
        steps = json.loads((captured / "reinjection-steps.json").read_text())
        if (len(steps) != 200 or [s["step_index"] for s in steps] != list(range(200))
                or any(s["step_count"] != 200 or s["post_overwrite_max_angstrom"] != 0 for s in steps)):
            raise ValueError("raw callback trace incomplete")
        selected = {item.atom for item in request.protein_context.fixed_conditioning} | set(request.generated_atoms)
        indices = tuple(i for i, atom in enumerate(axis) if atom in selected)
        if (len(set(axis)) != len(axis) or set(axis[i] for i in indices) != selected
                or list(indices) != native_axis["materialized_indices"]):
            raise ValueError("native ownership/transport axis mismatch")
        native = np.load(captured / "native-arrays.npz", allow_pickle=False)["coordinates"]
        if native.shape != (len(axis), 3) or not np.isfinite(native).all():
            raise ValueError("native coordinate array invalid")
        by_atom = {atom: native[i] for i, atom in enumerate(axis)}
        if any(not np.array_equal(by_atom[item.atom], np.asarray(item.coordinate_angstrom, dtype=np.float32))
               for item in request.protein_context.fixed_conditioning):
            raise ValueError("native represented observations differ before materialization")
        from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb

        reproduced = materialize_candidate_pdb(source.read_text(), tuple(axis[i] for i in indices), native[list(indices)],
            tuple(native_axis["residue_names"][i] for i in indices), request,
            elements=tuple(native_axis["elements"][i] for i in indices))
        text = (captured / "raw-materialized.pdb").read_text()
        if text != reproduced or fixed_records(text, set(request.fixed_atoms)) != fixed_records(source.read_text(), set(request.fixed_atoms)):
            raise ValueError("independent materialization/source custody failed")
        coordinates, _ = _coordinates(text)
        if set(coordinates) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("complete endpoint atom mask changed")
        assessment = assessment_adapter(request, source.read_text())(text)
        structure, all_l, chirality_failure = _structure(text, "raw-return.pdb"), True, None
        try:
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
        except ChiralityError as exc:
            all_l, chirality_failure = False, str(exc)
        execution = json.loads((work / "execution.json").read_text())
        rows.append({"entry_id": job["entry_id"], "seed": 11, "raw_endpoint_available": True,
            "complete_observations_byte_exact": True, "native_axis_independently_materialized": True,
            "updates": len(steps), "assessment": asdict(assessment), "global_all_l": all_l,
            "chirality_failure": chirality_failure,
            "raw_geometry_passed": assessment.geometry_passed and not assessment.admission_failures and all_l,
            "full_inference_resource_failure": execution["guard_stop_reason"],
            "full_inference_returncode": execution["returncode"], "repair_calls": 0,
            "no_publication_or_manual_donor_approval": True})
        for p in (manifest_path, request_path, source, work / "execution.json"):
            bindings[str(p)] = digest(p)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "raw-readback.json", {"rows": rows, "bindings": bindings,
        "full_inference_acceptance_not_implied": True, "auditor_sha256": digest(Path(__file__))})
    print(json.dumps({"cases": len(rows), "raw_endpoints": sum(r["raw_endpoint_available"] for r in rows),
                      "raw_geometry_passes": sum(r.get("raw_geometry_passed", False) for r in rows)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker", "audit"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        run(args.output.resolve())
    elif args.mode == "audit":
        audit(args.root.resolve(), args.output.resolve())
    else:
        worker(args.case.resolve(), args.output.resolve())
