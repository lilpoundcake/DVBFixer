#!/usr/bin/env python3
"""Declare and execute a bounded, serial MODELLER comparison, retaining failures."""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import time
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import audit, inventory
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save
from scripts.refine_mac_reconstruction_smoke import run as refine
from scripts.run_mac_reconstruction_guarded import resource_stop_reason

PROJECT = Path(__file__).resolve().parents[1]
HISTORICAL = PROJECT / "tmp/mac-reconstruction-2026-11-11"


def declare(output: Path, python: Path, recover_native: Path | None = None):
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    output.mkdir(parents=True, exist_ok=False)
    # Existing tasks: deliberately retrospective comparator, not validation.
    cases = [HISTORICAL / "smoke-input-03",
        HISTORICAL / "reference-panel-preparation/1A3N-A1-59-gap30-37",
        HISTORICAL / "reference-panel-preparation/1FR2-A3-85-gap24-31",
        HISTORICAL / "longer-reference-panel-preparation/1A3N-A1-59-gap30-45",
        HISTORICAL / "longer-reference-panel-preparation/1FR2-A3-85-gap24-39"]
    files = [PROJECT / f"scripts/{name}.py" for name in ("run_mac_modeller_comparison",
        "run_mac_modeller_proposal", "mac_reconstruction_comparison", "audit_mac_reconstruction_smoke",
        "refine_mac_reconstruction_smoke", "run_mac_reconstruction_guarded")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files.append(PROJECT / "src/dvbfixer/model/modeller_run.py")
    if recover_native is not None:
        files.append(recover_native)
    jobs = []
    for historical in cases:
        manifest = json.loads((historical / "manifest.json").read_text())
        for name, expected in manifest["files"].items():
            if digest(historical / name) != expected:
                raise ValueError(f"changed historical case: {historical / name}")
        case = output / "cases" / (historical.name if historical.name != "smoke-input-03" else "1UBQ-A1-72-gap10-13")
        case.mkdir(parents=True)
        for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb"):
            (case / name).write_bytes((historical / name).read_bytes())
        if (historical / "evaluation-mask.json").is_file():
            (case / "evaluation-mask.json").write_bytes((historical / "evaluation-mask.json").read_bytes())
        else:
            request = build_cli_diffusion_request(case / "source.pdb", {manifest.get("chain", "A"): manifest["sequence"]},
                seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
            (case / "evaluation-mask.json").write_text(request.to_json())
        manifest = {**manifest, "historical_manifest": str(historical / "manifest.json"),
            "historical_manifest_sha256": digest(historical / "manifest.json"),
            "files": {name: digest(case / name) for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb", "evaluation-mask.json")}}
        save(case / "manifest.json", manifest)
        files += [case / name for name in (*manifest["files"], "manifest.json")]
        for seed in DEVELOPMENT_SEEDS:
            job = {"case": str(case), "case_id": case.name, "seed": seed}
            if case.name == "1UBQ-A1-72-gap10-13" and seed == 7 and recover_native is not None:
                job["recover_native"] = str(recover_native)
            jobs.append(job)
    plan = {"protocol": "mac-reconstruction-modeller-comparison-v1", "phase": "retrospective baseline, not independent validation",
        "python": str(python), "python_sha256": digest(python), "jobs": jobs, "count": len(jobs),
        "generation": "MODELLER 10.8 full LoopModel; exact request PIR; rand_seed=-seed; 1 initial / 1 fast loop",
        "production_difference": "explicit seed, one loop rather than default two, no raw chirality rewrite, no failed-loop fallback",
        "selection": "single candidate per declared seed; no heldout-based selection",
        "observations": "native drift measured before restoration; source-restored stage reported separately",
        "repair": {"policy": "only source-restored geometry failures; all generated atoms together; fixed observations",
            "platform": "CPU", "threads": 1, "iterations": 1000, "restarts": 0, "optimizer_seed": 7,
            "geometry_preconditioning": False, "solvent": False},
        "metrics": ["native drift", "native/restored geometry and junctions", "final geometry", "global L",
            "heldout heavy/backbone RMSD", "generated movement during repair", "runtime and peak process RSS"],
        "native_vs_restored_vs_final_separate": True, "failures_in_denominator": True,
        "resource_policy": "pressure-v2 CPU arm, serial, 1800 seconds per complete job; 8 GiB process RSS; 10 GiB evidence cap",
        "reserved_validation_seeds_used": False, "automatic_retries": False, "qualification": False,
        "technical_amendment": None if recover_native is None else {
            "previous_failure": "attempt01 MODELLER Values/pdfterms JSON serialization after model.make returned",
            "recovery": "read back unchanged saved native loop PDB for seed7, no sampling; other 29 jobs not previously launched",
            "native_loop_status_metadata": "original opaque loop-output metadata unavailable; do not infer its score/status; independently audit complete saved proposal",
            "recovered_sha256": digest(recover_native)},
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}}
    save(output / "plan.json", plan)
    return plan


def worker(case: Path, output: Path, seed: int, recover_native: Path | None = None):
    from scripts.run_mac_modeller_proposal import run as proposal

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure

    if recover_native is None:
        proposal(case / "evaluation-mask.json", case / "source.pdb", output / "proposal", seed)
    else:
        from scripts.mac_reconstruction_comparison import normalize_modeller

        from dvbfixer.model.diffusion.contract import DiffusionRequest
        request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
        if seed != 7 or case.name != "1UBQ-A1-72-gap10-13":
            raise ValueError("recovery is declared only for the unchanged first native seed7 output")
        native = recover_native.read_bytes()
        normalized, transport = normalize_modeller(request, native.decode())
        (output / "proposal").mkdir()
        (output / "proposal/native.pdb").write_bytes(native)
        (output / "proposal/mapped-native.pdb").write_text(normalized)
        save(output / "proposal/mapping.json", transport)
        save(output / "proposal/recovery.json", {"original": str(recover_native), "sha256": digest(recover_native),
            "sampled_again": False, "opaque_loop_metadata_status": "unavailable; file independently audited"})
    report = audit(case, output / "proposal/mapped-native.pdb", output / "audit", seed=seed)
    refine(output / "audit", output / "shared-refinement", platform="CPU")
    final_path = output / "shared-refinement/refined.pdb"
    repaired = final_path.is_file()
    if not repaired:
        final_path = output / "audit/source-restored.pdb"
    final = inventory(final_path)
    reference = inventory(case / "reference-projection-NOT-CONDITIONING.pdb")
    source = inventory(case / "source.pdb")
    raw = inventory(output / "audit/source-restored.pdb")
    generated = sorted(set(reference)-set(source))
    backbone = [a for a in generated if a.atom_name in {"N", "CA", "C", "O"}]
    def rms(atoms, first, second):
        return float(np.sqrt(np.mean([np.sum((first[a][2]-second[a][2])**2) for a in atoms])))
    structure = _structure(final_path.read_text(), "final.pdb")
    global_l = True
    try:
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
    except Exception as exc:
        global_l = False
        save(output / "global-chirality-failure.json", {"type": type(exc).__name__, "message": str(exc)})
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter
    request = DiffusionRequest.from_json((output / "audit/evaluation-mask.json").read_text())
    assessed = assessment_adapter(request, (case / "source.pdb").read_text())(final_path.read_text())
    row = {"case": case.name, "seed": seed, "backend": "modeller", "repaired": repaired,
        "raw_restored_geometry_passed": report["assessments"]["source-restored"]["geometry_passed"],
        "native_observed_drift": report["native_observed_drift_after_observed_only_fit_angstrom"],
        "final_geometry_passed": assessed.geometry_passed, "final_admission_failures": list(assessed.admission_failures),
        "global_all_l": global_l, "final_source_geometry_and_l_passed": assessed.geometry_passed and not assessed.admission_failures and global_l,
        "raw_heldout_heavy_rmsd": report["heldout_generated_rmsd_angstrom"],
        "final_heldout_heavy_rmsd": rms(generated, final, reference),
        "final_heldout_backbone_rmsd": rms(backbone, final, reference),
        "generated_repair_rms_movement": rms(generated, final, raw),
        "final_path": str(final_path), "final_sha256": digest(final_path), "final_assessment": asdict(assessed),
        "qualification": False}
    save(output / "result.json", row)


def execute(output: Path):
    import fcntl

    import psutil

    plan = json.loads((output / "plan.json").read_text())
    for name, expected in plan["source_code_and_input_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"declared input/implementation changed: {name}")
    if digest(Path(plan["python"])) != plan["python_sha256"]:
        raise ValueError("declared interpreter changed")
    # Exclusive run reservation: a partial directory is never a completion record.
    (output / "execution").mkdir(exist_ok=False)
    rows = []
    with (PROJECT / ".artifacts/mac-reconstruction-2026-11-11/ml.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for job in plan["jobs"]:
            if sum(p.stat().st_size for p in output.rglob("*") if p.is_file()) > 10*1024**3:
                raise RuntimeError("10 GiB evidence cap exceeded; no more jobs launched")
            work = output / "execution" / f"{job['case_id']}-seed-{job['seed']}"
            work.mkdir()
            command = job.get("command") or [plan["python"], "-m", "scripts.run_mac_modeller_comparison", "worker",
                "--output", str(work), "--case", job["case"], "--seed", str(job["seed"])]
            if "recover_native" in job:
                command.extend(("--recover-native", job["recover_native"]))
            save(work / "command.json", command)
            env = {**os.environ, "PYTHONPATH": os.pathsep.join((str(PROJECT), str(PROJECT / "src"))),
                "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PYTHONNOUSERSITE": "1"}
            started = time.monotonic()
            peak, critical, low, reason = 0, None, None, None
            with (work / "output.log").open("xb") as log:
                process = subprocess.Popen(command, cwd=PROJECT, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    while process.poll() is None:
                        now = time.monotonic()
                        try:
                            parent = psutil.Process(process.pid)
                            rss = parent.memory_info().rss + sum(c.memory_info().rss for c in parent.children(recursive=True))
                            peak = max(peak, rss)
                        except (psutil.NoSuchProcess, psutil.ZombieProcess):
                            pass
                        pressure_result = subprocess.run(["sysctl", "-n", "kern.memorystatus_vm_pressure_level"], capture_output=True, text=True, timeout=5)
                        pressure = int(pressure_result.stdout.strip()) if pressure_result.returncode == 0 and pressure_result.stdout.strip().isdigit() else -1
                        critical = (critical if critical is not None else now) if pressure >= 4 else None
                        low = (low if low is not None else now) if psutil.virtual_memory().available < 1024**3 else None
                        reason = resource_stop_reason(elapsed=now-started, timeout=1800, peak_rss=peak, pressure=pressure,
                            critical_seconds=0 if critical is None else now-critical, low_headroom_seconds=0 if low is None else now-low)
                        if pressure == -1:
                            reason = "host-pressure-telemetry-unavailable"
                        if reason:
                            break
                        time.sleep(0.5)
                finally:
                    if process.poll() is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                    code = process.wait()
            row = {**job, "returncode": code, "elapsed_seconds": time.monotonic()-started,
                "peak_rss_bytes": peak, "guard_stop_reason": reason, "completed": (work / "result.json").is_file() and code == 0}
            if row["completed"]:
                row["result"] = json.loads((work / "result.json").read_text())
            save(work / "execution.json", row)
            rows.append(row)
            (output / "progress.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
            print(json.dumps({"case": job["case_id"], "seed": job["seed"], "completed": row["completed"], "returncode": code}), flush=True)
            if reason:
                raise RuntimeError(f"resource guard: {reason}; remaining jobs not launched")
            if code != 0 or not row["completed"]:
                raise RuntimeError("technical failure retained; remaining jobs not launched, no implicit retry")
    save(output / "completed.json", {"declared": len(plan["jobs"]), "executed": len(rows),
        "completed": sum(r["completed"] for r in rows), "qualification": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--python", type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--seed", type=int, choices=DEVELOPMENT_SEEDS)
    parser.add_argument("--recover-native", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.python is None:
            parser.error("declare requires --python")
        plan = declare(args.output.resolve(), args.python.absolute(), None if args.recover_native is None else args.recover_native.resolve())
        print(json.dumps({"count": plan["count"], "plan": str(args.output / "plan.json"), "amendment": plan["technical_amendment"]}, indent=2))
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if args.case is None or args.seed is None:
            parser.error("worker requires --case and --seed")
        try:
            worker(args.case.resolve(), args.output.resolve(), args.seed, None if args.recover_native is None else args.recover_native.resolve())
        except Exception:
            save(args.output / "failure.json", {"traceback": traceback.format_exc()})
            raise
