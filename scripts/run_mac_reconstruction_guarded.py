#!/usr/bin/env python3
"""Run ONE declared Mac scientific smoke with retained resource/device evidence.

Parent requires psutil. Child uses a small explicit entry wrapper to observe
PyTorch module tensors on MPS; no weights/algorithms or device APIs are replaced.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

CHILD = r'''
import atexit,json,runpy,sys,threading,time
from pathlib import Path
import torch
root=Path(sys.argv[1]); module=sys.argv[2]; arguments=sys.argv[3:]
if not torch.backends.mps.is_available():raise RuntimeError("MPS unavailable")
report={"torch":torch.__version__,"mps_available":True,"automatic_cpu_fallback":False,
        "model_mps_calls":0,"module_classes":{},"peak_mps_driver_bytes":0}
def tensors(x):
    if isinstance(x,torch.Tensor):yield x
    elif isinstance(x,dict):
        for v in x.values():yield from tensors(v)
    elif isinstance(x,(tuple,list)):
        for v in x:yield from tensors(v)
def hook(module,args,output):
    cls=type(module).__module__+"."+type(module).__name__
    if cls.startswith(("rfd3.","boltz.","protenix.")) and any(t.device.type=="mps" for t in tensors(output)):
        report["model_mps_calls"]+=1
        report["module_classes"][cls]=report["module_classes"].get(cls,0)+1
handle=torch.nn.modules.module.register_module_forward_hook(hook)
stop=threading.Event()
def snapshot():
    while not stop.wait(0.5):
        report["peak_mps_driver_bytes"]=max(report["peak_mps_driver_bytes"],torch.mps.driver_allocated_memory())
        (root/"device-progress.json").write_text(json.dumps(report,indent=2)+"\n")
        if report["peak_mps_driver_bytes"]>6*1024**3:
            (root/"device-limit-exceeded.json").write_text(json.dumps(report,indent=2)+"\n")
            import os,signal
            os.kill(os.getpid(),signal.SIGTERM)
thread=threading.Thread(target=snapshot,daemon=True);thread.start()
def save():
    stop.set();handle.remove()
    report["peak_mps_driver_bytes"]=max(report["peak_mps_driver_bytes"],torch.mps.driver_allocated_memory())
    (root/"device-final.json").write_text(json.dumps(report,indent=2)+"\n")
atexit.register(save)
sys.argv=[module,*arguments]
runpy.run_module(module,run_name="__main__")
if report["model_mps_calls"]==0:raise RuntimeError("No actual model MPS tensors observed")
'''


def resource_stop_reason(*, elapsed, timeout, peak_rss, pressure, critical_seconds, low_headroom_seconds):
    """Swap use alone is telemetry, not evidence of current memory pressure."""
    if elapsed > timeout:
        return "wall-time-limit"
    if peak_rss > 8 * 1024**3:
        return "aggregate-process-rss-limit"
    if pressure >= 4 and critical_seconds >= 5:
        return "sustained-critical-host-memory-pressure"
    if low_headroom_seconds >= 5:
        return "sustained-host-memory-headroom"
    return None


def main():
    import psutil

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--module", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--resource-policy", choices=("legacy-swap-v1", "pressure-v2"), default="pressure-v2")
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.output.absolute()
    lock_path = Path('.artifacts/mac-reconstruction-2026-11-11/ml.lock')
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        (root / "child-entry.py").write_text(CHILD)
        rest = args.args[1:] if args.args and args.args[0] == "--" else args.args
        command = [str(args.python.absolute()), str(root / "child-entry.py"), str(root), args.module, *rest]
        env = os.environ.copy()
        env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTORCH_ENABLE_MPS_FALLBACK="0",
                   WANDB_MODE="disabled", OMP_NUM_THREADS="1", TOKENIZERS_PARALLELISM="false")
        env.pop("PYTORCH_MPS_HIGH_WATERMARK_RATIO", None)
        env.pop("PYTORCH_MPS_LOW_WATERMARK_RATIO", None)
        (root / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        start = time.monotonic()
        swap_start = psutil.swap_memory().used
        peak_rss = 0
        reason = None
        critical_since = None
        low_since = None
        telemetry = []
        with (root / "output.log").open("xb") as output:
            process = subprocess.Popen(command, env=env, stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
            while process.poll() is None:
                try:
                    target = psutil.Process(process.pid)
                    rss = target.memory_info().rss + sum(p.memory_info().rss for p in target.children(recursive=True))
                    peak_rss = max(peak_rss, rss)
                except (psutil.NoSuchProcess, psutil.ZombieProcess):
                    pass
                now = time.monotonic()
                available = psutil.virtual_memory().available
                pressure_result = subprocess.run(["sysctl", "-n", "kern.memorystatus_vm_pressure_level"],
                                                 capture_output=True, text=True, timeout=5)
                if pressure_result.returncode or not pressure_result.stdout.strip().isdigit():
                    reason = "host-pressure-telemetry-unavailable"
                    pressure = -1
                else:
                    pressure = int(pressure_result.stdout.strip())
                    critical_since = (critical_since if critical_since is not None else now) if pressure >= 4 else None
                    low_since = (low_since if low_since is not None else now) if available < 1024**3 else None
                    reason = resource_stop_reason(elapsed=now-start, timeout=args.timeout, peak_rss=peak_rss,
                        pressure=pressure, critical_seconds=0 if critical_since is None else now-critical_since,
                        low_headroom_seconds=0 if low_since is None else now-low_since)
                if args.resource_policy == "legacy-swap-v1" and not reason:
                    if available < 1024**3:
                        reason = "host-memory-headroom"
                    elif psutil.swap_memory().used - swap_start > 512 * 1024**2:
                        reason = "host-swap-growth"
                telemetry.append({"elapsed_seconds": now-start, "rss_peak_bytes": peak_rss,
                    "available_bytes": available, "pressure_level": pressure,
                    "swap_growth_bytes": psutil.swap_memory().used - swap_start})
                if (root / "device-limit-exceeded.json").exists():
                    reason = "mps-driver-limit"
                if reason:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    break
                time.sleep(0.5)
            returncode = process.wait()
        report = {"returncode": returncode, "guard_stop_reason": reason,
            "resource_policy": args.resource_policy,
            "wall_seconds": time.monotonic() - start, "peak_aggregate_rss_bytes": peak_rss,
            "swap_growth_bytes": psutil.swap_memory().used - swap_start,
            "scientific_acceptance": False, "scope": "technical execution, independent scientific audit still required"}
        (root / "execution.json").write_text(json.dumps(report, indent=2) + "\n")
        (root / "host-telemetry.json").write_text(json.dumps(telemetry, indent=2) + "\n")
        print(json.dumps(report), flush=True)
        if returncode or reason:
            sys.exit(1)


if __name__ == "__main__":
    main()
