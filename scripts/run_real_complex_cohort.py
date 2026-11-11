#!/usr/bin/env python3
"""Execute admitted, prospectively frozen real-complex cases via managed GUI API.

No source refusals are replaced. The normal public CLI prepares each new input
and freezes its own final scientific request; the preflight request is evidence,
not a hidden GUI flag or a bypass around new-request preparation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    protocol = json.loads((root / "protocol.json").read_text())
    if protocol["seeds"] != [7, 101, 102, 103, 104] or protocol["profile"] != "protpardelle-1c-cc94-mps-v1":
        raise ValueError("unexpected cohort or seeds")
    if not json.loads((root / "runtime-preflight.json").read_text())["passed"]:
        raise ValueError("runtime preflight did not pass")
    preparation = json.loads((root / "preparation-summary.json").read_text())
    audits = {r["pdb"]: r for r in json.loads((root / "natural-gap-audit.json").read_text())}
    admitted = [r for r in preparation if r["status"] == "admitted"]
    for row in admitted:
        if not audits[row["pdb"]]["passed"] or sha(root / row["request"]) != row["request_sha256"]:
            raise ValueError("natural-gap audit or frozen preflight request binding failed")
    for relative, expected in protocol["source_code_sha256"].items():
        if sha(repo / relative) != expected:
            raise ValueError(f"source changed after cohort preparation: {relative}")
    run = root / "managed-execution"
    run.mkdir(exist_ok=False)
    shutil.copyfile(__file__, run / "run-script.py")
    # Outer HTTP/CLI limits encompass all five seeds and all per-gap samplings;
    # no 120-second interactive tool timeout can kill a long refinement.
    runner_seconds = 23 * 3600
    save(run / "execution-protocol.json", {
        "cohort_protocol_sha256": sha(root / "protocol.json"),
        "preparation_summary_sha256": sha(root / "preparation-summary.json"),
        "natural_gap_audit_sha256": sha(root / "natural-gap-audit.json"),
        "run_script_sha256": sha(__file__), "python": sys.executable,
        "built_server_sha256": sha(repo / "gui/dist-server/server.js"),
        "admitted_structures": [r["pdb"] for r in admitted],
        "refused_structures": [r["pdb"] for r in preparation if r["status"] != "admitted"],
        "seeds": protocol["seeds"], "candidate_budget": len(admitted) * len(protocol["seeds"]),
        "per_gap_sampler_budget": sum(r["minimum_denoiser_invocations"] for r in admitted),
        "runner_timeout_seconds": runner_seconds, "child_timeout_seconds": 24 * 3600,
        "entry_point": "POST /api/v1/workspaces/{workspaceId}/jobs",
        "fallback": False, "resampling": False, "parallel_sampling": False,
    })
    environment = {**os.environ, "DVBFIXER_HOST": "127.0.0.1", "DVBFIXER_PORT": "0",
        "DVBFIXER_GUI_DATA_DIR": str(run / "data"), "DVBFIXER_GUI_ROOT": str(repo / "gui"),
        "DVBFIXER_EXECUTABLE": sys.executable,
        "DVBFIXER_ARGS": json.dumps(["-c", "from dvbfixer.cli import main; main()"]),
        "DVBFIXER_TIMEOUT_MS": str(24 * 3600 * 1000),
        "DVBFIXER_MAX_CONCURRENT_PROCESSES": "1", "OPENMM_CPU_THREADS": "1",
        "OMP_NUM_THREADS": "1", "PYTHONNOUSERSITE": "1", "PYTHONPATH": str(repo / "src"),
        "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]}
    # Do not expose an inherited user-configured remote-auth credential through
    # this explicitly local acceptance service or its retained launch metadata.
    for key in ("DVBFIXER_AUTH_PRINCIPALS", "DVBFIXER_ALLOW_INSECURE_REMOTE"):
        environment.pop(key, None)
    node = shutil.which("node")
    if node is None:
        raise ValueError("Node executable is unavailable")
    process = subprocess.Popen([node, "dist-server/server.js"], cwd=repo / "gui", env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    ready = threading.Event()
    origin = []

    def drain(stream, name):
        with (run / name).open("x") as log:
            for line in stream:
                log.write(line)
                log.flush()
                match = re.search(r"listening on (http://127\.0\.0\.1:\d+)", line)
                if match:
                    origin.append(match.group(1))
                    ready.set()

    for stream, name in ((process.stdout, "server.stdout.log"), (process.stderr, "server.stderr.log")):
        threading.Thread(target=drain, args=(stream, name), daemon=True).start()

    def api(route, method="GET", body=None):
        data = None if body is None else json.dumps(body).encode()
        request = urllib.request.Request(origin[0] + route, data=data, method=method,
                                         headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)

    results = []
    active = None
    try:
        if not ready.wait(30):
            raise RuntimeError("standalone GUI service did not start")
        save(run / "service.json", {"origin": origin[0], "server_pid": process.pid,
             "orchestrator_pid": os.getpid(), "data_root": str(run / "data")})
        print("GUI acceptance service:", origin[0], flush=True)
        for row in preparation:
            pdb = row["pdb"]
            if row["status"] != "admitted":
                results.append({"pdb": pdb, "status": "admission_refused", "sampling_executed": False,
                                "error": row.get("error"), "requested_seeds": protocol["seeds"]})
                continue
            workspace = api("/api/workspaces", "POST", {"name": f"Real complex {pdb} · cc94 MPS"})
            directory = run / "data/projects" / workspace["id"]
            # Source is the unchanged pre-completion protein projection. Public
            # input preparation retains completion/repair provenance in this job.
            shutil.copyfile(root / "cases" / pdb / "protein.pdb", directory / f"{pdb}.pdb")
            body = {"command": "model", "inputFile": f"{pdb}.pdb",
                "fastaContent": (root / "cases" / pdb / "target.fasta").read_text(),
                "values": {"--backend": "diffusion", "--diffusion-model": "protpardelle",
                    "--diffusion-seed": "7,101,102,103,104", "--diffusion-timeout": runner_seconds}}
            save(run / f"{pdb}-submission.json", body)
            active = api(f"/api/v1/workspaces/{workspace['id']}/jobs", "POST", body)
            save(run / f"{pdb}-accepted-job.json", active)
            route = f"/api/v1/workspaces/{workspace['id']}/jobs/{active['id']}"

            def collect_events(event_route, name):
                try:
                    request = urllib.request.Request(origin[0] + event_route + "/events")
                    with urllib.request.urlopen(request, timeout=24 * 3600) as response, (run / name).open("xb") as stream:
                        for line in response:
                            stream.write(line)
                            stream.flush()
                except Exception as exc:
                    save(run / (name + ".error.json"), {"error": str(exc)})

            events = threading.Thread(target=collect_events, args=(route, f"{pdb}-events.sse"), daemon=True)
            events.start()
            print(pdb, "started", active["id"], "gaps", row["gap_count"], flush=True)
            deadline = time.monotonic() + 25 * 3600
            while active["status"] in {"queued", "running"}:
                if process.poll() is not None:
                    raise RuntimeError("GUI server exited while executing the cohort")
                if time.monotonic() > deadline:
                    raise TimeoutError("outer scientific-job deadline exceeded")
                time.sleep(3)
                active = api(route)
            events.join(timeout=10)
            save(run / f"{pdb}-finished-job.json", active)
            results.append({"pdb": pdb, "status": active["status"], "job_id": active["id"],
                "workspace_id": workspace["id"], "requested_seeds": protocol["seeds"],
                "diffusion": active.get("diffusion"), "error": active.get("error"),
                "scientific_acceptance": "inspect per-seed independent validation, not process exit"})
            temporary = run / "results.next.json"
            temporary.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n")
            temporary.replace(run / "results.json")
            print(pdb, "finished", active["status"], flush=True)
            active = None
        save(run / "completed.json", {"completed": True, "results": results,
            "scope": "current cc94 real-gap complexes; no reconstruction accuracy or broad qualification claim"})
    except Exception as exc:
        save(run / "execution-error.json", {"error_type": type(exc).__name__, "error": str(exc)})
        raise
    finally:
        if active and active["status"] in {"queued", "running"}:
            try:
                api(f"/api/v1/workspaces/{active['workspaceId']}/jobs/{active['id']}", "DELETE")
            except (OSError, urllib.error.URLError):
                pass
        process.terminate()
        # Server graceful shutdown drains managed scientific diagnostics. Do not
        # kill it with an interactive-shell timeout during result packaging.
        try:
            process.wait(timeout=300)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


if __name__ == "__main__":
    main()
