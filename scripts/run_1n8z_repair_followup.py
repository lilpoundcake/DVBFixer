#!/usr/bin/env python3
"""Run the approved five 1N8Z seeds after the unchanged baseline releases MPS."""
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
import urllib.request
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--gui-root", type=Path, required=True)
    args = parser.parse_args()
    cohort, gui = args.cohort.resolve(), args.gui_root.resolve()
    repo = Path(__file__).resolve().parents[1]
    work = cohort / "1N8Z-managed-followup"
    work.mkdir(exist_ok=False)
    source = cohort / "cases/1N8Z/protein.pdb"
    fasta = cohort / "cases/1N8Z/target.fasta"
    private = cohort / "1N8Z-followup-runtime/runners"
    launcher = private / "dvbfixer-diffusion-protpardelle-1c-cc94-mps-v1"
    preparation = cohort / "1N8Z-composed-repair-disulfide-followup"
    if not json.loads((preparation / "independent-current-audit.json").read_text())["passed"]:
        raise ValueError("independent preparation audit did not pass")
    if not json.loads((cohort / "1N8Z-followup-runtime/preflight.json").read_text())["passed"]:
        raise ValueError("private runtime preflight did not pass")
    if not (preparation / "python-fast-passed.log").is_file():
        raise ValueError("successful full-suite evidence is absent")
    hashes = {str(path.relative_to(repo)): sha(path) for path in [
        *(repo / "src/dvbfixer/model/diffusion").glob("*.py"),
        repo / "src/dvbfixer/model/diffusion_cli.py", repo / "deploy/protpardelle-1c/cc94_runner.py"]}
    save(work / "protocol.json", {"case": "1N8Z", "seeds": [7, 101, 102, 103, 104],
        "sampling_updates": 500, "step_scale": 1.0, "s_churn": 0, "gap_count": 4,
        "candidate_budget": 5, "per_gap_sampling_budget": 20,
        "source_sha256": sha(source), "fasta_sha256": sha(fasta),
        "run_script_sha256": sha(__file__),
        "code_sha256": hashes, "launcher_sha256": sha(launcher),
        "built_server_sha256": sha(gui / "dist-server/server.js"),
        "software_validation_log_sha256": sha(preparation / "python-fast-passed.log"),
        "original_refusal_retained": True, "same_software_as_baseline": False,
        "replacement_seeds": False, "fallback": False,
        "wait_for_baseline": str(cohort / "managed-execution/completed.json")})
    shutil.copyfile(__file__, work / "run-script.py")
    print("1N8Z: waiting for the original nine-complex cohort to finish; no new MPS sampling yet", flush=True)
    # This is durable orchestration, not a user-facing success claim. No sampling
    # is permitted before the baseline has published its terminal index and exited.
    deadline = time.monotonic() + 7 * 24 * 3600
    owner = json.loads((cohort / "managed-execution/service.json").read_text())["orchestrator_pid"]
    while True:
        if time.monotonic() > deadline:
            raise TimeoutError("baseline did not finish within the scheduling window")
        completed = cohort / "managed-execution/completed.json"
        is_complete = completed.is_file()
        if is_complete:
            if not json.loads(completed.read_text()).get("completed"):
                raise ValueError("baseline completion index is invalid")
        try:
            os.kill(owner, 0)
        except ProcessLookupError:
            if is_complete:
                break
            save(work / "scheduling-error.json", {"error": "baseline exited without a completion index; no sampling started"})
            raise RuntimeError("baseline exited without a completion index")
        time.sleep(30)
    frozen = json.loads((work / "protocol.json").read_text())
    if (any(sha(repo / relative) != digest for relative, digest in hashes.items())
            or sha(launcher) != frozen["launcher_sha256"] or sha(source) != frozen["source_sha256"]
            or sha(fasta) != frozen["fasta_sha256"] or sha(__file__) != frozen["run_script_sha256"]
            or sha(gui / "dist-server/server.js") != frozen["built_server_sha256"]):
        raise ValueError("follow-up software or inputs changed after its protocol was frozen")
    environment = {**os.environ, "PYTHONPATH": str(repo / "src"),
        "PATH": str(private) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"],
        "DVBFIXER_HOST": "127.0.0.1", "DVBFIXER_PORT": "0", "DVBFIXER_GUI_ROOT": str(gui),
        "DVBFIXER_GUI_DATA_DIR": str(work / "data"), "DVBFIXER_EXECUTABLE": sys.executable,
        "DVBFIXER_ARGS": json.dumps(["-c", "from dvbfixer.cli import main; main()"]),
        "DVBFIXER_TIMEOUT_MS": str(24 * 3600 * 1000), "DVBFIXER_MAX_CONCURRENT_PROCESSES": "1",
        "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1", "PYTHONNOUSERSITE": "1"}
    for key in ("DVBFIXER_AUTH_PRINCIPALS", "DVBFIXER_ALLOW_INSECURE_REMOTE"):
        environment.pop(key, None)
    process = subprocess.Popen(["node", "dist-server/server.js"], cwd=gui, env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    origin = []
    ready = threading.Event()

    def drain(stream, name):
        with (work / name).open("x") as log:
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
        request = urllib.request.Request(origin[0] + route, method=method,
            data=None if body is None else json.dumps(body).encode(), headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)

    job = None
    try:
        if not ready.wait(30):
            raise RuntimeError("follow-up GUI service did not start")
        save(work / "service.json", {"origin": origin[0], "server_pid": process.pid,
            "orchestrator_pid": os.getpid(), "data_root": str(work / "data")})
        workspace = api("/api/workspaces", "POST", {"name": "1N8Z · composed input repair · five seeds"})
        directory = work / "data/projects" / workspace["id"]
        shutil.copyfile(source, directory / "1N8Z.pdb")
        submission = {"command": "model", "inputFile": "1N8Z.pdb", "fastaContent": fasta.read_text(),
            "values": {"--backend": "diffusion", "--diffusion-model": "protpardelle",
                       "--diffusion-seed": "7,101,102,103,104", "--diffusion-timeout": 23 * 3600}}
        save(work / "submission.json", submission)
        job = api(f"/api/v1/workspaces/{workspace['id']}/jobs", "POST", submission)
        save(work / "accepted-job.json", job)
        route = f"/api/v1/workspaces/{workspace['id']}/jobs/{job['id']}"
        end = time.monotonic() + 25 * 3600
        with (work / "job-observations.jsonl").open("x") as log:
            while job["status"] in {"queued", "running"}:
                log.write(json.dumps(job) + "\n")
                log.flush()
                if time.monotonic() > end or process.poll() is not None:
                    raise RuntimeError("follow-up service exited or exceeded its outer deadline")
                time.sleep(3)
                job = api(route)
        save(work / "finished-job.json", job)
        save(work / "completed.json", {"completed": True, "job": job,
            "scientific_acceptance": "inspect independent per-seed reports; process success is not a geometry pass"})
        print("1N8Z follow-up finished:", job["status"], flush=True)
    finally:
        if job and job["status"] in {"queued", "running"}:
            try:
                api(f"/api/v1/workspaces/{job['workspaceId']}/jobs/{job['id']}", "DELETE")
            except OSError:
                pass
        process.terminate()
        try:
            process.wait(timeout=300)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


if __name__ == "__main__":
    main()
