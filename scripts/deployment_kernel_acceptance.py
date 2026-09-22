"""Privileged, destructive-on-disposable-host checks for the systemd profile.

Run only on a disposable Linux VM with a dedicated /var/lib/dvbfixer mount.
The CI job provisions that mount before invoking this script as root.
"""

from __future__ import annotations

import argparse
import errno
import os
import pathlib
import subprocess
import sys
import time


ROOT = pathlib.Path("/var/lib/dvbfixer")
SCRIPT = pathlib.Path(__file__).resolve()


def cgroup_file(name: str) -> pathlib.Path:
    entry = next(line for line in pathlib.Path("/proc/self/cgroup").read_text().splitlines() if line.startswith("0::"))
    return pathlib.Path("/sys/fs/cgroup") / entry[3:].lstrip("/") / name


def counter(name: str, key: str) -> int:
    values = dict(line.split() for line in cgroup_file(name).read_text().splitlines())
    return int(values[key])


def child(mode: str) -> None:
    if mode == "cpu":
        before = counter("cpu.stat", "nr_throttled")
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            pass
        assert counter("cpu.stat", "nr_throttled") > before, "CPU quota never throttled"
    elif mode == "memory":
        before = counter("memory.events", "oom_kill")
        workers = []
        for _ in range(2):
            workers.append(subprocess.Popen([
                sys.executable, "-c",
                "import time; b=bytearray(90*1024*1024); b[::4096]=b'x'*(len(b)//4096); time.sleep(5)",
            ]))
        for worker in workers:
            worker.wait(timeout=15)
        assert counter("memory.events", "oom_kill") > before, "combined children escaped MemoryMax"
    elif mode == "pids":
        before = counter("pids.events", "max")
        workers = []
        denied = False
        try:
            for _ in range(100):
                try:
                    workers.append(subprocess.Popen(["/usr/bin/sleep", "10"]))
                except OSError as exc:
                    assert exc.errno == errno.EAGAIN, exc
                    denied = True
                    break
            assert denied and counter("pids.events", "max") > before, "TasksMax did not block forks"
        finally:
            for worker in workers:
                worker.terminate()
            for worker in workers:
                worker.wait(timeout=5)
    elif mode in {"persistent", "temporary"}:
        target = ROOT / "quota-probe" if mode == "persistent" else pathlib.Path("/tmp/quota-probe")
        denied = False
        try:
            with target.open("wb", buffering=0) as stream:
                for _ in range(512):
                    stream.write(b"x" * 1024 * 1024)
        except OSError as exc:
            assert exc.errno == errno.ENOSPC, exc
            denied = True
        finally:
            target.unlink(missing_ok=True)
        assert denied, f"{mode} storage did not reach its boundary"
    elif mode == "detached":
        detached = subprocess.Popen(["/usr/bin/sleep", "300"], start_new_session=True)
        (ROOT / "detached.pid").write_text(str(detached.pid))
        while True:
            time.sleep(60)
    else:
        raise ValueError(mode)


def systemd_run(mode: str, *properties: str, wait: bool = True) -> subprocess.CompletedProcess[str]:
    command = ["systemd-run", "--quiet", f"--unit=dvbfixer-accept-{mode}",
               *[f"--property={item}" for item in properties]]
    if wait:
        command += ["--wait", "--collect", "--pipe"]
    command += [sys.executable, str(SCRIPT), "--child", mode]
    return subprocess.run(command, text=True, capture_output=True, timeout=45, check=wait)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--child", choices=["cpu", "memory", "pids", "persistent", "temporary", "detached"])
    args = parser.parse_args()
    if args.child:
        child(args.child)
        return
    assert os.geteuid() == 0, "run on a disposable VM as root"
    assert os.path.ismount(ROOT), "dedicated bounded filesystem is required"
    common = ["ProtectSystem=strict", "ReadWritePaths=/var/lib/dvbfixer", "KillMode=control-group"]
    systemd_run("cpu", *common, "CPUQuota=20%")
    print("CPU throttling: PASS", flush=True)
    systemd_run("memory", *common, "MemoryMax=128M", "MemorySwapMax=0")
    print("combined child memory: PASS", flush=True)
    systemd_run("pids", *common, "TasksMax=20")
    print("descendant task limit: PASS", flush=True)
    systemd_run("persistent", *common)
    print("persistent filesystem boundary: PASS", flush=True)
    systemd_run("temporary", *common, "TemporaryFileSystem=/tmp:rw,size=16M,mode=1777")
    print("private temporary filesystem boundary: PASS", flush=True)
    systemd_run("detached", *common, wait=False)
    pid_file = ROOT / "detached.pid"
    for _ in range(100):
        if pid_file.exists():
            break
        time.sleep(0.1)
    assert pid_file.exists(), "detached child did not start"
    pid = int(pid_file.read_text())
    subprocess.run(["systemctl", "stop", "dvbfixer-accept-detached.service"], check=True)
    for _ in range(100):
        if not pathlib.Path(f"/proc/{pid}").exists():
            break
        time.sleep(0.1)
    assert not pathlib.Path(f"/proc/{pid}").exists(), "session-detached child survived service stop"
    pid_file.unlink()
    print("session-detached descendant teardown: PASS", flush=True)


if __name__ == "__main__":
    main()
