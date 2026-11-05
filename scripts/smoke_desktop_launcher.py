#!/usr/bin/env python3
"""Exercise a built, installed GUI without opening a browser or changing user data.

Requires an existing test prefix containing Python/package, Node, built GUI and
Node dependencies. This is launcher acceptance, not scientific-profile acceptance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import secrets
import signal
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New directory for isolated data and retained evidence")
    parser.add_argument("--shortcut", action="store_true", help="Also execute a shortcut generated under the isolated test home")
    args = parser.parse_args()
    prefix, output = args.prefix.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = {key: value for key, value in os.environ.items() if not key.startswith("DVBFIXER_") and key not in {"PYTHONPATH", "PYTHONHOME", "XDG_DATA_HOME"}}
    env["PATH"] = "/usr/bin:/bin"
    python = str(prefix / "bin/python")
    checks = []
    subprocess.run([python, "-m", "dvbfixer.cli", "gui", "--check"], cwd=output, env=env, check=True)
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    shortcut = None
    if args.shortcut:
        result = subprocess.run([python, "-c",
            "from pathlib import Path; import sys; from dvbfixer.desktop_shortcuts import install_desktop_shortcut; "
            "print(install_desktop_shortcut(Path(sys.argv[1]), home=Path(sys.argv[2])))",
            str(prefix), str(output / "shortcut-home")], cwd=output, env=env, check=True, capture_output=True, text=True)
        target = Path(result.stdout.strip())
        shortcut = target / "Contents/MacOS/DVBFixer" if target.suffix == ".app" else prefix / "share/dvbfixer/desktop-launcher"
    for name in ("local", "authenticated", *(["shortcut"] if shortcut else [])):
        authenticated = name == "authenticated"
        arm = output / name
        arm.mkdir()
        child_env = dict(env)
        token = secrets.token_urlsafe(32)
        if authenticated:
            child_env["DVBFIXER_AUTH_PRINCIPALS"] = json.dumps({"version": 1, "principals": [
                {"id": "smoke", "tokenSha256": hashlib.sha256(token.encode()).hexdigest()},
            ]})
            child_env["DVBFIXER_LEGACY_WORKSPACE_OWNER"] = "smoke"
        base = ([str(shortcut)] if name == "shortcut" else [python, "-m", "dvbfixer.cli", "gui"])
        base.extend(["--data-dir", str(arm / "workspaces")])
        stdout, stderr = arm / "stdout.log", arm / "stderr.log"
        process = None
        try:
            with stdout.open("w") as out, stderr.open("w") as err:
                process = subprocess.Popen([*base, "--no-browser"], cwd=output, env=child_env, stdout=out, stderr=err)
                deadline = time.monotonic() + 40
                url = ""
                while time.monotonic() < deadline:
                    if process.poll() is not None:
                        raise RuntimeError(f"launcher exited; inspect {stderr}")
                    match = re.search(r"^DVBFixer GUI: (http://127\.0\.0\.1:\d+)$", stdout.read_text(), re.MULTILINE)
                    if match:
                        url = match[1]
                        break
                    time.sleep(0.05)
                if not url:
                    raise RuntimeError("launcher readiness timed out")
                if authenticated:
                    try:
                        opener.open(url + "/api/session", timeout=5)
                    except urllib.error.HTTPError as error:
                        assert error.code == 401
                    else:
                        raise AssertionError("anonymous API access was not rejected")
                request = urllib.request.Request(url + "/api/health", headers={"Authorization": f"Bearer {token}"} if authenticated else {})
                with opener.open(request, timeout=5) as response:
                    health = json.load(response)
                    assert health["status"] == "ready"
                    assert health["authenticationRequired"] is authenticated
                session = urllib.request.Request(url + "/api/session", headers={"Authorization": f"Bearer {token}"} if authenticated else {})
                with opener.open(session, timeout=5) as response:
                    assert json.load(response)["principal"]["id"] == ("smoke" if authenticated else "local")
                repeated = subprocess.run([*base, "--no-browser"], cwd=output, env=child_env, text=True, capture_output=True, check=True, timeout=15)
                assert f"DVBFixer GUI: {url}" in repeated.stdout
                if authenticated:
                    process.send_signal(signal.SIGTERM)
                else:
                    subprocess.run([*base, "--stop"], cwd=output, env=child_env, check=True, timeout=15)
                process.wait(timeout=25)
                assert process.returncode == 0
            assert token not in stdout.read_text() + stderr.read_text()
            checks.append({"arm": arm.name, "ready": True, "reused": True, "stopped": True,
                           "auth_preserved": authenticated, "browser_opened": False})
        finally:
            if process is not None and process.poll() is None:
                process.terminate()
                process.wait(timeout=25)
    (output / "summary.json").write_text(json.dumps({"kind": "installed-gui-launcher-smoke", "checks": checks,
        "scientific_profile_acceptance": False, "native_shortcut_acceptance": False}, indent=2) + "\n")
    print(f"Installed GUI launcher smoke passed: {output / 'summary.json'}")


if __name__ == "__main__":
    main()
