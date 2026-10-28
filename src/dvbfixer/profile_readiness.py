"""Executable profile prerequisites, isolated from the reporting Python process."""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from dvbfixer.installation_profiles import composition


def collect_profile_readiness(profile: str, addons: tuple[str, ...] = ()) -> dict[str, Any]:
    checks: list[dict[str, str]] = []

    def record(name: str, state: str, action: str = "") -> None:
        checks.append({"name": name, "state": state, "action": action})

    record("python-version", "ready" if (3, 11) <= sys.version_info[:2] < (3, 14) else "error",
           "Use Python >=3.11,<3.14.")
    for package in composition(profile, addons):
        module = package.get("module")
        if module:
            # Actual import catches missing shared libraries and ABI failures. Never
            # relay raw MODELLER output: third-party diagnostics may include secrets.
            try:
                result = subprocess.run([sys.executable, "-c", f"import {module}"],
                                        capture_output=True, timeout=30, check=False)
                state = "ready" if result.returncode == 0 else "error"
            except (OSError, subprocess.TimeoutExpired):
                state = "error"
            if state == "error" and module == "modeller":
                try:
                    installed = importlib.util.find_spec("modeller") is not None
                except (ImportError, ValueError):
                    installed = False
                if installed:
                    state = "action-required"
            record(module, state, package.get("manual", f"Reinstall {package['spec']} in the selected environment."))
        for executable in package.get("executables", []):
            record(executable, "ready" if shutil.which(executable) else "error",
                   f"Install {package['spec']} and put this environment's bin directory on PATH.")
    try:
        result = subprocess.run([sys.executable, "-c", "from pathlib import Path; "
                                "from anarci.anarci import HMM_path; "
                                "assert all(Path(HMM_path, 'ALL.hmm' + suffix).is_file() "
                                "for suffix in ('', '.h3f', '.h3i', '.h3m', '.h3p'))"],
                                capture_output=True, timeout=30, check=False)
        record("anarci-data", "ready" if result.returncode == 0 else "error",
               "Reinstall ANARCI with its pressed ALL.hmm database and HMMER.")
    except (OSError, subprocess.TimeoutExpired):
        record("anarci-data", "error", "ANARCI database check failed or timed out.")
    if "probe" in addons:
        from dvbfixer.diagnose.probe import SUPPORTED_VERSION

        try:
            result = subprocess.run(["probe", "-version"], capture_output=True,
                                    text=True, timeout=10, check=False)
            supported = result.returncode == 0 and result.stdout.strip() == SUPPORTED_VERSION
        except (OSError, subprocess.TimeoutExpired):
            supported = False
        record("probe-version", "ready" if supported else "error",
               f"Install the supported MolProbity {SUPPORTED_VERSION}; diagnose otherwise uses Python fallback.")
    # The catalog alone cannot prohibit unwanted *transitive* installs.
    conda_names = {p.name.split("-")[0] for p in Path(sys.prefix, "conda-meta").glob("*.json")}
    ml = any(importlib.util.find_spec(name) is not None for name in ("torch", "protenix"))
    record("runtime-isolation", "error" if ml else "ready",
           "Keep PyTorch/diffusion in a separate runtime; use openff-toolkit-base.")
    record("openbabel-provider", "error" if "openbabel" in conda_names else "ready",
           "Use only openbabel-wheel in this profile; do not overlay conda openbabel.")
    code = """import numpy as np
from scipy.linalg import solve
a = np.array([[3., 1.], [1., 2.]])
b = np.array([9., 8.])
assert np.allclose(np.linalg.solve(a,b), [2.,3.])
assert np.allclose(solve(a,b), [2.,3.])
"""
    try:
        result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                                timeout=30, check=False)
        record("native-linalg", "ready" if result.returncode == 0 else "error",
               "Recreate the environment with all three OpenBLAS wrappers from the profile.")
    except (OSError, subprocess.TimeoutExpired):
        record("native-linalg", "error", "NumPy/SciPy native solve failed or timed out.")
    if profile == "desktop":
        gui = Path(sys.prefix, "share/dvbfixer/gui")
        record("gui-build", "ready" if all((gui / p).is_file() for p in
               ("dist/index.html", "dist-server/server.js")) else "error",
               "Run the Desktop helper to install/build the GUI from its lockfile.")
    states = {check["state"] for check in checks}
    state = "error" if "error" in states else "action-required" if "action-required" in states else "ready"
    return {"profile": profile, "addons": list(addons), "state": state,
            "scope": "installation prerequisites; scientific acceptance is recorded separately",
            "checks": checks}


def main() -> None:
    print(json.dumps(collect_profile_readiness(sys.argv[1]), indent=2))


if __name__ == "__main__":
    main()
