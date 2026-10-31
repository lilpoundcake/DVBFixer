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


def pdbfixer_resources_ready(python: str | None = None) -> bool:
    """Importability alone does not prove that heavy-atom rebuild assets exist."""
    from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

    code = """import json,sys
from pathlib import Path
import pdbfixer
from openmm.app import ForceField
root = Path(pdbfixer.__file__).parent
soft = root / 'soft.xml'
assert soft.is_file() and soft.stat().st_size > 0
ForceField(str(soft))
for name in json.loads(sys.argv[1]):
    template = root / 'templates' / (name + '.pdb')
    assert template.is_file() and template.stat().st_size > 0
"""
    try:
        result = subprocess.run([python or sys.executable, "-c", code, json.dumps(sorted(CANONICAL_HEAVY_ATOMS))],
                                capture_output=True, timeout=30, check=False)
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def collect_profile_readiness(profile: str, addons: tuple[str, ...] = ()) -> dict[str, Any]:
    checks: list[dict[str, str]] = []

    def record(name: str, state: str, action: str = "") -> None:
        checks.append({"name": name, "state": state, "action": action})

    record("python-version", "ready" if (3, 11) <= sys.version_info[:2] < (3, 14) else "error",
           "Use Python >=3.11,<3.14.")
    if "ambertools-compat" in addons:
        specs = [package["spec"] for package in composition(profile, addons)
                 if package["group"] == "ambertools-compat"]
        script = '''import importlib.metadata,json,sys
from packaging.requirements import Requirement
for spec in json.loads(sys.argv[1]):
    requirement = Requirement(spec)
    installed = sys.version.split()[0] if requirement.name == "python" else importlib.metadata.version(requirement.name)
    if installed not in requirement.specifier:
        raise SystemExit(1)
'''
        try:
            result = subprocess.run([sys.executable, "-c", script, json.dumps(specs)],
                                    capture_output=True, timeout=30, check=False)
            compatible = result.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            compatible = False
        record("ambertools-compat-versions", "ready" if compatible else "error",
               "Create a new prefix with --addon ambertools-compat; do not downgrade an existing environment in place.")
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
            if module == "pdbfixer":
                record("pdbfixer-resources", "ready" if state == "ready" and pdbfixer_resources_ready() else "error",
                       "Recreate/reinstall the environment with complete PDBFixer soft.xml and canonical templates; do not patch missing assets piecemeal.")
        for executable in package.get("executables", []):
            record(executable, "ready" if shutil.which(executable) else "error",
                   f"Install {package['spec']} and put this environment's bin directory on PATH.")
    try:
        result = subprocess.run([sys.executable, "-c", "from dvbfixer.antibody import require_anarci; require_anarci()"],
                                capture_output=True, timeout=30, check=False)
        record("anarci-data", "ready" if result.returncode == 0 else "error",
               "Reinstall ANARCI with its pressed ALL.hmm database and HMMER.")
    except (OSError, subprocess.TimeoutExpired):
        record("anarci-data", "error", "ANARCI database check failed or timed out.")
    if "probe" in addons:
        from dvbfixer.diagnose.probe import SUPPORTED_VERSION

        try:
            probe_result = subprocess.run(["probe", "-version"], capture_output=True,
                                     text=True, timeout=10, check=False)
            supported = probe_result.returncode == 0 and probe_result.stdout.strip() == SUPPORTED_VERSION
        except (OSError, subprocess.TimeoutExpired):
            supported = False
        record("probe-version", "ready" if supported else "error",
               f"Install the supported MolProbity {SUPPORTED_VERSION}; diagnose otherwise uses Python fallback.")
    if "msa-alternatives" in addons:
        try:
            result = subprocess.run([sys.executable, "-c", "import shutil; "
                                    "from dvbfixer.msa import require_engine_version; "
                                    "require_engine_version('muscle', shutil.which('muscle') or 'muscle')"],
                                    capture_output=True, timeout=15, check=False)
            supported = result.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            supported = False
        record("muscle-version", "ready" if supported else "error", "Install MUSCLE 5 (muscle>=5,<6).")
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
        try:
            node_result = subprocess.run(["node", "--version"], capture_output=True, text=True,
                                    timeout=10, check=False)
            supported_node = node_result.returncode == 0 and int(node_result.stdout.strip().removeprefix("v").split(".")[0]) >= 22
        except (OSError, ValueError, subprocess.TimeoutExpired):
            supported_node = False
        record("node-version", "ready" if supported_node else "error", "Desktop requires Node.js >=22.")
        gui_roots = (Path(sys.prefix, "share/dvbfixer/gui"), Path(__file__).resolve().parents[2] / "gui")
        built = any(all((gui / p).is_file() for p in ("dist/index.html", "dist-server/server.js"))
                    for gui in gui_roots)
        record("gui-build", "ready" if built else "error",
               "Build Desktop using the helper or npm ci && npm run build in the editable source checkout's gui/.")
    states = {check["state"] for check in checks}
    state = "error" if "error" in states else "action-required" if "action-required" in states else "ready"
    return {"profile": profile, "addons": list(addons), "state": state,
            "scope": "installation prerequisites; scientific acceptance is recorded separately",
            "checks": checks}


def main() -> None:
    print(json.dumps(collect_profile_readiness(sys.argv[1]), indent=2))


if __name__ == "__main__":
    main()
