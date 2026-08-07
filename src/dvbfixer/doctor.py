"""Installation and external-tool capability report."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import shutil
import subprocess
import sys
from typing import Any

PYTHON_PACKAGES = {
    "OpenMM": "openmm", "PDBFixer": "pdbfixer", "Modeller": "modeller",
    "MDAnalysis": "MDAnalysis", "PROPKA": "propka", "ParmEd": "parmed",
    "ACPYPE": "acpype", "Open Babel Python": "openbabel", "RDKit": "rdkit",
}
EXECUTABLES = {
    "antechamber": "antechamber", "parmchk2": "parmchk2", "tleap": "tleap",
    "Reduce": "reduce", "Probe": "probe", "Open Babel": "obabel",
    "obminimize": "obminimize", "xTB": "xtb", "GROMACS": "gmx",
    "MAFFT": "mafft", "MUSCLE": "muscle", "Clustal Omega": "clustalo",
}


def _package_status(module: str) -> dict[str, Any]:
    try:
        spec = importlib.util.find_spec(module)
    except (ImportError, ModuleNotFoundError, ValueError):
        spec = None
    status: dict[str, Any] = {"available": spec is not None}
    if spec is not None:
        distribution = {"openmm": "OpenMM", "MDAnalysis": "MDAnalysis"}.get(module, module)
        try:
            status["version"] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            status["version"] = "unknown"
    return status


def collect_capabilities() -> dict[str, Any]:
    """Collect capabilities without importing heavyweight optional packages."""
    packages = {name: _package_status(module) for name, module in PYTHON_PACKAGES.items()}
    if packages["Modeller"]["available"]:
        try:
            probe = subprocess.run(
                [sys.executable, "-c", "import modeller"],
                capture_output=True, text=True, timeout=10, check=False,
            )
            packages["Modeller"]["usable"] = probe.returncode == 0
            if probe.returncode:
                message = (probe.stderr or probe.stdout).strip().splitlines()
                packages["Modeller"]["error"] = message[0] if message else "import failed"
        except (OSError, subprocess.TimeoutExpired) as exc:
            packages["Modeller"]["usable"] = False
            packages["Modeller"]["error"] = f"{type(exc).__name__}: {exc}"
    executables = {
        name: {"available": (path := shutil.which(executable)) is not None, "path": path}
        for name, executable in EXECUTABLES.items()
    }
    platforms: list[str] = []
    if packages["OpenMM"]["available"]:
        try:
            import openmm
            platforms = [
                openmm.Platform.getPlatform(index).getName()
                for index in range(openmm.Platform.getNumPlatforms())
            ]
        except Exception as exc:
            platforms = [f"unavailable: {type(exc).__name__}: {exc}"]
    return {"python_packages": packages, "executables": executables, "openmm_platforms": platforms}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer doctor",
        description="Report optional Python packages, external executables, and OpenMM platforms.",
    )
    parser.add_argument("--format", choices=["text", "json"], default="text")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    report = collect_capabilities()
    if args.format == "json":
        print(json.dumps(report, indent=2))
        return
    print("Python packages:")
    for name, status in report["python_packages"].items():
        usable = status["available"] and status.get("usable", True)
        detail = f" ({status.get('version')})" if status["available"] else ""
        if status.get("error"):
            detail += f": {status['error']}"
        print(f"  {'OK' if usable else 'MISSING':7s} {name}{detail}")
    print("External executables:")
    for name, status in report["executables"].items():
        detail = f" ({status['path']})" if status["available"] else ""
        print(f"  {'OK' if status['available'] else 'MISSING':7s} {name}{detail}")
    platforms = ", ".join(report["openmm_platforms"]) or "none"
    print(f"OpenMM platforms: {platforms}")


if __name__ == "__main__":
    main()
