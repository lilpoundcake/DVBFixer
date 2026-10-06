"""Installation and external-tool capability report."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.contract import DIFFUSION_SCHEMA_VERSION
from dvbfixer.model.diffusion.preflight import (
    DIFFUSION_PROFILES,
    AdapterPreflightCode,
    AdapterPreflightIssue,
    RunnerPreflightFacts,
    RunnerPreflightReport,
    invoke_runner_preflight,
)
from dvbfixer.model.diffusion.runner import DIFFUSION_RUNNER_PROTOCOL_VERSION
from dvbfixer.model.diffusion.runtime import resolve_diffusion_runtime

PYTHON_PACKAGES = {
    "OpenMM": "openmm", "PDBFixer": "pdbfixer", "Modeller": "modeller",
    "MDAnalysis": "MDAnalysis", "PROPKA": "propka", "ParmEd": "parmed",
    "ACPYPE": "acpype", "Open Babel Python": "openbabel", "RDKit": "rdkit",
    "Gemmi": "gemmi",
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
        distribution = {
            "openmm": "OpenMM", "MDAnalysis": "MDAnalysis",
            "openbabel": "openbabel-wheel",
        }.get(module, module)
        try:
            status["version"] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            status["version"] = "unknown"
    return status


def collect_capabilities(
    *,
    diffusion_model: str | None = None,
    diffusion_profile: str | None = None,
    diffusion_runner: str | None = None,
    diffusion_checkpoint: str | None = None,
    diffusion_timeout: float = 30.0,
) -> dict[str, Any]:
    """Collect capabilities without importing heavyweight optional packages."""
    if diffusion_model is not None:
        runtime = resolve_diffusion_runtime(diffusion_model)
        diffusion_profile = runtime.profile
        diffusion_runner = runtime.runner
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
    diffusion: dict[str, Any] = {
        "profiles": {
            name: metadata.to_dict() for name, metadata in DIFFUSION_PROFILES.items()
        },
        "selected_profile": None,
    }
    if diffusion_profile is not None:
        issues: list[AdapterPreflightIssue] = []
        runner_path = shutil.which(diffusion_runner) if diffusion_runner else None
        if diffusion_runner and (
            Path(diffusion_runner).is_absolute() or Path(diffusion_runner).parent != Path(".")
        ):
            candidate = Path(diffusion_runner).expanduser()
            try:
                executable = candidate.is_file() and bool(candidate.stat().st_mode & 0o111)
            except OSError:
                executable = False
            runner_path = str(candidate) if executable else None
        if not diffusion_runner or not runner_path:
            issues.append(AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_RUNNER,
                "the selected diffusion model has no installed executable launcher",
            ))
        checkpoint = Path(diffusion_checkpoint).expanduser() if diffusion_checkpoint else None
        if checkpoint is not None and (checkpoint.is_symlink() or not checkpoint.is_file()):
            issues.append(AdapterPreflightIssue(
                AdapterPreflightCode.MISSING_CHECKPOINT,
                "the diffusion checkpoint is not a regular file",
            ))
        if not issues:
            assert runner_path is not None
            preflight = invoke_runner_preflight(
                profile=diffusion_profile,
                runner=runner_path,
                checkpoint=checkpoint,
                timeout_seconds=diffusion_timeout,
            )
        else:
            preflight = RunnerPreflightReport(
                profile=diffusion_profile,
                runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
                contract_schema_version=DIFFUSION_SCHEMA_VERSION,
                facts=RunnerPreflightFacts(),
                issues=tuple(issues),
            )
        diffusion["selected_profile"] = preflight.to_dict()
    return {
        "python_packages": packages,
        "executables": executables,
        "openmm_platforms": platforms,
        "diffusion": diffusion,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer doctor",
        description="Report optional packages, executables, OpenMM platforms, and diffusion profiles.",
    )
    parser.add_argument("--format", choices=["text", "json"], default="text",
                        help="Report format (default: text)")
    diffusion = parser.add_argument_group("Diffusion options")
    diffusion.add_argument(
        "--diffusion-model",
        choices=("protpardelle", "protenix"),
        help="Optionally preflight an installed diffusion model",
    )
    diffusion.add_argument(
        "--diffusion-timeout",
        type=float,
        default=30.0,
        metavar="SECONDS",
        help="Preflight handshake timeout in seconds (default: 30)",
    )
    from dvbfixer.batch import add_runtime_help
    add_runtime_help(parser)
    args = parser.parse_args(argv)
    if not math.isfinite(args.diffusion_timeout) or args.diffusion_timeout <= 0:
        parser.error("--diffusion-timeout must be a positive finite number")
    return args


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    report = collect_capabilities(
        diffusion_model=args.diffusion_model,
        diffusion_timeout=args.diffusion_timeout,
    )
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
    print("Diffusion profiles:")
    for name, metadata in report["diffusion"]["profiles"].items():
        labels = ", ".join(metadata["evidence_labels"])
        print(f"  {name}: {metadata['status']} [{labels}]")
    selected = report["diffusion"]["selected_profile"]
    if selected is not None:
        state = "OK" if selected["passed"] else "BLOCKED"
        print(f"Selected diffusion profile: {selected['profile']} ({state})")
        for issue in selected["issues"]:
            print(f"  {issue['code']}: {issue['message']}")


if __name__ == "__main__":
    main()
