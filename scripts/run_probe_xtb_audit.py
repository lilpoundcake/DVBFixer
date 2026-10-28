"""Portable macOS/Linux replay of the bounded Probe/xTB audit protocol.

Run in the scientific environment, providing the independently built binaries:
  python scripts/run_probe_xtb_audit.py --output .artifacts/audit-linux \
    --probe /path/to/probe --xtb installed=/path/to/xtb --xtb fixed=/path/to/fixed-xtb

No packages are installed. Outputs must not exist. Default mode prepares fresh
baselines from checksummed tracked inputs. --baselines replays existing coordinates
and is explicitly NOT an independent preparation/platform qualification.
Audit completion does not imply scientific acceptance: inspect report.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "hinge": "tests/fixtures/hinge_CH3_glycosylated.pdb",
    "7x35": "tests/fixtures/lipid/7x35_r_u.pdb",
    "ubq": "tests/fixtures/protein_context/1UBQ.pdb",
    "asn": "tests/fixtures/ASN.pdb",
}


def verified_sources():
    inventory = {}
    for line in (ROOT / "tests/fixtures/MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        inventory[name] = digest
    result = {}
    for name in SOURCES.values():
        digest = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if digest != inventory[name]:
            raise ValueError(f"Source fixture checksum mismatch: {name}")
        result[name] = digest
    return result


def binary_bindings(values):
    result = {}
    for value in values:
        label, separator, path = value.partition("=")
        if not separator or not re.fullmatch(r"[A-Za-z0-9_-]+", label) or label in result:
            raise ValueError("--xtb requires unique safe labels: LABEL=/absolute/path")
        executable = Path(path).resolve()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f"Executable is unavailable: {executable}")
        result[label] = str(executable)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--xtb", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--baselines", type=Path, help="Replay hinge-openmm.pdb, 7x35-openmm.pdb, ubq-prepared.pdb")
    parser.add_argument("--extended-analysis", action="store_true", help="Also compare separate glycan components and measure geometry deltas")
    args = parser.parse_args()
    try:
        bindings, sources = binary_bindings(args.xtb), verified_sources()
    except ValueError as exc:
        parser.error(str(exc))
    args.probe = args.probe.resolve()
    if not args.probe.is_file() or not os.access(args.probe, os.X_OK):
        parser.error("Probe executable is unavailable")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    baseline_dir = args.baselines.resolve() if args.baselines else args.output / "baselines"
    if not args.baselines:
        baseline_dir.mkdir()
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    report = {"protocol": "probe-xtb-macos-linux-v2-removed-production-xtb", "platform": platform.platform(),
              "production_xtb_adapter": "removed; native research experiments remain independent",
              "python": sys.version, "source_sha256": sources, "tasks": [],
              "baseline_mode": "frozen coordinate replay" if args.baselines else "fresh preparation and OpenMM",
              "binary_sha256": {k: hashlib.sha256(Path(v).read_bytes()).hexdigest() for k, v in bindings.items()},
              "probe_sha256": hashlib.sha256(args.probe.read_bytes()).hexdigest(),
              "acceptance_note": "Completed execution is not an assertion of chemical quality or adapter safety"}
    scripts = ["run_probe_xtb_audit.py", "audit_probe_contacts.py", "audit_xtb_controls.py", "audit_probe_xtb_contracts.py",
               "audit_xtb_geometry_deltas.py"]
    snapshot = args.output / "scripts"
    snapshot.mkdir()
    for name in scripts:
        (snapshot / name).write_bytes((ROOT / "scripts" / name).read_bytes())
    guard = args.output / "no-xtb-discovery"
    guard.mkdir()
    marker = guard / "unexpected-lookup.txt"
    (guard / "sitecustomize.py").write_text(
        "import shutil\nfrom pathlib import Path\n_original = shutil.which\n"
        "def _guard(command, *args, **kwargs):\n"
        "    if command == 'xtb':\n"
        f"        Path({str(marker)!r}).write_text('unexpected xtb lookup')\n"
        "        raise RuntimeError('xTB discovery without an explicit refinement request')\n"
        "    return _original(command, *args, **kwargs)\nshutil.which = _guard\n"
    )

    def run(label, command, timeout=300, guard_discovery=False):
        command = [str(x) for x in command]
        started = time.monotonic()
        row = {"label": label, "command": command}
        try:
            child_env = dict(env)
            if guard_discovery:
                child_env["PYTHONPATH"] = str(guard) + os.pathsep + child_env.get("PYTHONPATH", "")
            result = subprocess.run(command, cwd=ROOT, env=child_env, capture_output=True, timeout=timeout)
            (args.output / f"{label}.stdout").write_bytes(result.stdout)
            (args.output / f"{label}.stderr").write_bytes(result.stderr)
            row["returncode"] = result.returncode
        except subprocess.TimeoutExpired as exc:
            (args.output / f"{label}.stdout").write_bytes(exc.stdout or b"")
            (args.output / f"{label}.stderr").write_bytes(exc.stderr or b"")
            row["timeout"] = True
        row["seconds"] = time.monotonic() - started
        report["tasks"].append(row)
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(label, row.get("returncode", "timeout"), flush=True)
        return row.get("returncode") == 0

    cli = Path(sys.executable).with_name("dvbfixer")
    if not args.baselines:
        for label in ("hinge", "7x35", "ubq"):
            prepare = [cli, "prepare", ROOT / SOURCES[label], "--atom-naming", "standard",
                       "-o", baseline_dir / f"{label}-prepared.pdb", "--dat", baseline_dir / f"{label}.dat"]
            if label == "7x35":
                prepare += ["--ph", "7.4", "--smiles", "PLM=CCCCCCCCCCCCCCCC(=O)[O-]"]
            if run(f"prepare-{label}", prepare, guard_discovery=True) and label != "ubq":
                minimize = [cli, "minimize", baseline_dir / f"{label}-prepared.pdb", "--dat", baseline_dir / f"{label}.dat",
                            "--no-solvent", "--max-iter", "200", "--platform", "CPU", "--atom-naming", "standard",
                            "-o", baseline_dir / f"{label}-openmm.pdb"]
                if label == "7x35":
                    minimize.append("--parametrize-ligands")
                run(f"minimize-{label}", minimize, guard_discovery=True)
    baseline_paths = [baseline_dir / name for name in ("hinge-openmm.pdb", "7x35-openmm.pdb", "ubq-prepared.pdb")]
    report["baseline_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in baseline_paths if p.is_file()}
    missing = [str(p) for p in baseline_paths if not p.is_file()]
    report["missing_baselines"] = missing
    report["no_xtb_discovery_verified"] = (
        not args.baselines and not missing and not marker.exists()
        and all(row.get("returncode") == 0 for row in report["tasks"])
    )
    run("native-contracts", [sys.executable, snapshot / "audit_probe_xtb_contracts.py", "--output", args.output / "native-contracts"])
    # Execute repository scripts (not snapshot copies, whose ROOT would differ).
    for radii in ("electron", "nuclear"):
        command = [sys.executable, ROOT / "scripts/audit_probe_contacts.py", "--probe", args.probe,
                   "--output", args.output / f"probe-{radii}", "--align-atom-fields", "--input", ROOT / SOURCES["asn"]]
        if radii == "nuclear":
            command.append("--probe-nuclear")
        for path in baseline_paths:
            if path.is_file():
                command += ["--input", path]
        run(f"contacts-{radii}", command)
    for label, binary in bindings.items():
        common = [sys.executable, ROOT / "scripts/audit_xtb_controls.py", "--xtb", binary, "--timeout", "90"]
        run(f"{label}-controls", [*common, "--output", args.output / f"{label}-controls",
                                  "--fragmented-control", "--engines", "default", "rf"], timeout=600)
        if baseline_paths[0].is_file():
            for arm, cycles in (("production", "200"), ("extended", "1000")):
                run(f"{label}-glycan-{arm}", [*common, "--output", args.output / f"{label}-glycan-{arm}",
                    "--pdb", baseline_paths[0], "--glycam-bonds", "--case", "pdb-0-hinge-openmm",
                    "--engines", "default", "rf", "--cycles", cycles], timeout=240)
            if args.extended_analysis:
                for component in (0, 1):
                    run(f"{label}-glycan-separate-{component}", [*common, "--output", args.output / f"{label}-glycan-separate-{component}",
                        "--pdb", baseline_paths[0], "--glycam-bonds", "--case", "pdb-0-hinge-openmm",
                        "--component", str(component), "--engines", "default", "--cycles", "1000"], timeout=240)
        if baseline_paths[1].is_file():
            for arm, flags in [("production", []), ("charge", ["--charge", "-1", "--uhf", "0"]),
                               ("water", ["--charge", "-1", "--uhf", "0", "--solvent", "water"]),
                               ("qm-vacuum", ["--method", "gfn2", "--charge", "-1", "--uhf", "0"]),
                               ("qm-water", ["--method", "gfn2", "--charge", "-1", "--uhf", "0", "--solvent", "water"])]:
                run(f"{label}-plm-{arm}", [*common, "--output", args.output / f"{label}-plm-{arm}",
                    "--pdb", baseline_paths[1], "--case", "pdb-0-7x35-openmm", "--engines", "default", *flags])
    if args.extended_analysis:
        for directory in sorted(args.output.iterdir()):
            if not directory.is_dir() or not (directory / "summary.json").exists():
                continue
            glycan = "-glycan-" in directory.name
            if not glycan and "-plm-" not in directory.name:
                continue
            run(f"geometry-{directory.name}", [sys.executable, ROOT / "scripts/audit_xtb_geometry_deltas.py",
                "--baseline", baseline_paths[0 if glycan else 1], "--run", directory,
                "--output", args.output / f"geometry-{directory.name}.json", *(["--glycam-bonds"] if glycan else [])])
    report["execution_complete"] = not missing and all(task.get("returncode") == 0 for task in report["tasks"])
    contract_report = args.output / "native-contracts/report.json"
    report["adapter_contracts_met"] = (
        json.loads(contract_report.read_text())["all_requirements_met"] if contract_report.is_file() else None
    )
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = {str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if not report["execution_complete"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
