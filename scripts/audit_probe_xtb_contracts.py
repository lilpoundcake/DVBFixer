"""Exercise current Probe contracts with bounded native fixtures.

The executables generated here are protocol fixtures, NOT scientific engines.
Only discovery and subprocess timeout are overridden; subprocess outputs are real.
The deleted production xTB adapter is recorded as removed, not reimplemented.
Historical xTB failures remain in commit 3b4ac0b and the archived audit evidence.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from openmm import Vec3
from openmm.app import Element, PDBFile, Topology
from openmm.unit import Quantity, nanometer

from dvbfixer.diagnose import steric
from dvbfixer.diagnose.probe import SUPPORTED_VERSION


def executable_body(engine, mode):
    body = "import sys,time\nfrom pathlib import Path\n"
    version = "probe.9.99" if mode == "unsupported-version" else SUPPORTED_VERSION
    body += f"if '-version' in sys.argv:\n    print({version!r})\n    sys.exit(0)\n"
    if mode == "nonzero":
        return body + "print('declared failure',file=sys.stderr)\nsys.exit(7)\n"
    if mode == "timeout":
        return body + "time.sleep(5)\n"
    if engine == "probe":
        output = {
            "clean": "",
            "invalid-format": "unsupported successful output\n",
            "valid-clash": ":1->1:bo: A   1 LIG  C1  : a   2 LIG  C2  :11:-1.400:-0.8:0:0:0:0:0:C:C:0:0:0:0:0\n",
            "nonfinite": ":1->1:bo: A   1 LIG  C1  : a   2 LIG  C2  :11:nan:-0.8:0:0:0:0:0:C:C:0:0:0:0:0\n",
        }.get(mode, "")
        return body + f"sys.stdout.write({output!r})\n"
    raise ValueError("The production xTB adapter was removed; use archived evidence for its historical contracts")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "audit-script.py").write_bytes(Path(__file__).read_bytes())
    top = Topology()
    for chain_id, resid, name in [("A", "1", "C1"), ("a", "2", "C2")]:
        residue = top.addResidue("LIG", top.addChain(chain_id), resid)
        top.addAtom(name, Element.getBySymbol("C"), residue)
    positions = Quantity([Vec3(0, 0, 0), Vec3(.2, 0, 0)], nanometer)
    input_path = args.output / "protocol-input.pdb"
    with input_path.open("w") as handle:
        PDBFile.writeFile(top, positions, handle, keepIds=True)
    real_run = subprocess.run
    real_python = steric.clashes_python
    report = {"scope": "native subprocess protocol fixtures, not scientific acceptance",
              "timeout_override_seconds_for_timeout_fixture": .5,
              "removed_adapters": {"xtb": {"status": "removed", "historical_commit": "3b4ac0b"}}, "runs": []}
    for engine, modes in [
        ("probe", ["missing", "nonzero", "timeout", "invalid-format", "nonfinite", "unsupported-version", "clean", "valid-clash"]),
    ]:
        for mode in modes:
            work = args.output / f"{engine}-{mode}"
            work.mkdir()
            binary = work / engine
            binary.write_text(f"#!{sys.executable}\n" + executable_body(engine, mode))
            binary.chmod(0o700)
            events, fallback = [], []

            def bounded_run(*a, _mode=mode, _events=events, **kw):
                kw["timeout"] = .5 if _mode == "timeout" else 30
                try:
                    result = real_run(*a, **kw)
                except subprocess.TimeoutExpired:
                    _events.append({"timeout": True})
                    raise
                _events.append({"command": result.args, "returncode": result.returncode,
                               "stdout": result.stdout, "stderr": result.stderr})
                return result

            def python_engine(*a, _fallback=fallback, **kw):
                _fallback.append(True)
                return real_python(*a, **kw)

            module = steric
            log = io.StringIO()
            record = {"engine": engine, "mode": mode, "native_events": events}
            with patch.object(module, "_find_binary", return_value=None if mode == "missing" else str(binary.resolve())), \
                    patch.object(subprocess, "run", bounded_run), \
                    patch.object(steric, "clashes_python", python_engine), \
                    contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                try:
                    provenance = {}
                    result = steric.run_all(top, positions, input_path, clash_warn_a=.4, provenance=provenance)
                    record.update({"python_fallback": bool(fallback), "findings": len(result), "provenance": provenance})
                    record["requirement_met"] = (
                        bool(fallback) if mode not in ("clean", "valid-clash")
                        else not fallback and (len(result) == 0 if mode == "clean" else len(result) > 0)
                    )
                except Exception as exc:
                    record.update({"exception": type(exc).__name__, "message": str(exc), "requirement_met": False})
            (work / "adapter.stdout").write_text(log.getvalue())
            record["fixture_valid"] = (mode == "timeout") == any(event.get("timeout") for event in events)
            record["requirement_met"] = bool(record["requirement_met"] and record["fixture_valid"])
            report["runs"].append(record)
            print(engine, mode, "contract met", record["requirement_met"], flush=True)
    report["all_requirements_met"] = all(row["requirement_met"] for row in report["runs"])
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = {str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
