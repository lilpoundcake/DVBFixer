"""Exercise production adapters against real, bounded subprocess failure fixtures.

The executables generated here are protocol fixtures, NOT scientific engines.
Only discovery and subprocess timeout are overridden; subprocess outputs are real.
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

import numpy as np
from openmm import Vec3
from openmm.app import Element, Topology
from openmm.unit import Quantity, nanometer

from dvbfixer.diagnose import steric
from dvbfixer.minimize import refine


def executable_body(engine, mode):
    body = "import sys,time\nfrom pathlib import Path\n"
    if mode == "nonzero":
        return body + "print('declared failure',file=sys.stderr)\nsys.exit(7)\n"
    if mode == "timeout":
        return body + "time.sleep(5)\n"
    if engine == "probe":
        output = {
            "clean": "",
            "invalid-format": "unsupported successful output\n",
            "valid-clash": ":1->1:bo: A   1 LIG  C1  : a   2 LIG  C2  :11:-1.400:-0.8:0:0:0:0:0:C:C:0:0:0:0:0\n",
        }.get(mode, "")
        return body + f"sys.stdout.write({output!r})\n"
    xyz = {
        "valid": "2\ncontrol\nC 0 0 0\nC 2 0 0\n",
        "unconverged": "2\ncontrol\nC 0 0 0\nC 2 0 0\n",
        "bad-count": "1\ncontrol\nC 0 0 0\n",
        "bad-elements": "2\ncontrol\nO 0 0 0\nC 2 0 0\n",
        "nonfinite": "2\ncontrol\nC nan 0 0\nC 2 0 0\n",
        "truncated": "2\ncontrol\nC 0 0 0\n",
    }.get(mode)
    if xyz is not None:
        body += f"Path('xtbopt.xyz').write_text({xyz!r})\n"
    message = "FAILED TO CONVERGE" if mode == "unconverged" else "GEOMETRY OPTIMIZATION CONVERGED"
    return body + f"print({message!r})\n"


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
    real_run = subprocess.run
    real_python = steric.clashes_python
    report = {"scope": "native subprocess protocol fixtures, not scientific acceptance",
              "timeout_override_seconds_for_timeout_fixture": .5, "runs": []}
    for engine, modes in [
        ("probe", ["missing", "nonzero", "timeout", "invalid-format", "clean", "valid-clash"]),
        ("xtb", ["missing", "nonzero", "timeout", "no-output", "bad-count", "bad-elements",
                 "nonfinite", "truncated", "unconverged", "valid"]),
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

            module = steric if engine == "probe" else refine
            log = io.StringIO()
            record = {"engine": engine, "mode": mode, "native_events": events}
            with patch.object(module, "_find_binary", return_value=None if mode == "missing" else str(binary.resolve())), \
                    patch.object(subprocess, "run", bounded_run), \
                    patch.object(steric, "clashes_python", python_engine), contextlib.redirect_stdout(log):
                try:
                    if engine == "probe":
                        result = steric.run_all(top, positions, "unused.pdb", clash_warn_a=.4)
                        record.update({"python_fallback": bool(fallback), "findings": len(result)})
                        record["requirement_met"] = (
                            bool(fallback) if mode in ("missing", "nonzero", "timeout", "invalid-format")
                            else not fallback and (len(result) == 0 if mode == "clean" else len(result) > 0)
                        )
                    else:
                        result = refine.refine_with_xtb(top, positions, cycles=2)
                        array = np.asarray(result.value_in_unit(nanometer))
                        unchanged = result is positions
                        record.update({"skipped_preserving_input": unchanged, "atom_count": len(result),
                                       "coordinates_finite": bool(np.isfinite(array).all())})
                        record["requirement_met"] = (
                            len(result) == 2 and np.isfinite(array).all() and not unchanged
                            if mode == "valid" else unchanged
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
