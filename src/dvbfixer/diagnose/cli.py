"""CLI parsing for ``dvbfixer diagnose``.

Argparse only — no OpenMM / PDBFixer imports at module scope.
"""

from __future__ import annotations

import argparse


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="dvbfixer diagnose",
        description="Inspect a PDB file and report structure-quality "
                    "issues (missing atoms, coincident atoms, valence "
                    "violations, steric clashes, chain breaks, etc.). "
                    "Report-only — never mutates the input. Use "
                    "`dvbfixer prepare` / `minimize` / `protonate` / "
                    "`pull` to apply repairs.",
    )

    io = p.add_argument_group("Input / output")
    io.add_argument("input", help="Input PDB file")
    io.add_argument("-o", "--output",
                    help="Write report to file (default: stdout)")

    checks = p.add_argument_group("Check selection")
    checks.add_argument(
        "--only",
        choices=["all", "structural", "chemistry", "steric"],
        default="all",
        help="Restrict to one category of checks (default: all).",
    )
    checks.add_argument(
        "--severity",
        choices=["ERROR", "WARNING", "INFO"],
        default="INFO",
        help="Minimum severity to include in the report (default: INFO). "
             "Set to ERROR for pre-commit / CI usage.",
    )

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument(
        "-v", "--verbose", action="store_true",
        help="Include per-check timing and expanded per-atom detail.",
    )

    return p.parse_args(argv)
