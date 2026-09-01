"""CLI parsing for ``dvbfixer diagnose``.

Argparse only — no OpenMM / PDBFixer imports at module scope.
"""

from __future__ import annotations

import argparse

from dvbfixer.diagnose.steric import CLASH_MODE_PRESETS, DEFAULT_CLASH_MODE


def _parse_clash_cutoff(s: str) -> tuple[float, float]:
    """Parse ``--clash-cutoff WARN,ERROR`` into a (warn_a, error_a) pair.

    Raises ``argparse.ArgumentTypeError`` on malformed input or when
    the warn threshold isn't ≤ the error threshold.
    """
    parts = s.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"expected 'WARN,ERROR' in Å (e.g. '0.4,0.5'), got {s!r}"
        )
    try:
        warn = float(parts[0])
        err = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"could not parse floats from {s!r}: {exc}"
        ) from exc
    if warn <= 0 or err <= 0:
        raise argparse.ArgumentTypeError(
            f"cutoffs must be positive (got warn={warn}, err={err})"
        )
    if warn > err:
        raise argparse.ArgumentTypeError(
            f"WARN cutoff ({warn}) must be ≤ ERROR cutoff ({err})"
        )
    return (warn, err)


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
    io.add_argument("input", help="Input PDB, PDBx/mmCIF, or crystallographic CIF file")
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
    checks.add_argument(
        "--include-water", action="store_true",
        help="Include water residues (HOH/WAT/TIP3/SOL) in checks. "
             "Off by default — crystallographic waters generate massive "
             "chain-break and steric noise.",
    )
    _presets = ", ".join(
        f"{name} ({warn}/{err} Å)"
        for name, (warn, err) in CLASH_MODE_PRESETS.items()
    )
    checks.add_argument(
        "--clash-mode",
        choices=sorted(CLASH_MODE_PRESETS.keys()),
        default=DEFAULT_CLASH_MODE,
        help=f"Preset clash overlap thresholds (WARN / ERROR in Å). "
             f"Default: {DEFAULT_CLASH_MODE}. Available: {_presets}.",
    )
    checks.add_argument(
        "--clash-cutoff", type=_parse_clash_cutoff, default=None,
        metavar="WARN,ERROR",
        help="Explicit clash overlap cutoffs in Å (overrides --clash-mode). "
             "Example: --clash-cutoff 0.35,0.45 for extra-strict validation.",
    )

    fmt = p.add_argument_group("Output format")
    fmt.add_argument(
        "--format", choices=["text", "json"], default="text",
        dest="output_format",
        help="Output format. `text` is the plain-text report; `json` "
             "emits a machine-readable list of findings (usable for CI).",
    )

    diag = p.add_argument_group("Diagnostics")
    diag.add_argument(
        "-v", "--verbose", action="store_true",
        help="Include per-check timing on stderr.",
    )

    from dvbfixer.batch import add_runtime_help
    add_runtime_help(p, batch=True)
    return p.parse_args(argv)
