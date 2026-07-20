"""Orchestrator for ``dvbfixer diagnose``.

Loads the input PDB once, prepares a shared ``PDBFixer`` for the
structural checks that need it, dispatches to the three check
families, filters findings by severity, formats the report, and
exits with code 1 if any ERROR finding remains.
"""

from __future__ import annotations

import sys
from pathlib import Path

from dvbfixer.diagnose import chemistry, steric, structural
from dvbfixer.diagnose.cli import parse_args
from dvbfixer.diagnose.report import (
    Finding,
    Severity,
    format_report,
    sev_at_least,
)


def _load_topology_and_fixer(pdb_path: Path):
    """Load the input via OpenMM ``PDBFile`` (for topology + positions)
    and via ``PDBFixer`` (for missing-atom / missing-terminals checks).
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    pdb = PDBFile(str(pdb_path))
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    return pdb.topology, pdb.positions, fixer


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(2)

    topology, positions, fixer = _load_topology_and_fixer(input_path)

    findings: list[Finding] = []
    if args.only in ("all", "structural"):
        findings.extend(structural.run_all(topology, positions, input_path, fixer))
    if args.only in ("all", "chemistry"):
        findings.extend(chemistry.run_all(topology, positions))
    if args.only in ("all", "steric"):
        findings.extend(steric.run_all(topology, positions, input_path))

    # Severity filter.
    threshold = Severity(args.severity)
    findings = [f for f in findings if sev_at_least(f.severity, threshold)]

    n_atoms = sum(1 for _ in topology.atoms())
    n_residues = sum(1 for _ in topology.residues())
    n_chains = sum(1 for _ in topology.chains())
    report_text = format_report(
        str(input_path), n_atoms, n_residues, n_chains, findings,
    )

    if args.output:
        Path(args.output).write_text(report_text + "\n")
    else:
        print(report_text)

    # Exit 1 if any ERROR remains (after severity filter). Matches
    # ruff's convention — usable in shell gates.
    if any(f.severity == Severity.ERROR for f in findings):
        sys.exit(1)
    sys.exit(0)
