"""Orchestrator for ``dvbfixer diagnose``.

Loads the input PDB once, prepares a shared ``PDBFixer`` for the
structural checks that need it, dispatches to the three check
families, filters findings by severity, formats the report, and
exits with code 1 if any ERROR finding remains.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

from dvbfixer.diagnose import chemistry, steric, structural
from dvbfixer.diagnose.cli import parse_args
from dvbfixer.diagnose.report import (
    Finding,
    Severity,
    findings_to_dict_list,
    format_report,
    sev_at_least,
)


def _detect_multi_model(pdb_path: Path) -> int:
    """Return the number of MODEL records in the file. 0 means single-model."""
    count = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("MODEL "):
                count += 1
    return count


def _extract_first_model(pdb_path: Path, work_dir: Path) -> Path:
    """Copy MODEL 1 (through its ENDMDL) to a temp file. Non-ATOM/HETATM
    headers before the first MODEL are preserved so PDBFixer still sees
    SEQRES / CONECT / HELIX etc.
    """
    out = work_dir / "diagnose-model1.pdb"
    in_model = False
    seen_first = False
    with open(pdb_path) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith("MODEL "):
                if not seen_first:
                    seen_first = True
                    in_model = True
                    continue
                # second MODEL — stop
                break
            if line.startswith("ENDMDL"):
                if in_model:
                    in_model = False
                    seen_first = True  # done with first model
                continue
            if seen_first and not in_model:
                # Emit trailing CONECT / END records that live after
                # the first ENDMDL.
                if line.startswith(("CONECT", "END", "MASTER")):
                    fout.write(line)
                continue
            fout.write(line)
    return out


def _load_topology_and_fixer(pdb_path: Path) -> tuple[Any, Any, Any]:
    """Load the input via OpenMM ``PDBFile`` (for topology + positions)
    and via ``PDBFixer`` (for missing-atom / missing-terminals checks).
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    pdb = PDBFile(str(pdb_path))
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    # findMissingAtoms populates missingTerminals as a side effect in
    # current PDBFixer, but call it explicitly for future-proofing.
    if hasattr(fixer, "findMissingTerminals"):
        fixer.findMissingTerminals()
    return pdb.topology, pdb.positions, fixer


def _multi_model_finding(n_models: int) -> Finding:
    return Finding(
        severity=Severity.WARNING,
        category="multi_model",
        chain="*",
        resid="*",
        resname="*",
        message=f"input has {n_models} MODEL records; diagnose analysed "
                f"MODEL 1 only. Split the ensemble with `dvbfixer split` "
                f"to inspect other frames.",
        fix_hint="dvbfixer split (if per-frame analysis is desired)",
    )


def _run_family(
    name: str,
    fn: Any,
    verbose: bool,
) -> list[Finding]:
    if not verbose:
        return list(fn())
    t0 = time.time()
    out = list(fn())
    print(f"[diagnose] {name}: {len(out)} findings in "
          f"{(time.time() - t0):.2f}s", file=sys.stderr)
    return out


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(2)

    # Multi-model handling: extract MODEL 1 so every downstream check
    # (raw-text and OpenMM) sees a single, coherent frame.
    n_models = _detect_multi_model(input_path)
    working_input = input_path
    tmp_pdb: Path | None = None
    prelim: list[Finding] = []
    if n_models > 1:
        import tempfile
        tmp_dir = Path(tempfile.mkdtemp(prefix="dvbfixer-diagnose-"))
        tmp_pdb = _extract_first_model(input_path, tmp_dir)
        working_input = tmp_pdb
        prelim.append(_multi_model_finding(n_models))

    try:
        topology, positions, fixer = _load_topology_and_fixer(working_input)
    except Exception as exc:
        print(f"Failed to load {input_path}: {exc}", file=sys.stderr)
        sys.exit(2)

    include_water: bool = args.include_water
    verbose: bool = args.verbose

    # Clash cutoffs: --clash-cutoff wins over --clash-mode; if neither
    # is set, fall through to the preset default.
    from dvbfixer.diagnose.steric import CLASH_MODE_PRESETS
    if args.clash_cutoff is not None:
        clash_warn_a, clash_error_a = args.clash_cutoff
    else:
        clash_warn_a, clash_error_a = CLASH_MODE_PRESETS[args.clash_mode]

    findings: list[Finding] = list(prelim)
    if args.only in ("all", "structural"):
        findings.extend(_run_family(
            "structural",
            lambda: structural.run_all(
                topology, positions, working_input, fixer,
                include_water=include_water,
            ),
            verbose,
        ))
    if args.only in ("all", "chemistry"):
        findings.extend(_run_family(
            "chemistry",
            lambda: chemistry.run_all(topology, positions),
            verbose,
        ))
    if args.only in ("all", "steric"):
        findings.extend(_run_family(
            "steric",
            lambda: steric.run_all(
                topology, positions, working_input,
                include_water=include_water,
                clash_warn_a=clash_warn_a,
                clash_error_a=clash_error_a,
            ),
            verbose,
        ))

    # Severity filter.
    threshold = Severity(args.severity)
    findings = [f for f in findings if sev_at_least(f.severity, threshold)]

    n_atoms = sum(1 for _ in topology.atoms())
    n_residues = sum(1 for _ in topology.residues())
    n_chains = sum(1 for _ in topology.chains())

    if args.output_format == "json":
        payload = {
            "input": str(input_path),
            "n_models": n_models if n_models > 1 else 1,
            "n_atoms": n_atoms,
            "n_residues": n_residues,
            "n_chains": n_chains,
            "findings": findings_to_dict_list(findings),
            "summary": {
                sev.value: sum(1 for f in findings if f.severity == sev)
                for sev in Severity
            },
        }
        report_text = json.dumps(payload, indent=2)
    else:
        report_text = format_report(
            str(input_path), n_atoms, n_residues, n_chains, findings,
        )

    if args.output:
        Path(args.output).write_text(report_text + "\n")
    else:
        print(report_text)

    # Clean up temp MODEL-1 file.
    if tmp_pdb is not None:
        try:
            tmp_pdb.unlink()
            tmp_pdb.parent.rmdir()
        except OSError:
            pass

    if any(f.severity == Severity.ERROR for f in findings):
        sys.exit(1)
    sys.exit(0)
