"""Structural-integrity checks for ``dvbfixer diagnose``.

Reuses existing detection helpers:
- :func:`dvbfixer.ffutils.geometry.detect_coincident_atoms` — H/heavy
  pairs coincident in the same residue.
- :func:`dvbfixer.ffutils.geometry.detect_misplaced_hydrogens` — H with
  wrong bond geometry.
- PDBFixer's ``findMissingResidues`` / ``findMissingAtoms`` /
  ``findMissingTerminals`` — missing structural pieces.
- :func:`dvbfixer.split_chains.find_chain_breaks` — chain-break
  heuristics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dvbfixer.diagnose._common import is_water
from dvbfixer.diagnose.report import Finding, Severity


def _resid_str(res_id: str, icode: str = "") -> str:
    icode = (icode or "").strip()
    return f"{res_id}{icode}" if icode else str(res_id)


def _res_loc(res: Any) -> tuple[str, str]:
    """Return (chain_id, resid+icode) for a residue."""
    icode = getattr(res, "insertionCode", "").strip()
    return res.chain.id, _resid_str(res.id, icode)


def check_missing_atoms(fixer: Any) -> list[Finding]:
    """PDBFixer's ``findMissingAtoms`` — missing heavy atoms in
    existing residues. Terminal atoms (OXT) are reported by
    :func:`check_missing_terminals` instead.
    """
    findings: list[Finding] = []
    for res, atoms in fixer.missingAtoms.items():
        chain, resid = _res_loc(res)
        names = ", ".join(a.name for a in atoms)
        findings.append(Finding(
            severity=Severity.ERROR,
            category="missing_atoms",
            chain=chain,
            resid=resid,
            resname=res.name,
            message=f"missing heavy atoms: {names}",
            fix_hint="dvbfixer prepare",
        ))
    return findings


def check_missing_residues(fixer: Any) -> list[Finding]:
    """SEQRES-vs-ATOM gaps that PDBFixer would fill via
    ``addMissingResidues``.
    """
    findings: list[Finding] = []
    chains_list = list(fixer.topology.chains())
    for (chain_idx, res_idx), resnames in fixer.missingResidues.items():
        chain_id = chains_list[chain_idx].id if chain_idx < len(chains_list) else "?"
        run = " ".join(resnames)
        findings.append(Finding(
            severity=Severity.ERROR,
            category="missing_residues",
            chain=chain_id,
            resid=str(res_idx),
            resname=resnames[0] if resnames else "?",
            message=f"gap between existing residues: {len(resnames)} residues "
                    f"({run}) not present in ATOM records",
            fix_hint="dvbfixer model (uses Modeller LoopModel)",
        ))
    return findings


def check_missing_terminals(fixer: Any) -> list[Finding]:
    """Terminal atoms (OXT on C-terminus, H2/H3 on N-terminus) missing."""
    findings: list[Finding] = []
    for res, atoms in fixer.missingTerminals.items():
        chain, resid = _res_loc(res)
        names = ", ".join(a if isinstance(a, str) else a.name for a in atoms)
        findings.append(Finding(
            severity=Severity.ERROR,
            category="missing_terminals",
            chain=chain,
            resid=resid,
            resname=res.name,
            message=f"missing terminal atom(s): {names}",
            fix_hint="dvbfixer prepare",
        ))
    return findings


def check_coincident_atoms(topology: Any, positions: Any) -> list[Finding]:
    """H atoms placed on top of a heavy atom in the same residue.

    Reuses :func:`dvbfixer.ffutils.geometry.detect_coincident_atoms`.
    Reported case: SER HG dumped at OXT position by an upstream tool
    (test/broken_SER/SER.pdb).
    """
    from dvbfixer.ffutils.geometry import detect_coincident_atoms
    findings: list[Finding] = []
    for coinc in detect_coincident_atoms(topology, positions):
        chain, resid = _res_loc(coinc.residue)
        findings.append(Finding(
            severity=Severity.ERROR,
            category="coincident_atoms",
            chain=chain,
            resid=resid,
            resname=coinc.residue.name,
            atom=coinc.hydrogen.name,
            message=f"coincident with {coinc.heavy_atom.name} "
                    f"({coinc.distance_nm * 10:.3f} Å apart)",
            fix_hint="dvbfixer prepare (0.4.1 pre-strip re-places both)",
        ))
    return findings


def check_misplaced_hydrogens(topology: Any, positions: Any) -> list[Finding]:
    """H atoms placed too far from their bonded heavy-atom parent.

    Reuses :func:`dvbfixer.ffutils.geometry.detect_misplaced_hydrogens`.
    """
    from dvbfixer.ffutils.geometry import detect_misplaced_hydrogens
    findings: list[Finding] = []
    for m in detect_misplaced_hydrogens(topology, positions):
        # Coincident cases are already reported by check_coincident_atoms —
        # skip the duplicate.
        if m.reason == "coincident":
            continue
        chain, resid = _res_loc(m.h_atom.residue)
        actual_a = m.parent_distance_nm * 10
        expected_a = m.expected_bond_nm * 10
        findings.append(Finding(
            severity=Severity.ERROR,
            category="misplaced_hydrogen",
            chain=chain,
            resid=resid,
            resname=m.h_atom.residue.name,
            atom=m.h_atom.name,
            message=f"{actual_a:.2f} Å from parent {m.parent.name} "
                    f"(expected ~{expected_a:.2f} Å)",
            fix_hint="dvbfixer prepare (0.4.2 post-check re-places)",
        ))
    return findings


def check_altloc_conflicts(pdb_path: str | Path) -> list[Finding]:
    """Report any residue whose input PDB carries alt-loc labels.

    Parses the raw PDB text (altLoc is column 17, one-based; index 16
    zero-based). Any non-blank altLoc on any atom flags the residue.
    """
    findings: list[Finding] = []
    seen: set[tuple[str, str, str]] = set()  # (chain, resid, altlocs-joined)
    per_res: dict[tuple[str, str, str], set[str]] = {}
    per_res_names: dict[tuple[str, str, str], str] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 17:
            continue
        altloc = line[16]
        if altloc == " ":
            continue
        chain = line[21]
        resid = line[22:26].strip()
        icode = line[26].strip() if len(line) > 26 else ""
        resname = line[17:20].strip()
        key = (chain, resid, icode)
        per_res.setdefault(key, set()).add(altloc)
        per_res_names[key] = resname

    for (chain, resid, icode), altlocs in per_res.items():
        joined = ", ".join(sorted(altlocs))
        findings.append(Finding(
            severity=Severity.WARNING,
            category="altloc_conflict",
            chain=chain,
            resid=_resid_str(resid, icode),
            resname=per_res_names[(chain, resid, icode)],
            message=f"alternative locations: {joined}",
            fix_hint="manual: pick one altLoc (e.g. via pymol / PDB editor)",
        ))
    return findings


def check_chain_breaks(
    pdb_path: str | Path,
    include_water: bool = False,
) -> list[Finding]:
    """Report inter-residue distance jumps that split a chain.

    Uses :func:`dvbfixer.split_chains.find_chain_breaks` heuristics
    (C→N > 2.5 Å, resSeq backwards, nearest-atom > 15 Å). Water
    residues generate massive noise (every ordered water is a chain
    break vs its neighbour) and are excluded unless ``include_water``.
    """
    findings: list[Finding] = []
    try:
        from dvbfixer.split_chains import (
            DEFAULT_DISTANCE_CUTOFF,
            DEFAULT_GAP_CUTOFF,
            find_chain_breaks,
        )
    except Exception:
        return findings

    def _keep(line: str) -> bool:
        if not line.startswith(("ATOM  ", "HETATM")):
            return False
        if not include_water and is_water(line[17:20]):
            return False
        return True

    # Diagnose only internal breaks. Existing PDB chain-ID transitions are
    # already explicit boundaries and must never be compared as consecutive
    # residues (e.g. A/GLY86 followed by B/MET1). Keep file order within each
    # chain because the shared split helper intentionally ignores chain IDs
    # for its separate empirical-splitting use case.
    lines_by_chain: dict[str, list[str]] = {}
    with open(pdb_path) as f:
        for line in f:
            if _keep(line):
                lines_by_chain.setdefault(line[21], []).append(line)

    for chain_lines in lines_by_chain.values():
        try:
            breaks = find_chain_breaks(
                chain_lines,
                distance_cutoff=DEFAULT_DISTANCE_CUTOFF,
                gap_cutoff=DEFAULT_GAP_CUTOFF,
                use_distance=True,
            )
        except Exception:
            continue

        # `find_chain_breaks` returns atom-index break points; index zero is
        # the normal start of this chain, while later indices are findings.
        for idx in breaks:
            if idx == 0 or idx >= len(chain_lines):
                continue
            prev = chain_lines[idx - 1]
            curr = chain_lines[idx]
            chain = curr[21]
            resid = curr[22:26].strip()
            icode = curr[26].strip() if len(curr) > 26 else ""
            resname = curr[17:20].strip()
            prev_chain = prev[21]
            prev_resid = prev[22:26].strip()
            prev_resname = prev[17:20].strip()
            findings.append(Finding(
                severity=Severity.WARNING,
                category="chain_break",
                chain=chain,
                resid=_resid_str(resid, icode),
                resname=resname,
                message=f"chain break after {prev_chain}/{prev_resname}{prev_resid} "
                        f"(prev-C to curr-N distance jump or resSeq reset)",
                fix_hint="dvbfixer model (if the gap is real)",
            ))
    return findings


def check_duplicate_chain_coordinates(pdb_path: str | Path) -> list[Finding]:
    """Warn about coordinate-identical protein chains in one PDB model."""
    from dvbfixer.pdbutils.duplicates import duplicate_protein_chain_coordinates

    return [
        Finding(
            severity=Severity.WARNING,
            category="duplicate_chain_coordinates",
            chain=right,
            resid="*",
            resname="*",
            message=(f"protein chains {left} and {right} have identical coordinates "
                     f"for all {atom_count} protein atoms; frames may have been "
                     "concatenated without MODEL/ENDMDL separators"),
            fix_hint="restore MODEL/ENDMDL records or remove the duplicate chain before minimize",
        )
        for left, right, atom_count in duplicate_protein_chain_coordinates(pdb_path)
    ]


def check_insertion_codes(pdb_path: str | Path) -> list[Finding]:
    """Report residues carrying insertion codes.

    Common on antibody structures (Kabat / Chothia CDR numbering) —
    INFO only, not an error.
    """
    findings: list[Finding] = []
    per_chain: dict[str, list[str]] = {}
    per_chain_first: dict[str, tuple[str, str]] = {}  # (resname, resid+icode)

    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 27:
            continue
        icode = line[26].strip()
        if not icode:
            continue
        chain = line[21]
        resid = line[22:26].strip()
        resname = line[17:20].strip()
        key = _resid_str(resid, icode)
        per_chain.setdefault(chain, []).append(key)
        if chain not in per_chain_first:
            per_chain_first[chain] = (resname, key)

    for chain, keys in per_chain.items():
        n = len(set(keys))
        resname, first = per_chain_first[chain]
        findings.append(Finding(
            severity=Severity.INFO,
            category="insertion_codes",
            chain=chain,
            resid=first,
            resname=resname,
            message=f"{n} residue(s) with insertion codes on chain {chain} "
                    f"(often antibody CDR numbering — Kabat/Chothia)",
            fix_hint=None,
        ))
    return findings


def check_seqres_vs_atom(pdb_path: str | Path) -> list[Finding]:
    """Report residues present in SEQRES but absent from ATOM records.

    PDBFixer's ``findMissingResidues`` catches INTERNAL gaps only —
    terminal truncations (a Fab His-tag in SEQRES but not resolved in
    the ATOM section) are silently ignored. This check reads SEQRES
    directly, computes each chain's SEQRES length, compares to the
    number of unique ATOM residues on that chain, and flags the
    difference.
    """
    findings: list[Finding] = []
    seqres_counts: dict[str, int] = {}
    atom_resids: dict[str, set[str]] = {}

    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("SEQRES"):
            if len(line) < 20:
                continue
            chain = line[11]
            # Column 13-17 = total residues on this chain (per PDB spec).
            try:
                total = int(line[13:17].strip())
            except ValueError:
                continue
            seqres_counts[chain] = total
            continue
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
            chain = line[21]
            resid = line[22:27]  # includes iCode
            atom_resids.setdefault(chain, set()).add(resid)

    for chain, total in seqres_counts.items():
        observed = len(atom_resids.get(chain, ()))
        if observed >= total:
            continue
        missing = total - observed
        findings.append(Finding(
            severity=Severity.WARNING,
            category="seqres_gap",
            chain=chain,
            resid="*",
            resname="*",
            message=f"SEQRES declares {total} residues on chain {chain}; "
                    f"ATOM records cover {observed}. {missing} residue(s) "
                    f"missing — likely truncated terminal(s) (e.g. His-tag, "
                    f"disordered N/C-terminus).",
            fix_hint="dvbfixer model (LoopModel can rebuild terminal gaps)",
        ))
    return findings


def run_all(
    topology: Any,
    positions: Any,
    pdb_path: str | Path,
    fixer: Any | None = None,
    include_water: bool = False,
) -> list[Finding]:
    """Execute every structural check and return the concatenated findings.

    ``fixer`` is an already-loaded ``PDBFixer`` with
    ``findMissingResidues`` + ``findMissingAtoms`` +
    ``findMissingTerminals`` already called. ``pipeline.main`` builds
    it once and hands it to us to avoid re-parsing the PDB.
    """
    findings: list[Finding] = []
    if fixer is not None:
        findings.extend(check_missing_residues(fixer))
        findings.extend(check_missing_atoms(fixer))
        findings.extend(check_missing_terminals(fixer))
    findings.extend(check_coincident_atoms(topology, positions))
    findings.extend(check_misplaced_hydrogens(topology, positions))
    findings.extend(check_altloc_conflicts(pdb_path))
    findings.extend(check_chain_breaks(pdb_path, include_water=include_water))
    findings.extend(check_duplicate_chain_coordinates(pdb_path))
    findings.extend(check_insertion_codes(pdb_path))
    findings.extend(check_seqres_vs_atom(pdb_path))
    return findings
