"""Steric clash detection for ``dvbfixer diagnose``.

Two engines:

1. **Pure Python** (default): ``scipy.spatial.cKDTree.query_pairs``
   for O(N log N) neighbor search + vdW radii lookup. No new deps
   (scipy is already a dvbfixer dep).

2. **MolProbity ``probe``** (auto-used if the binary is on ``PATH``):
   shells out to ``probe`` for a MolProbity-quality all-atom contact
   analysis. Falls back to the Python engine on any error.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from dvbfixer.diagnose.report import Finding, Severity

# van der Waals radii (Å). Standard values from the MolProbity dataset
# (contact_dots.py) — matches most other structure-quality tools.
_VDW_A: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "SE": 1.90,
    "MG": 1.73,
    "CA": 2.31,
    "NA": 2.27,
    "K": 2.75,
    "ZN": 1.39,
    "FE": 2.00,
}

# MolProbity-standard clash cutoffs (Å of vdW overlap).
_CLASH_ERROR_OVERLAP_A = 0.5
_CLASH_WARNING_OVERLAP_A = 0.2

# Neighbor-list cutoff (Å). Should be > 2 * max vdW radius to catch all
# possible overlaps. 3.0 Å covers everything up to S-S contact.
_NEIGHBOR_CUTOFF_A = 3.0


def _resid_str(res_id: str, icode: str = "") -> str:
    icode = (icode or "").strip()
    return f"{res_id}{icode}" if icode else str(res_id)


def _res_loc(res: Any) -> tuple[str, str]:
    icode = getattr(res, "insertionCode", "").strip()
    return res.chain.id, _resid_str(res.id, icode)


def _find_binary(name: str) -> str | None:
    """Locate ``name`` on PATH or in the current Python env's bin dir."""
    import os
    import shutil
    import sys
    found = shutil.which(name)
    if found:
        return found
    py_bin = os.path.dirname(sys.executable)
    candidate = os.path.join(py_bin, name)
    if os.access(candidate, os.X_OK):
        return candidate
    return None


def _build_bond_exclusion(topology: Any) -> set[frozenset[int]]:
    """Return {frozenset({i, j})} for every 1-2, 1-3, and 1-4 pair in
    the topology — pairs we must EXCLUDE from clash detection because
    they aren't sterically independent.

    1-4 exclusion matches AMBER's default fudgeLJ=0.5 treatment (which
    already discounts 1-4 interactions) and prevents false-positive
    clashes on tight but chemically-valid rotamers (e.g. SER CB-HB2
    vs OG-HG at ~2.1 Å).
    """
    neighbors: dict[int, set[int]] = {}
    for b1, b2 in topology.bonds():
        neighbors.setdefault(b1.index, set()).add(b2.index)
        neighbors.setdefault(b2.index, set()).add(b1.index)

    excluded: set[frozenset[int]] = set()
    # 1-2 (direct bonds)
    for i, js in neighbors.items():
        for j in js:
            excluded.add(frozenset({i, j}))
    # 1-3 (via one bond hop)
    for i, js in neighbors.items():
        for k in js:
            for j in neighbors.get(k, ()):
                if j == i:
                    continue
                excluded.add(frozenset({i, j}))
    # 1-4 (via two bond hops)
    for i, js in neighbors.items():
        for k1 in js:
            for k2 in neighbors.get(k1, ()):
                if k2 == i:
                    continue
                for j in neighbors.get(k2, ()):
                    if j == i or j == k1:
                        continue
                    excluded.add(frozenset({i, j}))
    return excluded


def clashes_python(topology: Any, positions: Any) -> list[Finding]:
    """Pure-Python clash detection via scipy cKDTree.

    For every non-bonded, non-1-3 pair within ``_NEIGHBOR_CUTOFF_A``,
    computes vdW overlap = (r_i + r_j) - d. Overlap ≥ 0.5 Å → ERROR;
    0.2 ≤ overlap < 0.5 → WARNING; below 0.2 Å → not reported (below
    MolProbity's noise floor).
    """
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except ImportError:
        return []
    from openmm.unit import nanometer

    atoms = list(topology.atoms())
    n = len(atoms)
    if n == 0:
        return []

    coords_nm = np.array([positions[a.index].value_in_unit(nanometer) for a in atoms])
    coords_a = coords_nm * 10.0

    tree = cKDTree(coords_a)
    # Query all pairs within the cutoff.
    pairs = tree.query_pairs(r=_NEIGHBOR_CUTOFF_A)
    if not pairs:
        return []

    excluded = _build_bond_exclusion(topology)

    findings: list[Finding] = []
    for i, j in pairs:
        ai = atoms[i]
        aj = atoms[j]
        if ai.element is None or aj.element is None:
            continue
        if frozenset({ai.index, aj.index}) in excluded:
            continue
        ri = _VDW_A.get(ai.element.symbol.upper(), 1.7)
        rj = _VDW_A.get(aj.element.symbol.upper(), 1.7)
        d = float(np.linalg.norm(coords_a[i] - coords_a[j]))
        overlap = (ri + rj) - d
        if overlap < _CLASH_WARNING_OVERLAP_A:
            continue
        severity = (Severity.ERROR
                    if overlap >= _CLASH_ERROR_OVERLAP_A
                    else Severity.WARNING)
        chain_i, resid_i = _res_loc(ai.residue)
        chain_j, resid_j = _res_loc(aj.residue)
        loc_j = f"{chain_j}/{aj.residue.name}{resid_j}:{aj.name}"
        findings.append(Finding(
            severity=severity,
            category="clash",
            chain=chain_i,
            resid=resid_i,
            resname=ai.residue.name,
            atom=ai.name,
            message=f"clashes with {loc_j} — overlap {overlap:.2f} Å "
                    f"(d={d:.2f} Å, vdW sum {ri + rj:.2f} Å)",
            fix_hint="dvbfixer minimize (relaxes non-bonded interactions)",
        ))
    return findings


def _run_probe(pdb_path: str | Path, binary: str) -> str | None:
    """Invoke ``probe`` and return its stdout, or None on error."""
    try:
        result = subprocess.run(
            [binary, "-unformated", "-condensed", "-Q", "-mc", "-het",
             "ALL", "ALL", str(pdb_path)],
            capture_output=True, text=True, timeout=60, check=False,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def clashes_probe(pdb_path: str | Path, binary: str) -> list[Finding] | None:
    """MolProbity ``probe`` clash detection.

    Parses probe's condensed output (one line per contact). Only "bad
    overlaps" (bumps and clashes) are surfaced — H-bonds and wide
    contacts are skipped since diagnose is about problems, not
    interactions.

    Returns ``None`` on any probe failure so :func:`run_all` can fall
    back to the Python engine. Returns an empty list when probe ran
    cleanly and found no clashes.
    """
    raw = _run_probe(pdb_path, binary)
    if raw is None:
        return None

    findings: list[Finding] = []
    for line in raw.splitlines():
        # Probe's condensed format is colon-separated; the first token
        # is a contact-type label. "bo" = bad overlap; "so" = small
        # overlap. We treat both as clashes.
        parts = line.strip().split(":")
        if len(parts) < 13:
            continue
        contact_type = parts[2].strip() if len(parts) > 2 else ""
        if contact_type not in ("bo", "so"):
            continue
        # Fields: :name:type:srcAtom:trgAtom:mingap:...
        # srcAtom / trgAtom fields include chain, resnum, resname, atomname.
        try:
            src_field = parts[3].strip()  # e.g. " A 123 SER  HA "
            trg_field = parts[4].strip()
            mingap_str = parts[5].strip()
            overlap = -float(mingap_str)  # probe reports gap (negative on overlap)
        except (IndexError, ValueError):
            continue
        if overlap < _CLASH_WARNING_OVERLAP_A:
            continue
        severity = (Severity.ERROR
                    if overlap >= _CLASH_ERROR_OVERLAP_A
                    else Severity.WARNING)

        def _parse(field: str) -> tuple[str, str, str, str]:
            # Fixed-width: cols 0=chain, 1-4=resid, 5-7=resname, 8-11=atomname
            pad = field.ljust(20)
            return pad[0].strip(), pad[1:5].strip(), pad[5:9].strip(), pad[10:15].strip()

        chain_i, resid_i, resname_i, atom_i = _parse(src_field)
        chain_j, resid_j, resname_j, atom_j = _parse(trg_field)
        loc_j = f"{chain_j}/{resname_j}{resid_j}:{atom_j}"
        findings.append(Finding(
            severity=severity,
            category="clash",
            chain=chain_i or "?",
            resid=resid_i or "?",
            resname=resname_i or "?",
            atom=atom_i or "?",
            message=f"clashes with {loc_j} — overlap {overlap:.2f} Å "
                    f"(via MolProbity probe)",
            fix_hint="dvbfixer minimize",
        ))
    return findings


def run_all(
    topology: Any,
    positions: Any,
    pdb_path: str | Path,
    prefer_probe: bool = True,
) -> list[Finding]:
    """Run clash detection. Prefers ``probe`` when the binary is
    available; falls back to the Python engine otherwise.
    """
    if prefer_probe:
        binary = _find_binary("probe")
        if binary is not None:
            findings = clashes_probe(pdb_path, binary)
            if findings is not None:  # None signals probe run failure
                return findings
    return clashes_python(topology, positions)
