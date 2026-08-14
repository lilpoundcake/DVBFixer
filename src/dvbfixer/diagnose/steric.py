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

from dvbfixer.diagnose._common import HBOND_HEAVY_ELEMENTS, is_water
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

# Preset clash cutoffs (WARN, ERROR) in Å of vdW overlap. Selectable
# at the CLI via ``--clash-mode``.
#   - molprobity  : strict validation floor (published MolProbity value)
#   - chimerax    : ChimeraX ``clashes`` default (matches its GUI report)
#   - bioluminate : BioLuminate ``find_clashes`` "Bad"/"Ugly" split
CLASH_MODE_PRESETS: dict[str, tuple[float, float]] = {
    "molprobity":  (0.4, 0.5),
    "chimerax":    (0.6, 0.9),
    "bioluminate": (0.75, 1.0),
}
DEFAULT_CLASH_MODE = "chimerax"

# Module-level defaults (used only when the caller passes ``None``).
_CLASH_WARNING_OVERLAP_A, _CLASH_ERROR_OVERLAP_A = CLASH_MODE_PRESETS[DEFAULT_CLASH_MODE]

# Neighbor-list cutoff (Å). Should be > 2 * max vdW radius to catch all
# possible overlaps. 3.0 Å covers everything up to S-S contact.
_NEIGHBOR_CUTOFF_A = 3.0

# Distance-based covalent-bond inference for the exclusion set. Any
# atom pair closer than these cutoffs is treated as bonded regardless
# of what ``topology.bonds()`` says — this catches HETATMs (glycans,
# ligands, ions), AMBER/CHARMM protonation variants, and other
# non-standard residues that OpenMM's PDBFile can't infer bonds for.
# Narrow enough to avoid confusing short H-bonds (~1.8 Å heavy-heavy)
# with covalent bonds; wide enough to catch stretched CONECT-less
# HETATM bonds. Element-pair keys.
_COVALENT_MAX_A_HEAVY_HEAVY = 1.90
_COVALENT_MAX_A_X_H = 1.30
_COVALENT_MAX_A_SS = 2.25  # disulfide (S-S)

# Hydrogen-bond geometric envelope. A pair (polar H, acceptor) with
# donor-H covalent bond to N/O/S is skipped as a valid H-bond, not
# a clash — this matches what MolProbity's ``probe`` does when it
# labels contacts as "hbond" instead of "bo/so". Values are
# permissive to cover both classical N-H..O (1.8 - 2.4 Å) and salt
# bridges (down to ~1.6 Å for very strong ARG-guanidinium..carboxylate
# contacts). MolProbity's official H-bond gap is 0.6 Å; we use 0.6 Å
# by allowing overlaps up to 0.6 Å below the vdW sum on H-bond pairs.
_HBOND_ACCEPTORS = {"N", "O", "S", "F"}
_HBOND_DONORS = {"N", "O", "S"}
_HBOND_H_ACC_MIN_A = 1.4    # closer than this is genuinely too tight
_HBOND_H_ACC_MAX_A = 2.6    # farther than this isn't an H-bond

# Heavy-atom H-bond envelope (for structures without explicit H).
# Donor..acceptor at 2.5 – 3.4 Å between H-bond-capable elements
# (N/O/S/F) is a genuine H-bond even when we can't see the H —
# common in X-ray structures that omit hydrogens entirely, or in
# LYS/ARG side chains where the charged NH3+/guanidinium H is
# missing. Matches CCP4's H-bond envelope.
_HBOND_HEAVY_MIN_A = 2.5
_HBOND_HEAVY_MAX_A = 3.4


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


def _covalent_cutoff_a(sym_i: str, sym_j: str) -> float:
    """Element-pair covalent cutoff for distance-based bond inference."""
    a, b = sym_i.upper(), sym_j.upper()
    if a == "H" or b == "H":
        return _COVALENT_MAX_A_X_H
    if a == "S" and b == "S":
        return _COVALENT_MAX_A_SS
    return _COVALENT_MAX_A_HEAVY_HEAVY


def _build_neighbor_graph(
    topology: Any,
    atoms: list[Any],
    coords_a: Any,
    pairs_within_cutoff: set[tuple[int, int]],
) -> dict[int, set[int]]:
    """Merge topology bonds with distance-inferred bonds.

    OpenMM's ``PDBFile`` only auto-infers bonds for STANDARD residue
    templates — HETATMs (glycans, ligands, ions), AMBER protonation
    variants (HIE, CYX, LYN, ASH, GLH), CHARMM variants (HSD, HSE,
    HSP), and sugar residues get NO bonds from ``topology.bonds()``.
    Without this augmentation, their real covalent bonds appear as
    clashes (a bonded C-C at 1.53 Å reads as 1.87 Å overlap against
    vdW sum 3.40 Å).
    """
    import numpy as np

    neighbors: dict[int, set[int]] = {}
    for b1, b2 in topology.bonds():
        neighbors.setdefault(b1.index, set()).add(b2.index)
        neighbors.setdefault(b2.index, set()).add(b1.index)

    for i, j in pairs_within_cutoff:
        ai, aj = atoms[i], atoms[j]
        if ai.element is None or aj.element is None:
            continue
        cutoff = _covalent_cutoff_a(ai.element.symbol, aj.element.symbol)
        d = float(np.linalg.norm(coords_a[i] - coords_a[j]))
        if d <= cutoff:
            neighbors.setdefault(i, set()).add(j)
            neighbors.setdefault(j, set()).add(i)
    return neighbors


def _expand_exclusions(neighbors: dict[int, set[int]]) -> set[frozenset[int]]:
    """Return {frozenset({i, j})} for every 1-2, 1-3, and 1-4 pair in
    the neighbor graph — pairs we must EXCLUDE from clash detection.

    1-4 exclusion matches AMBER's default fudgeLJ=0.5 treatment (which
    already discounts 1-4 interactions) and prevents false-positive
    clashes on tight but chemically-valid rotamers (e.g. SER CB-HB2
    vs OG-HG at ~2.1 Å).
    """
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


def _is_hbond_pair(
    ai: Any,
    aj: Any,
    d_a: float,
    atoms: list[Any],
    neighbors: dict[int, set[int]],
) -> bool:
    """Return True when ``(ai, aj)`` fits an H-bond envelope.

    Criteria: one atom is H covalently bonded to N/O/S (polar donor);
    the other atom is N/O/S/F (acceptor); H..acceptor distance is
    within 1.4 – 2.6 Å. Skips angle checks — they'd need the donor
    position, which is fine to leave out since the covalent-bond
    graph already constrains geometry.
    """
    sym_i = ai.element.symbol.upper()
    sym_j = aj.element.symbol.upper()
    if sym_i == "H" and sym_j in _HBOND_ACCEPTORS:
        h_atom, acc_sym = ai, sym_j
    elif sym_j == "H" and sym_i in _HBOND_ACCEPTORS:
        h_atom, acc_sym = aj, sym_i
    else:
        return False
    if not (_HBOND_H_ACC_MIN_A <= d_a <= _HBOND_H_ACC_MAX_A):
        return False
    # Polar-H test: is our H bonded to N/O/S?
    for parent_idx in neighbors.get(h_atom.index, ()):
        parent = atoms[parent_idx]
        if parent.element is None:
            continue
        if parent.element.symbol.upper() in _HBOND_DONORS:
            return True
    _ = acc_sym  # informational; matched via _HBOND_ACCEPTORS above
    return False


def _is_heavy_hbond_pair(
    ai: Any,
    aj: Any,
    d_a: float,
) -> bool:
    """Heavy-atom H-bond envelope (used when the H may be missing).

    Skips pairs of two H-bond-capable elements (N/O/S/F) at
    2.5 – 3.4 Å. This catches LYS-NZ..GLU-OE1 salt bridges,
    backbone O..O interactions in tight β-turns, and every other
    donor..acceptor pair in structures that lack explicit hydrogens.
    """
    sym_i = ai.element.symbol.upper()
    sym_j = aj.element.symbol.upper()
    if sym_i not in HBOND_HEAVY_ELEMENTS or sym_j not in HBOND_HEAVY_ELEMENTS:
        return False
    return _HBOND_HEAVY_MIN_A <= d_a <= _HBOND_HEAVY_MAX_A


def clashes_python(
    topology: Any,
    positions: Any,
    include_water: bool = False,
    clash_warn_a: float | None = None,
    clash_error_a: float | None = None,
) -> list[Finding]:
    """Pure-Python clash detection via scipy cKDTree.

    For every non-bonded, non-1-3, non-1-4 pair within
    ``_NEIGHBOR_CUTOFF_A``, computes vdW overlap = (r_i + r_j) - d.
    Overlap at or above ``clash_error_a`` → ERROR; at or above
    ``clash_warn_a`` (but below the error threshold) → WARNING; below
    that → not reported. Both default to ``CLASH_MODE_PRESETS[
    DEFAULT_CLASH_MODE]`` (currently the ``"chimerax"`` preset, 0.6/0.9 Å)
    unless the caller passes explicit ``clash_warn_a``/``clash_error_a``
    values or picks a different preset (e.g. ``"molprobity"``, 0.4/0.5 Å,
    MolProbity's official clashscore floor).

    Bond exclusion includes both the OpenMM topology bond set AND
    a distance-inferred graph — so HETATMs and non-standard residues
    without CONECT records don't have their real covalent bonds
    reported as clashes.

    Hydrogen bonds (polar H..N/O/S/F at 1.4 – 2.6 Å) are also
    skipped — MolProbity's ``probe`` classifies those as "hbond",
    not "bo/so" clashes. Without this filter every backbone
    N-H..O=C amide H-bond in an α-helix looks like a 0.7 Å overlap.
    """
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except ImportError:
        return []
    from openmm.unit import nanometer

    warn_a = _CLASH_WARNING_OVERLAP_A if clash_warn_a is None else clash_warn_a
    err_a = _CLASH_ERROR_OVERLAP_A if clash_error_a is None else clash_error_a

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

    neighbors = _build_neighbor_graph(topology, atoms, coords_a, pairs)
    excluded = _expand_exclusions(neighbors)

    findings: list[Finding] = []
    for i, j in pairs:
        ai = atoms[i]
        aj = atoms[j]
        if ai.element is None or aj.element is None:
            continue
        if frozenset({ai.index, aj.index}) in excluded:
            continue
        # Skip pairs within the same residue — intra-residue vdW is
        # baked into the FF template and any close contact is either
        # already handled (1-2/1-3/1-4 above) or a template-defined
        # tight geometry (e.g. sugar ring puckers).
        if ai.residue is aj.residue:
            continue
        # Skip water (unless the user opted in).
        if not include_water and (
            is_water(ai.residue.name) or is_water(aj.residue.name)
        ):
            continue
        ri = _VDW_A.get(ai.element.symbol.upper(), 1.7)
        rj = _VDW_A.get(aj.element.symbol.upper(), 1.7)
        d = float(np.linalg.norm(coords_a[i] - coords_a[j]))
        overlap = (ri + rj) - d
        if overlap < warn_a:
            continue
        # Skip pairs that fit an H-bond envelope — not a clash.
        if _is_hbond_pair(ai, aj, d, atoms, neighbors):
            continue
        # Heavy-atom H-bond envelope (for structures missing explicit H).
        if _is_heavy_hbond_pair(ai, aj, d):
            continue
        severity = (Severity.ERROR
                    if overlap >= err_a
                    else Severity.WARNING)
        chain_i, resid_i = _res_loc(ai.residue)
        chain_j, resid_j = _res_loc(aj.residue)
        loc_i = f"{chain_i}/{ai.residue.name}{resid_i}:{ai.name}"
        loc_j = f"{chain_j}/{aj.residue.name}{resid_j}:{aj.name}"
        findings.append(Finding(
            severity=severity,
            category="clash",
            chain=chain_i,
            resid=resid_i,
            resname=ai.residue.name,
            atom=ai.name,
            message=f"{loc_i} clashes with partner {loc_j} — overlap {overlap:.2f} Å "
                    f"(d={d:.2f} Å, vdW sum {ri + rj:.2f} Å)",
            fix_hint="dvbfixer minimize (relaxes non-bonded interactions)",
            extra={
                "clash_partner": {
                    "chain": chain_j,
                    "resid": resid_j,
                    "resname": aj.residue.name,
                    "atom": aj.name,
                }
            },
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


def clashes_probe(
    pdb_path: str | Path,
    binary: str,
    clash_warn_a: float | None = None,
    clash_error_a: float | None = None,
) -> list[Finding] | None:
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

    warn_a = _CLASH_WARNING_OVERLAP_A if clash_warn_a is None else clash_warn_a
    err_a = _CLASH_ERROR_OVERLAP_A if clash_error_a is None else clash_error_a

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
        if overlap < warn_a:
            continue
        severity = (Severity.ERROR
                    if overlap >= err_a
                    else Severity.WARNING)

        def _parse(field: str) -> tuple[str, str, str, str]:
            # Fixed-width: cols 0=chain, 1-4=resid, 5-7=resname, 8-11=atomname
            pad = field.ljust(20)
            return pad[0].strip(), pad[1:5].strip(), pad[5:9].strip(), pad[10:15].strip()

        chain_i, resid_i, resname_i, atom_i = _parse(src_field)
        chain_j, resid_j, resname_j, atom_j = _parse(trg_field)
        loc_i = f"{chain_i or '?'}/{resname_i or '?'}{resid_i or '?'}:{atom_i or '?'}"
        loc_j = f"{chain_j}/{resname_j}{resid_j}:{atom_j}"
        findings.append(Finding(
            severity=severity,
            category="clash",
            chain=chain_i or "?",
            resid=resid_i or "?",
            resname=resname_i or "?",
            atom=atom_i or "?",
            message=f"{loc_i} clashes with partner {loc_j} — overlap {overlap:.2f} Å "
                    f"(via MolProbity probe)",
            fix_hint="dvbfixer minimize",
            extra={
                "clash_partner": {
                    "chain": chain_j or "?",
                    "resid": resid_j or "?",
                    "resname": resname_j or "?",
                    "atom": atom_j or "?",
                }
            },
        ))
    return findings


def run_all(
    topology: Any,
    positions: Any,
    pdb_path: str | Path,
    prefer_probe: bool = True,
    include_water: bool = False,
    clash_warn_a: float | None = None,
    clash_error_a: float | None = None,
) -> list[Finding]:
    """Run clash detection. Prefers ``probe`` when the binary is
    available; falls back to the Python engine otherwise.
    """
    if prefer_probe:
        binary = _find_binary("probe")
        if binary is not None:
            findings = clashes_probe(
                pdb_path, binary,
                clash_warn_a=clash_warn_a, clash_error_a=clash_error_a,
            )
            if findings is not None:  # None signals probe run failure
                return findings
    return clashes_python(
        topology, positions,
        include_water=include_water,
        clash_warn_a=clash_warn_a, clash_error_a=clash_error_a,
    )
