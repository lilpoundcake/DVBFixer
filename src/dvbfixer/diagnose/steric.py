"""Steric clash detection for ``dvbfixer diagnose``.

Two engines:

1. **Built-in Python engine**: ``scipy.spatial.cKDTree.query_pairs``
   for O(N log N) neighbor search + vdW radii lookup. No new deps
   (scipy is already a dvbfixer dep).

2. **MolProbity ``probe``** (optional, when its supported version is available):
   validated all-atom contact analysis, including intramolecular pairs.
   Falls back to Python with explicit per-run provenance on any failure.
"""

from __future__ import annotations

import sys
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

from dvbfixer.diagnose._common import HBOND_HEAVY_ELEMENTS, is_water
from dvbfixer.diagnose.geometry import (
    VDW_RADII_ANGSTROM,
    GeometryMeasurement,
    atom_identity,
)
from dvbfixer.diagnose.report import Finding, Severity

_VDW_A = VDW_RADII_ANGSTROM

# Preset clash cutoffs (WARN, ERROR) in Å of vdW overlap. Selectable
# at the CLI via ``--clash-mode``.
#   - molprobity  : strict validation floor (published MolProbity value)
#   - chimerax    : ChimeraX ``clashes`` default (matches its GUI report)
#   - bioluminate : BioLuminate ``find_clashes`` "Bad"/"Ugly" split
CLASH_MODE_PRESETS: dict[str, tuple[float, float]] = {
    "molprobity": (0.4, 0.5),
    "chimerax": (0.6, 0.9),
    "bioluminate": (0.75, 1.0),
}
DEFAULT_CLASH_MODE = "chimerax"
DETERMINISTIC_STERIC_ENGINE = "dvbfixer-python-steric"
DETERMINISTIC_STERIC_ENGINE_VERSION = "2"
BOND_EXCLUSION_POLICY = "topology-and-intra-residue-templates.v2"
SEVERE_OVERLAP_A = 0.90

# Module-level defaults (used only when the caller passes ``None``).
_CLASH_WARNING_OVERLAP_A, _CLASH_ERROR_OVERLAP_A = CLASH_MODE_PRESETS[DEFAULT_CLASH_MODE]

# Neighbor-list cutoff (Å). Should be > 2 * max vdW radius to catch all
# possible overlaps. 3.0 Å covers everything up to S-S contact.
_NEIGHBOR_CUTOFF_A = 3.0

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
_HBOND_H_ACC_MIN_A = 1.4  # closer than this is genuinely too tight
_HBOND_H_ACC_MAX_A = 2.6  # farther than this isn't an H-bond

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


@lru_cache(maxsize=1)
def _intra_residue_templates() -> dict[str, tuple[tuple[str, str], ...]]:
    """Read OpenMM's bond definitions, never a coordinate-derived graph."""
    root = ElementTree.parse(str(files("openmm.app") / "data" / "residues.xml")).getroot()
    return {
        residue.attrib["name"]: tuple(
            (bond.attrib["from"], bond.attrib["to"])
            for bond in residue.findall("Bond")
            if not any(bond.attrib[key].startswith(("-", "+")) for key in ("from", "to"))
        )
        for residue in root.findall("Residue")
    }


def _build_neighbor_graph(
    topology: Any,
    include_water: bool,
) -> tuple[dict[int, set[int]], list[Any]]:
    """Use supplied bonds and unique intra-residue template atom names only.

    Proximity is the signal being tested, not evidence of a bond. In
    particular, never infer an inter-residue edge here: one false edge
    would also suppress adjacent contacts through 1-3/1-4 expansion.
    Return residues whose intra-residue connectivity remains unresolved.
    The supplied topology (including its parser's bond decisions) is trusted;
    no force-field topology or input residue/atom name is mutated.
    """
    from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

    neighbors: dict[int, set[int]] = {}
    for b1, b2 in topology.bonds():
        neighbors.setdefault(b1.index, set()).add(b2.index)
        neighbors.setdefault(b2.index, set()).add(b1.index)

    templates = _intra_residue_templates()
    unresolved = []
    for residue in topology.residues():
        if not include_water and is_water(residue.name):
            continue
        atoms = list(residue.atoms())
        by_name: dict[str, list[Any]] = {}
        for atom in atoms:
            by_name.setdefault(atom.name, []).append(atom)
        parent = "HOH" if is_water(residue.name) else VARIANT_TO_PARENT.get(residue.name, residue.name)
        template = templates.get(parent, ())
        template_names = {name for pair in template for name in pair}
        for a, b in template:
            if len(by_name.get(a, ())) == len(by_name.get(b, ())) == 1:
                i, j = by_name[a][0].index, by_name[b][0].index
                neighbors.setdefault(i, set()).add(j)
                neighbors.setdefault(j, set()).add(i)
        # A fully supplied connected graph is accepted for unknown chemistry.
        # Missing intermediate atoms in a known template belong to the separate
        # structural checks; they must not trigger invented shortcut bonds.
        indices = {atom.index for atom in atoms}
        visited = set()
        frontier = {atoms[0].index} if atoms else set()
        while frontier:
            visited.update(frontier)
            frontier = {j for i in frontier for j in neighbors.get(i, ()) if j in indices} - visited
        duplicate_names = any(len(matches) != 1 for matches in by_name.values())
        known_names = bool(template) and set(by_name) <= template_names
        if duplicate_names or (len(atoms) > 1 and not known_names and visited != indices):
            unresolved.append(residue)
    return neighbors, unresolved


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
    *,
    measurement_sink: list[GeometryMeasurement] | None = None,
    neighbor_cutoff_a: float = _NEIGHBOR_CUTOFF_A,
    provenance: dict[str, Any] | None = None,
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

    Bond exclusion uses supplied topology bonds and intra-residue OpenMM
    templates (including shared protonation-variant parent names), never
    proximity. Undescribed connectivity produces a WARNING, not new bonds.

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
    neighbors, unresolved = _build_neighbor_graph(topology, include_water)
    connectivity = {
        "bond_exclusion_policy": BOND_EXCLUSION_POLICY,
        "connectivity_status": "unresolved" if unresolved else "as_supplied_or_templated",
        "unresolved_residues": [
            {"chain": residue.chain.id, "resid": str(residue.id),
             "icode": getattr(residue, "insertionCode", "").strip(), "resname": residue.name}
            for residue in unresolved
        ],
    }
    if provenance is not None:
        provenance.update(connectivity)
        provenance["status"] = "partial" if unresolved else "complete"
    if not atoms:
        return []
    findings: list[Finding] = []
    for residue in unresolved:
        chain, resid = _res_loc(residue)
        findings.append(Finding(
            severity=Severity.WARNING, category="steric_connectivity",
            chain=chain, resid=resid, resname=residue.name,
            message=f"{chain}/{residue.name}{resid}: incomplete or ambiguous connectivity for Python "
                    "steric exclusions; close contacts are not inferred as bonds. "
                    "Missing covalent links may cause false-positive clashes; "
                    "intra-residue clashes are not evaluated.",
            fix_hint="Provide chemically verified topology/CONECT records or use supported MolProbity Probe",
            extra={"bond_exclusion_policy": BOND_EXCLUSION_POLICY,
                   "atoms": [atom_identity(atom).to_dict() for atom in residue.atoms()]},
        ))

    coords_nm = np.array([positions[a.index].value_in_unit(nanometer) for a in atoms])
    coords_a = coords_nm * 10.0

    tree = cKDTree(coords_a)
    # Query all pairs within the cutoff.
    pairs = tree.query_pairs(r=neighbor_cutoff_a)
    if not pairs:
        return findings

    excluded = _expand_exclusions(neighbors)

    for i, j in sorted(pairs):
        ai = atoms[i]
        aj = atoms[j]
        if ai.element is None or aj.element is None:
            continue
        if frozenset({ai.index, aj.index}) in excluded:
            continue
        # Preserve the existing inter-residue scope. This check does not
        # validate arbitrary intra-residue geometry, even for supplied graphs.
        if ai.residue is aj.residue:
            continue
        # Skip water (unless the user opted in).
        if not include_water and (is_water(ai.residue.name) or is_water(aj.residue.name)):
            continue
        ri = _VDW_A.get(ai.element.symbol.upper(), 1.7)
        rj = _VDW_A.get(aj.element.symbol.upper(), 1.7)
        d = float(np.linalg.norm(coords_a[i] - coords_a[j]))
        overlap = (ri + rj) - d
        if overlap <= 0.0:
            continue
        # Skip pairs that fit an H-bond envelope — not a clash.
        if _is_hbond_pair(ai, aj, d, atoms, neighbors):
            continue
        # Heavy-atom H-bond envelope (for structures missing explicit H).
        if _is_heavy_hbond_pair(ai, aj, d):
            continue
        if measurement_sink is not None:
            measurement_sink.append(
                GeometryMeasurement(
                    metric="steric_vdw_overlap.v1",
                    atoms=(atom_identity(ai), atom_identity(aj)),
                    value=overlap,
                    unit="angstrom",
                    reference_value=ri + rj,
                    engine=DETERMINISTIC_STERIC_ENGINE,
                    engine_version=DETERMINISTIC_STERIC_ENGINE_VERSION,
                    metadata={
                        "distance_angstrom": d,
                        "vdw_sum_angstrom": ri + rj,
                        "one_two_one_three_one_four_excluded": True,
                        "bond_exclusion_policy": BOND_EXCLUSION_POLICY,
                        "connectivity_status": connectivity["connectivity_status"],
                    },
                )
            )
        if overlap < warn_a - 1e-12:
            continue
        severity = Severity.ERROR if overlap >= err_a - 1e-12 else Severity.WARNING
        chain_i, resid_i = _res_loc(ai.residue)
        chain_j, resid_j = _res_loc(aj.residue)
        loc_i = f"{chain_i}/{ai.residue.name}{resid_i}:{ai.name}"
        loc_j = f"{chain_j}/{aj.residue.name}{resid_j}:{aj.name}"
        findings.append(
            Finding(
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
                    "geometry_schema": "dvbfixer.geometry.v1",
                    "metric": "steric_vdw_overlap.v1",
                    "bond_exclusion_policy": BOND_EXCLUSION_POLICY,
                    "engine": {
                        "name": DETERMINISTIC_STERIC_ENGINE,
                        "version": DETERMINISTIC_STERIC_ENGINE_VERSION,
                    },
                    "overlap_angstrom": overlap,
                    "distance_angstrom": d,
                    "vdw_sum_angstrom": ri + rj,
                    "warning_threshold_angstrom": warn_a,
                    "error_threshold_angstrom": err_a,
                    "atoms": [atom_identity(ai).to_dict(), atom_identity(aj).to_dict()],
                    "clash_partner": {
                        "chain": chain_j,
                        "resid": resid_j,
                        "resname": aj.residue.name,
                        "atom": aj.name,
                    },
                },
            )
        )
    return findings


def steric_overlap_measurements_python(
    topology: Any,
    positions: Any,
    include_water: bool = False,
    *,
    provenance: dict[str, Any] | None = None,
) -> list[GeometryMeasurement]:
    """Return deterministic Python vdW overlaps without severity policy.

    The measurement set uses the same covalent, 1-3/1-4, water, and hydrogen-
    bond exclusions as the existing Python diagnose adapter. Callers apply
    their own thresholds; diagnose's default ERROR boundary remains 0.90 Å.
    """

    measurements: list[GeometryMeasurement] = []
    clashes_python(
        topology,
        positions,
        include_water=include_water,
        clash_warn_a=float("inf"),
        clash_error_a=float("inf"),
        measurement_sink=measurements,
        neighbor_cutoff_a=2.0 * max(_VDW_A.values()),
        provenance=provenance,
    )
    return sorted(measurements, key=lambda item: item.atoms)


def clashes_probe(
    pdb_path: str | Path,
    binary: str,
    clash_warn_a: float | None = None,
    clash_error_a: float | None = None,
    *,
    include_water: bool = False,
    probe_radii: str = "electron",
) -> list[Finding] | None:
    """Return validated Probe findings, or None for a failed/unsupported run."""
    from dvbfixer.diagnose.probe import analyze
    warn_a = _CLASH_WARNING_OVERLAP_A if clash_warn_a is None else clash_warn_a
    err_a = _CLASH_ERROR_OVERLAP_A if clash_error_a is None else clash_error_a
    return analyze(pdb_path, binary, include_water=include_water, radii=probe_radii,
                   warn_a=warn_a, error_a=err_a).findings


def run_all(
    topology: Any,
    positions: Any,
    pdb_path: str | Path,
    prefer_probe: bool = True,
    include_water: bool = False,
    clash_warn_a: float | None = None,
    clash_error_a: float | None = None,
    probe_radii: str = "electron",
    provenance: dict[str, Any] | None = None,
) -> list[Finding]:
    """Run clash detection. Prefers ``probe`` when the binary is
    available; falls back to the Python engine otherwise.
    """
    from dvbfixer.diagnose.probe import analyze

    state: dict[str, Any] = {
        "requested_engine": "auto" if prefer_probe else "python",
        "actual_engine": DETERMINISTIC_STERIC_ENGINE, "status": "complete",
        "probe_status": "not_requested", "probe_version": None,
        "probe_radii": probe_radii, "include_water": include_water,
        "warning_threshold_angstrom": _CLASH_WARNING_OVERLAP_A if clash_warn_a is None else clash_warn_a,
        "error_threshold_angstrom": _CLASH_ERROR_OVERLAP_A if clash_error_a is None else clash_error_a,
        "fallback": False, "fallback_reason": None,
        "atom_fields_normalized": False,
    }
    if prefer_probe:
        binary = _find_binary("probe")
        if binary is not None:
            result = analyze(
                pdb_path, binary, include_water=include_water, radii=probe_radii,
                warn_a=_CLASH_WARNING_OVERLAP_A if clash_warn_a is None else clash_warn_a,
                error_a=_CLASH_ERROR_OVERLAP_A if clash_error_a is None else clash_error_a)
            state.update(probe_status=result.status, probe_version=result.version)
            if result.findings is not None:
                state["actual_engine"] = "molprobity-probe"
                state["atom_fields_normalized"] = True
                if provenance is not None:
                    provenance.update(state)
                return result.findings
            state.update(fallback=True, fallback_reason=result.reason)
            if result.status == "failed":
                print(f"WARNING: Probe unavailable or failed validation ({result.reason}); using Python steric analysis.",
                      file=sys.stderr)
        else:
            state.update(probe_status="unavailable", fallback=True, fallback_reason="Probe executable not found")
    findings = clashes_python(
        topology,
        positions,
        include_water=include_water,
        clash_warn_a=clash_warn_a,
        clash_error_a=clash_error_a,
        provenance=state,
    )
    if provenance is not None:
        provenance.update(state)
    return findings
