"""Post-addHydrogens geometry sanity check.

Purpose: detect hydrogens that ``Modeller.addHydrogens`` placed at a
physically impossible position and re-place them with a coarse
linear-anti geometry. Downstream minimize relaxes the linear placement
to proper sp3 tetrahedral geometry, so this is sufficient for the
pipeline — the important property is that no H is coincident with
another atom and every H sits within a normal bond length of its
parent heavy atom.

The bug this catches
--------------------

Reported against a C-terminal SER whose input PDB has ``OXT`` present
but ``HG`` missing. OpenMM's ``addHydrogens`` (via the CSER template
path) then places the new ``HG`` right on top of ``OXT``. HG ends up
1.7 Å from its own OG (correct O-H is ~0.97 Å) and 0.001 Å from
OXT — the geometry is chemically nonsense and downstream tools
either produce garbage or crash. The pre-existing coincident-atom
strip in ``prepare/pipeline.py`` doesn't help because the failing
input has no coincident pair to detect — the coincidence is CREATED
by addHydrogens itself.

The fix
-------

After every ``modeller.addHydrogens(...)`` call, invoke
:func:`repair_misplaced_hydrogens` on the resulting topology. It
walks every H, verifies proximity to a bonded heavy-atom parent, and
repairs any violations in place.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Bond-length targets (nm). Matches AMBER14-standard covalent bond
# lengths for the given parent element. Used as both the tolerance
# for "misplaced" detection (with a 1.6x safety factor) and the
# replacement bond length when we repair.
_BOND_LEN_NM = {
    "O": 0.097,   # O-H
    "N": 0.101,   # N-H
    "C": 0.109,   # C-H
    "S": 0.134,   # S-H (CYS)
    "P": 0.142,   # P-H
}

# Coincident tolerance (nm). Any H within this of another atom in
# the same residue is considered misplaced regardless of its
# distance to its bonded parent.
_COINCIDENT_TOL_NM = 0.05  # 0.5 Å

# Misplaced-parent-distance factor: an H is misplaced if its
# distance to its bonded parent exceeds this multiple of the
# expected covalent bond length.
_PARENT_DIST_FACTOR = 1.6

# Fallback bond length when the parent element isn't in _BOND_LEN_NM.
_DEFAULT_BOND_LEN_NM = 0.1


def _bonded_heavy_parent(topology: Any, h_atom: Any) -> Any | None:
    """Return the heavy-atom parent of ``h_atom`` as recorded in the
    topology's bond list, or None if no bonded heavy neighbor exists.
    """
    idx = h_atom.index
    for b1, b2 in topology.bonds():
        if b1.index == idx and b2.element.symbol != "H":
            return b2
        if b2.index == idx and b1.element.symbol != "H":
            return b1
    return None


def _parents_other_heavy_neighbor(topology: Any, parent: Any, h_atom: Any) -> Any | None:
    """Return the parent atom's OTHER non-H bonded neighbor (ignoring
    the given ``h_atom``), used to fix the linear-anti direction for
    the replacement H position. Returns None if no such neighbor exists.
    """
    idx_p = parent.index
    idx_h = h_atom.index
    for b1, b2 in topology.bonds():
        if b1.index == idx_p and b2.element.symbol != "H" and b2.index != idx_h:
            return b2
        if b2.index == idx_p and b1.element.symbol != "H" and b1.index != idx_h:
            return b1
    return None


@dataclass
class MisplacedHydrogen:
    """Detection record for a hydrogen at a physically impossible position.

    Returned by :func:`detect_misplaced_hydrogens`. Consumed by
    :func:`repair_misplaced_hydrogens` (which re-places the H) and by
    :mod:`dvbfixer.diagnose.structural` (which reports without
    repairing).
    """

    h_atom: Any                 # OpenMM Atom
    parent: Any                 # OpenMM Atom — the H's bonded heavy neighbor
    parent_distance_nm: float   # actual d(H, parent) in nm
    expected_bond_nm: float     # canonical covalent bond length in nm
    reason: str                 # 'parent_distance' or 'coincident'
    coincident_with: Any | None  # OpenMM Atom, only set when reason == 'coincident'


@dataclass
class CoincidentAtoms:
    """Detection record for two atoms placed at the same position in the
    same residue.

    Returned by :func:`detect_coincident_atoms`. Consumed by
    :mod:`dvbfixer.prepare.pipeline` (which strips both atoms so
    PDBFixer's ``addMissingAtoms`` can re-place them cleanly) and by
    :mod:`dvbfixer.diagnose.structural` (report-only).
    """

    residue: Any                # OpenMM Residue
    hydrogen: Any               # OpenMM Atom, H-symbol
    heavy_atom: Any             # OpenMM Atom, coincident with the H
    distance_nm: float          # d(hydrogen, heavy_atom) in nm


def _pos(positions: Any, atom: Any) -> tuple[float, float, float]:
    """Extract atom coordinates in nm as a plain 3-tuple."""
    from openmm.unit import nanometer
    p = positions[atom.index].value_in_unit(nanometer)
    return float(p[0]), float(p[1]), float(p[2])


def _dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def detect_coincident_atoms(
    topology: Any,
    positions: Any,
    tol_nm: float = _COINCIDENT_TOL_NM,
) -> list[CoincidentAtoms]:
    """Return every (H, non-H atom) pair in the same residue placed
    within ``tol_nm`` of each other.

    Pure detection — no mutation. Used by:
    - ``dvbfixer.prepare.pipeline`` (before H-strip: strip the heavy
      partner too so PDBFixer's ``addMissingAtoms`` re-places it).
    - ``dvbfixer.diagnose.structural`` (report an ERROR finding).
    """
    findings: list[CoincidentAtoms] = []
    for res in topology.residues():
        atoms = list(res.atoms())
        for h_atom in atoms:
            if h_atom.element is None or h_atom.element.symbol != "H":
                continue
            pa = _pos(positions, h_atom)
            for other in atoms:
                if other.index == h_atom.index or other.element is None:
                    continue
                if other.element.symbol == "H":
                    continue
                pb = _pos(positions, other)
                d2 = _dist2(pa, pb)
                if d2 < tol_nm * tol_nm:
                    findings.append(
                        CoincidentAtoms(
                            residue=res,
                            hydrogen=h_atom,
                            heavy_atom=other,
                            distance_nm=d2 ** 0.5,
                        )
                    )
                    break
    return findings


def detect_misplaced_hydrogens(
    topology: Any,
    positions: Any,
) -> list[MisplacedHydrogen]:
    """Return every hydrogen whose position is chemically impossible.

    A hydrogen is misplaced when EITHER:
    - distance to its bonded heavy-atom parent exceeds
      ``_PARENT_DIST_FACTOR`` * canonical covalent bond length, OR
    - it sits within ``_COINCIDENT_TOL_NM`` of any other atom in the
      same residue.

    Pure detection — no mutation. Used by:
    - :func:`repair_misplaced_hydrogens` (repairs each finding).
    - :mod:`dvbfixer.diagnose.structural` (report-only).
    """
    findings: list[MisplacedHydrogen] = []
    for res in topology.residues():
        atoms = list(res.atoms())
        for h_atom in atoms:
            if h_atom.element is None or h_atom.element.symbol != "H":
                continue

            parent = _bonded_heavy_parent(topology, h_atom)
            if parent is None:
                continue

            bond_len = _BOND_LEN_NM.get(parent.element.symbol, _DEFAULT_BOND_LEN_NM)
            parent_pos = _pos(positions, parent)
            h_pos = _pos(positions, h_atom)
            d2_parent = _dist2(h_pos, parent_pos)
            parent_distance_nm = d2_parent ** 0.5

            if d2_parent > (bond_len * _PARENT_DIST_FACTOR) ** 2:
                findings.append(
                    MisplacedHydrogen(
                        h_atom=h_atom,
                        parent=parent,
                        parent_distance_nm=parent_distance_nm,
                        expected_bond_nm=bond_len,
                        reason="parent_distance",
                        coincident_with=None,
                    )
                )
                continue

            for other in atoms:
                if other.index == h_atom.index or other.element is None:
                    continue
                d2 = _dist2(h_pos, _pos(positions, other))
                if d2 < _COINCIDENT_TOL_NM * _COINCIDENT_TOL_NM:
                    findings.append(
                        MisplacedHydrogen(
                            h_atom=h_atom,
                            parent=parent,
                            parent_distance_nm=parent_distance_nm,
                            expected_bond_nm=bond_len,
                            reason="coincident",
                            coincident_with=other,
                        )
                    )
                    break

    return findings


def repair_misplaced_hydrogens(
    topology: Any,
    positions: Any,
    verbose: bool = False,
) -> int:
    """Detect misplaced hydrogens (via :func:`detect_misplaced_hydrogens`)
    and re-place each one at ``bond_length`` from its parent in the
    linear-anti direction from the parent's other heavy neighbor.

    Chemically coarse but topologically sane — downstream ``minimize``
    relaxes the anti placement to proper sp3.

    Args:
        topology: OpenMM ``Topology``.
        positions: Modeller / Simulation positions (list-like of
            ``Vec3`` in nanometers). Mutated in place.
        verbose: Print each repair.

    Returns:
        Number of hydrogens repaired.
    """
    from openmm import Vec3
    from openmm.unit import nanometer

    repairs = 0
    for finding in detect_misplaced_hydrogens(topology, positions):
        h_atom = finding.h_atom
        parent = finding.parent
        res = h_atom.residue

        # Place H at bond_len from parent, anti to the parent's other
        # heavy neighbor.
        neighbor = _parents_other_heavy_neighbor(topology, parent, h_atom)
        if neighbor is None:
            if verbose:
                print(f"  [geom] {res.chain.id}/{res.name}{res.id}/{h_atom.name}: "
                      f"misplaced but parent {parent.name} has no other heavy "
                      f"neighbor — leaving as-is")
            continue

        parent_pos = _pos(positions, parent)
        nx, ny, nz = _pos(positions, neighbor)
        dx = parent_pos[0] - nx
        dy = parent_pos[1] - ny
        dz = parent_pos[2] - nz
        norm = (dx * dx + dy * dy + dz * dz) ** 0.5
        if norm < 1e-9:
            if verbose:
                print(f"  [geom] {res.chain.id}/{res.name}{res.id}/{h_atom.name}: "
                      f"parent and neighbor coincident — leaving as-is")
            continue

        scale = finding.expected_bond_nm / norm
        new_h = (
            parent_pos[0] + dx * scale,
            parent_pos[1] + dy * scale,
            parent_pos[2] + dz * scale,
        )
        positions[h_atom.index] = Vec3(*new_h) * nanometer
        repairs += 1
        if verbose:
            reason = (
                f"coincident with {finding.coincident_with.name}"
                if finding.coincident_with is not None
                else f"parent {parent.name} distance "
                     f"{finding.parent_distance_nm * 10:.2f} Å > "
                     f"{finding.expected_bond_nm * _PARENT_DIST_FACTOR * 10:.2f} Å"
            )
            print(f"  [geom] repaired {res.chain.id}/{res.name}{res.id}/{h_atom.name} "
                  f"({reason}) → linear-anti from {neighbor.name} at "
                  f"{finding.expected_bond_nm * 10:.2f} Å from {parent.name}")

    return repairs


# ---------------------------------------------------------------------------
# Cα chirality repair
# ---------------------------------------------------------------------------
#
# ``PDBFixer.addMissingAtoms`` rebuilds missing sidechain heavy atoms
# via ideal AMBER template alignment. For residues where only the
# backbone (N, CA, C, O) was present in the input, the template
# alignment can pick the D configuration — placing CB on the wrong
# side of the Cα-N-C plane. Branched-Cβ residues (VAL, ILE, THR) are
# most vulnerable because their CB carries a chiral centre once the
# rest of the sidechain lands. Nothing in the pipeline detects this
# without ``fix_ca_chirality``.
#
# Detection: scalar triple product (N - CA) × (C - CA) · (CB - CA).
# Positive → L-amino acid (right-handed). Negative → D. GLY (no CB)
# is skipped.
#
# Repair: reflect CB through the plane spanned by (N - CA) and
# (C - CA). No other atoms are moved.

def fix_ca_chirality(
    topology: Any,
    positions: Any,
    verbose: bool = False,
) -> int:
    """Detect D-amino acid Cα stereochemistry and reflect CB back to L.

    Returns the number of residues repaired.
    """
    from openmm import Vec3
    from openmm.unit import nanometer

    _skip_names = frozenset({"GLY", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
                             "SOL", "SPC", "SPCE"})

    def _find(res: Any, name: str) -> Any | None:
        for a in res.atoms():
            if a.name == name:
                return a
        return None

    repairs = 0
    for res in topology.residues():
        if res.name in _skip_names:
            continue
        n = _find(res, "N")
        ca = _find(res, "CA")
        c = _find(res, "C")
        cb = _find(res, "CB")
        if not (n and ca and c and cb):
            continue

        ca_p = _pos(positions, ca)
        vn = tuple(a - b for a, b in zip(_pos(positions, n), ca_p))
        vc = tuple(a - b for a, b in zip(_pos(positions, c), ca_p))
        vcb = tuple(a - b for a, b in zip(_pos(positions, cb), ca_p))

        # Plane normal n̂ = (N-CA) × (C-CA) / |...|
        cross = (
            vn[1] * vc[2] - vn[2] * vc[1],
            vn[2] * vc[0] - vn[0] * vc[2],
            vn[0] * vc[1] - vn[1] * vc[0],
        )
        triple = cross[0] * vcb[0] + cross[1] * vcb[1] + cross[2] * vcb[2]
        # Deadband around zero avoids false positives on marginal
        # geometry (mirrors ``check_ca_chirality``'s threshold).
        if triple >= -1e-4:
            continue

        norm2 = cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2
        if norm2 < 1e-18:
            # Degenerate — N, CA, C colinear. Nothing sensible to do.
            continue
        # CB_new = CB - 2 · ((CB - CA) · n̂) · n̂ where n̂ = cross / |cross|
        # dot(vcb, n̂) · n̂ = (dot(vcb, cross) / norm2) · cross
        proj = triple / norm2  # = dot(vcb, cross) / |cross|²
        new_cb = (
            ca_p[0] + vcb[0] - 2.0 * proj * cross[0],
            ca_p[1] + vcb[1] - 2.0 * proj * cross[1],
            ca_p[2] + vcb[2] - 2.0 * proj * cross[2],
        )
        positions[cb.index] = Vec3(*new_cb) * nanometer
        repairs += 1
        if verbose:
            print(f"  [geom] flipped Cα chirality on "
                  f"{res.chain.id}/{res.name}{res.id} "
                  f"(triple {triple:.5f} nm³ → reflected CB)")
    return repairs
