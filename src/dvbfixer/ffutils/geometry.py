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


class ChiralityError(RuntimeError):
    """Raised when D-amino-acid Cα stereochemistry survives all repair
    attempts.

    ``residues`` is a list of ``(chain_id, resid, resname, triple)``
    tuples describing every residue that failed the L-configuration
    check. Consumers (``dvbfixer.zbs``, ``dvbfixer.protonate``,
    ``dvbfixer.minimize``) print this list and exit with a non-zero
    status rather than write a PDB that carries D-Cα geometry into
    downstream MD.
    """

    def __init__(self, residues: list[tuple[str, str, str, float]]) -> None:
        self.residues = residues
        lines = [f"{c}/{name}{rid}: triple={t:+.5f} nm³"
                 for (c, rid, name, t) in residues]
        super().__init__(
            f"D-amino-acid Cα chirality detected on {len(residues)} residue"
            f"{'s' if len(residues) != 1 else ''} after all repair passes:\n"
            + "\n".join(f"  {ln}" for ln in lines)
        )

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

# Sibling-H clash tolerance (nm). Two H atoms bonded to the same
# heavy parent should sit at proper sp3 spacing — 2·bond_len·sin(54.75°)
# ≈ 1.78 Å for methylene/methyl. Anything under 0.15 nm (1.5 Å) is
# a placement bug that produces astronomical LJ 1/r^12 energies at
# minimize startup, which minimize then "relieves" by flipping the
# parent's chirality. Detect and re-place tetrahedrally.
_SIBLING_H_MIN_NM = 0.15  # 1.5 Å

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

    A hydrogen is misplaced when ANY of:
    - distance to its bonded heavy-atom parent exceeds
      ``_PARENT_DIST_FACTOR`` * canonical covalent bond length, OR
    - it sits within ``_COINCIDENT_TOL_NM`` of any other atom in the
      same residue, OR
    - it sits within ``_SIBLING_H_MIN_NM`` of a sibling H (same heavy
      parent). Sibling-H spacing must be ≥ sp3 minimum (~1.78 Å for
      methylene at 1.09 Å bond); anything closer produces astronomical
      LJ energies at minimize startup which drives the parent's
      chirality flip during relaxation.

    Pure detection — no mutation. Used by:
    - :func:`repair_misplaced_hydrogens` (repairs each finding).
    - :mod:`dvbfixer.diagnose.structural` (report-only).
    """
    findings: list[MisplacedHydrogen] = []
    seen: set[int] = set()

    for res in topology.residues():
        atoms = list(res.atoms())
        for h_atom in atoms:
            if h_atom.index in seen:
                continue
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
                seen.add(h_atom.index)
                continue

            coincident_hit = False
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
                    seen.add(h_atom.index)
                    coincident_hit = True
                    break
            if coincident_hit:
                continue

            # Sibling-H clash: too close to another H bonded to the same
            # heavy parent. Threshold 0.15 nm (1.5 Å) — well below proper
            # sp3 (1.78 Å) but well above the 0.05 nm coincident-tol.
            # Flag both siblings so repair rebuilds the parent's H cloud
            # as one tetrahedral group.
            siblings = _hydrogens_bonded_to(topology, parent)
            for sib in siblings:
                if sib.index == h_atom.index:
                    continue
                d2s = _dist2(h_pos, _pos(positions, sib))
                if d2s < _SIBLING_H_MIN_NM * _SIBLING_H_MIN_NM:
                    findings.append(
                        MisplacedHydrogen(
                            h_atom=h_atom,
                            parent=parent,
                            parent_distance_nm=parent_distance_nm,
                            expected_bond_nm=bond_len,
                            reason="sibling_clash",
                            coincident_with=sib,
                        )
                    )
                    seen.add(h_atom.index)
                    # Also flag the sibling so the group is rebuilt.
                    if sib.index not in seen:
                        sib_bond_len = _BOND_LEN_NM.get(
                            parent.element.symbol, _DEFAULT_BOND_LEN_NM,
                        )
                        sib_d = _dist2(_pos(positions, sib), parent_pos) ** 0.5
                        findings.append(
                            MisplacedHydrogen(
                                h_atom=sib,
                                parent=parent,
                                parent_distance_nm=sib_d,
                                expected_bond_nm=sib_bond_len,
                                reason="sibling_clash",
                                coincident_with=h_atom,
                            )
                        )
                        seen.add(sib.index)
                    break

    return findings


def _hydrogens_bonded_to(topology: Any, parent: Any) -> list[Any]:
    """Return every H atom bonded to ``parent``."""
    result = []
    for b1, b2 in topology.bonds():
        if b1.index == parent.index and b2.element is not None \
                and b2.element.symbol == "H":
            result.append(b2)
        elif b2.index == parent.index and b1.element is not None \
                and b1.element.symbol == "H":
            result.append(b1)
    return result


def _heavy_neighbors_of(topology: Any, atom: Any) -> list[Any]:
    """Return every non-H atom bonded to ``atom``."""
    result = []
    for b1, b2 in topology.bonds():
        if b1.index == atom.index and b2.element is not None \
                and b2.element.symbol != "H":
            result.append(b2)
        elif b2.index == atom.index and b1.element is not None \
                and b1.element.symbol != "H":
            result.append(b1)
    return result


def _tetrahedral_h_positions(
    parent_pos: tuple[float, float, float],
    heavy_neighbor_positions: list[tuple[float, float, float]],
    bond_len: float,
    n_h: int,
) -> list[tuple[float, float, float]] | None:
    """Compute ideal H positions on ``parent`` given its heavy neighbours.

    Handles the geometry cases needed by the ``repair_misplaced_hydrogens``
    caller — every combination that arises in a standard AA:

    * (1 heavy, 3 H): methyl / NH3+ / SH3 — three H at 109.5° from the
      neighbour direction, 120° apart around it.
    * (1 heavy, 2 H): amide NH2 / -NH2 — two H at 120° apart in the plane
      containing parent & neighbour.
    * (1 heavy, 1 H): OH / SH — one H at 109.5° from the neighbour in an
      arbitrary perpendicular direction. The old "linear-anti (180°)"
      choice is chemically strained; 109.5° is sp3-correct.
    * (2 heavy, 2 H): methylene CH2 — two H symmetric about the plane of
      the two heavy neighbours (β=54.75° from the bisector to give
      H-C-H ≈ 109.5°).
    * (2 heavy, 1 H): aromatic CH — one H at the negative bisector of
      the two heavy neighbours.
    * (3 heavy, 1 H): sp3 CH like HA on the backbone Cα — one H at the
      fourth tetrahedral vertex (negative sum of heavy unit vectors).

    Returns ``None`` for unhandled combinations (caller falls back to the
    per-H "rotate 60°" heuristic).
    """
    import math

    n_heavy = len(heavy_neighbor_positions)
    if n_heavy == 0 or n_h < 1:
        return None

    def _unit(vec: tuple[float, float, float]) -> tuple[float, float, float] | None:
        n = math.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
        if n < 1e-9:
            return None
        return (vec[0] / n, vec[1] / n, vec[2] / n)

    heavy_dirs: list[tuple[float, float, float]] = []
    for h_pos in heavy_neighbor_positions:
        u = _unit((h_pos[0] - parent_pos[0], h_pos[1] - parent_pos[1],
                   h_pos[2] - parent_pos[2]))
        if u is None:
            return None
        heavy_dirs.append(u)

    def _perp_to(u: tuple[float, float, float]) -> tuple[float, float, float]:
        seed = (1.0, 0.0, 0.0) if abs(u[0]) < 0.9 else (0.0, 1.0, 0.0)
        dot = seed[0] * u[0] + seed[1] * u[1] + seed[2] * u[2]
        r = (seed[0] - dot * u[0], seed[1] - dot * u[1], seed[2] - dot * u[2])
        r = _unit(r)
        # `seed` is not parallel to `u` so `r` is well-defined.
        assert r is not None
        return r

    def _cross(a: tuple[float, float, float],
               b: tuple[float, float, float]) -> tuple[float, float, float]:
        return (a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    if n_heavy == 1 and n_h == 3:
        u = heavy_dirs[0]
        p1 = _perp_to(u)
        p2 = _cross(u, p1)
        # cos(109.5°) = -0.334 → component along u is negative (away from
        # the heavy neighbour), which is what we want.
        ca = math.cos(math.radians(109.5))
        sa = math.sin(math.radians(109.5))
        out = []
        for i in range(3):
            th = math.radians(120 * i)
            ct, st = math.cos(th), math.sin(th)
            hx = parent_pos[0] + bond_len * (ca * u[0] + sa * (ct * p1[0] + st * p2[0]))
            hy = parent_pos[1] + bond_len * (ca * u[1] + sa * (ct * p1[1] + st * p2[1]))
            hz = parent_pos[2] + bond_len * (ca * u[2] + sa * (ct * p1[2] + st * p2[2]))
            out.append((hx, hy, hz))
        return out

    if n_heavy == 1 and n_h == 2:
        u = heavy_dirs[0]
        p1 = _perp_to(u)
        # Two H's at 120° from the parent-neighbour direction, symmetric
        # about that axis. cos(120°) = -0.5 gives the -u component,
        # sin(120°) = √3/2 the perpendicular component.
        c = -0.5
        s = math.sqrt(3) / 2
        out = []
        for sign in (+1.0, -1.0):
            hx = parent_pos[0] + bond_len * (c * u[0] + sign * s * p1[0])
            hy = parent_pos[1] + bond_len * (c * u[1] + sign * s * p1[1])
            hz = parent_pos[2] + bond_len * (c * u[2] + sign * s * p1[2])
            out.append((hx, hy, hz))
        return out

    if n_heavy == 1 and n_h == 1:
        u = heavy_dirs[0]
        p1 = _perp_to(u)
        ca = math.cos(math.radians(109.5))
        sa = math.sin(math.radians(109.5))
        return [(parent_pos[0] + bond_len * (ca * u[0] + sa * p1[0]),
                 parent_pos[1] + bond_len * (ca * u[1] + sa * p1[1]),
                 parent_pos[2] + bond_len * (ca * u[2] + sa * p1[2]))]

    if n_heavy == 2 and n_h == 2:
        u1, u2 = heavy_dirs
        s_vec = (u1[0] + u2[0], u1[1] + u2[1], u1[2] + u2[2])
        s_n = math.sqrt(s_vec[0] ** 2 + s_vec[1] ** 2 + s_vec[2] ** 2)
        if s_n < 1e-6:
            return None
        bisector = (-s_vec[0] / s_n, -s_vec[1] / s_n, -s_vec[2] / s_n)
        perp = _cross(u1, u2)
        perp_u = _unit(perp)
        if perp_u is None:
            return None
        # β = 54.75° gives H-C-H = 109.5°.
        ca = math.cos(math.radians(54.75))
        sa = math.sin(math.radians(54.75))
        out = []
        for sign in (+1.0, -1.0):
            hx = parent_pos[0] + bond_len * (ca * bisector[0] + sign * sa * perp_u[0])
            hy = parent_pos[1] + bond_len * (ca * bisector[1] + sign * sa * perp_u[1])
            hz = parent_pos[2] + bond_len * (ca * bisector[2] + sign * sa * perp_u[2])
            out.append((hx, hy, hz))
        return out

    if n_heavy == 2 and n_h == 1:
        u1, u2 = heavy_dirs
        s_vec = (u1[0] + u2[0], u1[1] + u2[1], u1[2] + u2[2])
        s_n = math.sqrt(s_vec[0] ** 2 + s_vec[1] ** 2 + s_vec[2] ** 2)
        if s_n < 1e-6:
            return None
        return [(parent_pos[0] - bond_len * s_vec[0] / s_n,
                 parent_pos[1] - bond_len * s_vec[1] / s_n,
                 parent_pos[2] - bond_len * s_vec[2] / s_n)]

    if n_heavy == 3 and n_h == 1:
        s_vec = (heavy_dirs[0][0] + heavy_dirs[1][0] + heavy_dirs[2][0],
                 heavy_dirs[0][1] + heavy_dirs[1][1] + heavy_dirs[2][1],
                 heavy_dirs[0][2] + heavy_dirs[1][2] + heavy_dirs[2][2])
        s_n = math.sqrt(s_vec[0] ** 2 + s_vec[1] ** 2 + s_vec[2] ** 2)
        if s_n < 1e-6:
            return None
        return [(parent_pos[0] - bond_len * s_vec[0] / s_n,
                 parent_pos[1] - bond_len * s_vec[1] / s_n,
                 parent_pos[2] - bond_len * s_vec[2] / s_n)]

    return None


def repair_misplaced_hydrogens(
    topology: Any,
    positions: Any,
    verbose: bool = False,
) -> int:
    """Detect misplaced hydrogens (via :func:`detect_misplaced_hydrogens`)
    and re-place them at proper sp3 tetrahedral positions.

    When ``Modeller.addHydrogens`` occasionally places two or three H
    atoms on the same parent (methylene CH2, methyl / NH3+ / SH3) at
    identical coordinates, we can't repair each H independently to a
    single "linear-anti" target — they'd end up coincident again. This
    function groups misplaced findings by parent atom and calls
    :func:`_tetrahedral_h_positions` to compute a proper set of
    distinct sp3 positions, then assigns each misplaced H to a slot
    not already claimed by a well-placed sibling.

    Fallback (unusual heavy-neighbour count): rotate a linear-anti
    placement 60°·i around the parent-neighbour axis for the i-th
    misplaced H so at least the coincidence is broken.

    Downstream ``minimize`` relaxes the placements to the FF's own
    equilibrium — this function only needs to give it a non-degenerate
    starting point.

    Args:
        topology: OpenMM ``Topology``.
        positions: Modeller / Simulation positions (list-like of
            ``Vec3`` in nanometers). Mutated in place.
        verbose: Print each repair.

    Returns:
        Number of hydrogens repaired.
    """
    import math
    from collections import defaultdict

    from openmm import Vec3
    from openmm.unit import nanometer

    findings = detect_misplaced_hydrogens(topology, positions)
    if not findings:
        return 0

    by_parent: dict[int, list[MisplacedHydrogen]] = defaultdict(list)
    for f in findings:
        by_parent[f.parent.index].append(f)

    repairs = 0
    for parent_findings in by_parent.values():
        parent = parent_findings[0].parent
        bond_len = parent_findings[0].expected_bond_nm
        res = parent.residue
        all_h = _hydrogens_bonded_to(topology, parent)
        heavy_neighbors = _heavy_neighbors_of(topology, parent)
        parent_pos = _pos(positions, parent)

        h_targets = _tetrahedral_h_positions(
            parent_pos,
            [_pos(positions, hn) for hn in heavy_neighbors],
            bond_len,
            n_h=len(all_h),
        )

        if h_targets is not None:
            # Mark slots already occupied by well-placed sibling H's.
            misplaced_ids = {f.h_atom.index for f in parent_findings}
            used = set()
            for h in all_h:
                if h.index in misplaced_ids:
                    continue
                h_p = _pos(positions, h)
                best = min(range(len(h_targets)),
                           key=lambda i: _dist2(h_p, h_targets[i]))
                used.add(best)
            for finding in parent_findings:
                slot = next((i for i in range(len(h_targets)) if i not in used),
                            None)
                if slot is None:
                    if verbose:
                        print(f"  [geom] {res.chain.id}/{res.name}{res.id}/"
                              f"{finding.h_atom.name}: no tetrahedral slot "
                              f"left, skipping")
                    continue
                positions[finding.h_atom.index] = Vec3(*h_targets[slot]) * nanometer
                used.add(slot)
                repairs += 1
                if verbose:
                    print(f"  [geom] repaired {res.chain.id}/{res.name}{res.id}"
                          f"/{finding.h_atom.name} → tetrahedral slot "
                          f"{slot} @ {bond_len * 10:.2f} Å from {parent.name}")
            continue

        # Fallback: unhandled heavy-neighbour count → rotate a linear-anti
        # placement 60°·i around the parent-neighbour axis for each
        # subsequent misplaced H so at least the coincidence is broken.
        neighbor = _parents_other_heavy_neighbor(topology, parent,
                                                  parent_findings[0].h_atom)
        if neighbor is None:
            if verbose:
                print(f"  [geom] {res.chain.id}/{res.name}{res.id}/"
                      f"{parent.name}: no heavy neighbour — leaving H atom(s) "
                      f"as-is")
            continue
        nx, ny, nz = _pos(positions, neighbor)
        dx, dy, dz = parent_pos[0] - nx, parent_pos[1] - ny, parent_pos[2] - nz
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        if norm < 1e-9:
            continue
        ux, uy, uz = dx / norm, dy / norm, dz / norm
        # Any perpendicular seed vector.
        if abs(ux) < 0.9:
            sx, sy, sz = 1.0, 0.0, 0.0
        else:
            sx, sy, sz = 0.0, 1.0, 0.0
        dot = sx * ux + sy * uy + sz * uz
        sx, sy, sz = sx - dot * ux, sy - dot * uy, sz - dot * uz
        sn = math.sqrt(sx * sx + sy * sy + sz * sz)
        sx, sy, sz = sx / sn, sy / sn, sz / sn
        for i, finding in enumerate(parent_findings):
            theta = math.radians(60 * i)
            # Rodrigues-style rotation of the "linear-anti + small perp"
            # target around u by theta. For i=0, θ=0 → pure linear-anti.
            offset_scale = 0.08 if i > 0 else 0.0  # 0.8 Å tangential kick
            hx = parent_pos[0] + bond_len * ux + offset_scale * (
                math.cos(theta) * sx)
            hy = parent_pos[1] + bond_len * uy + offset_scale * (
                math.cos(theta) * sy)
            hz = parent_pos[2] + bond_len * uz + offset_scale * (
                math.cos(theta) * sz)
            positions[finding.h_atom.index] = Vec3(hx, hy, hz) * nanometer
            repairs += 1
            if verbose:
                print(f"  [geom] repaired {res.chain.id}/{res.name}{res.id}"
                      f"/{finding.h_atom.name} → fallback rotated {60 * i}°")

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

# Fixed backbone atom names that never move during chirality repair.
# Everything else in the residue is a sidechain atom and gets reflected
# alongside CB.
_BACKBONE_ATOM_NAMES = frozenset({
    "N", "CA", "C", "O",
    "H", "HA", "HA2", "HA3",   # amide H + Cα H (GLY has HA2/HA3)
    "H1", "H2", "H3",           # N-terminal H
    "OXT",                       # C-terminal OXT
})


def fix_ca_chirality(
    topology: Any,
    positions: Any,
    verbose: bool = False,
) -> int:
    """Detect D-amino acid Cα stereochemistry and reflect the sidechain
    back to L.

    Reflects every non-backbone atom in the residue (CB, CG, OG, ND1,
    etc.) through the plane spanned by (N - CA) and (C - CA). The
    backbone (N, CA, C, O, HA, H, terminal H1/H2/H3, OXT) is left in
    place. Reflecting the whole sidechain (rather than only CB)
    preserves the sidechain's internal geometry after the flip —
    crucial when the input already has SER OG, VAL CG1/CG2, etc.
    placed on the wrong face.

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
        # Deadband -1e-6 nm³ ≈ -1e-3 Å³ — well above float noise for
        # nm-scale vectors and ~0.3% of a typical L triple (~0.35 Å³),
        # so marginally D geometries no longer slip.
        if triple >= -1e-6:
            continue

        norm2 = cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2
        if norm2 < 1e-18:
            # Degenerate — N, CA, C colinear. Nothing sensible to do.
            continue
        # Near-degeneracy gate: |cross| / (|vn|·|vc|) is sin(angle
        # between N-CA and C-CA). If that's < 0.1 the plane is
        # ill-defined and reflecting would just flip noise. Skip.
        vn_norm2 = vn[0] ** 2 + vn[1] ** 2 + vn[2] ** 2
        vc_norm2 = vc[0] ** 2 + vc[1] ** 2 + vc[2] ** 2
        if norm2 < 0.01 * vn_norm2 * vc_norm2:
            continue

        # Reflect every sidechain atom (anything not in the fixed
        # backbone set) through the CA-N-C plane.
        # For a point P, P_new = P - 2 · ((P - CA) · n̂) · n̂
        #                     = P - (2 · dot(P-CA, cross) / |cross|²) · cross
        n_moved = 0
        for atom in res.atoms():
            if atom.name in _BACKBONE_ATOM_NAMES:
                continue
            p = _pos(positions, atom)
            vp = (p[0] - ca_p[0], p[1] - ca_p[1], p[2] - ca_p[2])
            dot_vp_cross = vp[0] * cross[0] + vp[1] * cross[1] + vp[2] * cross[2]
            proj = dot_vp_cross / norm2
            new_p = (
                p[0] - 2.0 * proj * cross[0],
                p[1] - 2.0 * proj * cross[1],
                p[2] - 2.0 * proj * cross[2],
            )
            positions[atom.index] = Vec3(*new_p) * nanometer
            n_moved += 1

        repairs += 1
        if verbose:
            print(f"  [geom] flipped Cα chirality on "
                  f"{res.chain.id}/{res.name}{res.id} "
                  f"(triple {triple:.5f} nm³ → reflected {n_moved} sidechain atoms)")
    return repairs


# Residue names skipped by both the detector and the restraint builder.
# Same set as ``fix_ca_chirality``: GLY has no CB; waters/ions have no
# Cα chirality to defend.
_CHIRALITY_SKIP_RESNAMES = frozenset({
    "GLY", "HOH", "WAT", "TIP3", "TIP4", "TIP5", "SOL", "SPC", "SPCE",
})


def _find_ncacb(res: Any) -> tuple[Any, Any, Any, Any] | None:
    """Return (N, CA, C, CB) atoms from ``res`` or None if any is missing."""
    n = ca = c = cb = None
    for a in res.atoms():
        if a.name == "N":
            n = a
        elif a.name == "CA":
            ca = a
        elif a.name == "C":
            c = a
        elif a.name == "CB":
            cb = a
    if not (n and ca and c and cb):
        return None
    return n, ca, c, cb


def _ca_triple(
    positions: Any, n: Any, ca: Any, c: Any, cb: Any,
) -> tuple[float, float, float, float]:
    """Return (triple, |cross|², |vn|², |vc|²) for the CA chirality test."""
    ca_p = _pos(positions, ca)
    vn = tuple(a - b for a, b in zip(_pos(positions, n), ca_p))
    vc = tuple(a - b for a, b in zip(_pos(positions, c), ca_p))
    vcb = tuple(a - b for a, b in zip(_pos(positions, cb), ca_p))
    cross = (
        vn[1] * vc[2] - vn[2] * vc[1],
        vn[2] * vc[0] - vn[0] * vc[2],
        vn[0] * vc[1] - vn[1] * vc[0],
    )
    triple = cross[0] * vcb[0] + cross[1] * vcb[1] + cross[2] * vcb[2]
    norm2 = cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2
    vn_norm2 = vn[0] ** 2 + vn[1] ** 2 + vn[2] ** 2
    vc_norm2 = vc[0] ** 2 + vc[1] ** 2 + vc[2] ** 2
    return triple, norm2, vn_norm2, vc_norm2


def find_d_residues(
    topology: Any, positions: Any,
) -> list[tuple[str, str, str, float]]:
    """Return every D-amino-acid residue in the topology.

    Uses the same triple-product test and deadband as
    :func:`fix_ca_chirality` plus the near-degeneracy gate, so the
    detector and the repair function agree on what counts as D.

    Each entry is ``(chain_id, resid, resname, triple_nm3)``.
    """
    offenders: list[tuple[str, str, str, float]] = []
    for res in topology.residues():
        if res.name in _CHIRALITY_SKIP_RESNAMES:
            continue
        found = _find_ncacb(res)
        if found is None:
            continue
        n, ca, c, cb = found
        triple, norm2, vn2, vc2 = _ca_triple(positions, n, ca, c, cb)
        if triple >= -1e-6:
            continue
        # Skip near-degenerate geometries the repair also skips.
        if norm2 < 1e-18 or norm2 < 0.01 * vn2 * vc2:
            continue
        offenders.append((res.chain.id, str(res.id), res.name, triple))
    return offenders


def assert_all_l(topology: Any, positions: Any) -> None:
    """Raise :class:`ChiralityError` if any residue is still D.

    Wrapper around :func:`find_d_residues`. Call after
    :func:`fix_ca_chirality` (and any subsequent OpenMM minimization
    pass) to guarantee no D geometry survives into the output.
    """
    offenders = find_d_residues(topology, positions)
    if offenders:
        raise ChiralityError(offenders)


# Matches PDBFixer's own internal ``addMissingAtoms()`` clash-escape
# threshold (``_findNearestDistance``'s 0.13 nm cutoff) — used to
# verify that escape mechanism actually worked, since it's an
# unseeded, fixed-iteration-budget Langevin walk and doesn't always
# succeed within that budget.
_CLASH_CUTOFF_NM = 0.13


def find_clashing_atoms(
    topology: Any,
    positions: Any,
    atom_indices: Any,
    cutoff_nm: float = _CLASH_CUTOFF_NM,
) -> list[int]:
    """Return the subset of ``atom_indices`` within ``cutoff_nm`` of any
    other atom in a DIFFERENT residue.

    Used right after ``PDBFixer.addMissingAtoms()`` to check whether
    its own internal clash-escape (unseeded Langevin dynamics,
    triggered when a freshly-rebuilt atom lands within this same 0.13
    nm cutoff of a neighbor) actually resolved the clash. That escape
    is a fixed-budget random walk — it can fail to fully separate the
    atoms, especially in a densely packed environment, and unlike a
    hard crash this failure is silent: the rebuilt residue is simply
    left overlapping its neighbor.

    ``atom_indices`` is any iterable of ``atom.index`` ints (the
    caller already knows which atoms were newly added — no new
    tracking needed). Exclusion is by WHOLE RESIDUE, not just direct
    bonds — matching PDBFixer's own exclusion set in this same check
    (``exclusions = {a.index for a in atom.residue.atoms()}``, plus
    any cross-residue bond partner). A flexible sidechain can legally
    place two of its own atoms (e.g. NZ and a mid-chain HG) within this
    cutoff in a gauche conformation; only an INTER-residue contact this
    close is a real clash.
    """
    atom_list = list(topology.atoms())
    same_residue: dict[int, set[int]] = {
        a.index: {other.index for other in a.residue.atoms()} for a in atom_list
    }
    for bond in topology.bonds():
        i1, i2 = bond.atom1.index, bond.atom2.index
        same_residue[i1].add(i2)
        same_residue[i2].add(i1)

    all_positions = [_pos(positions, a) for a in atom_list]
    cutoff2 = cutoff_nm * cutoff_nm

    clashing: list[int] = []
    for idx in atom_indices:
        p_idx = all_positions[idx]
        skip = same_residue[idx]
        for j, p_j in enumerate(all_positions):
            if j == idx or j in skip:
                continue
            if _dist2(p_idx, p_j) < cutoff2:
                clashing.append(idx)
                break
    return clashing


_DEFAULT_REBUILD_ATTEMPTS = 5


def rebuild_missing_atoms_with_retry(
    fixer: Any,
    verbose: bool = False,
    max_attempts: int = _DEFAULT_REBUILD_ATTEMPTS,
    log_prefix: str = "",
) -> None:
    """Call ``fixer.addMissingAtoms()`` with an explicit, retried seed
    until the rebuild is clean (L-chirality, clash-free) — instead of
    accepting whatever PDBFixer's unseeded default hands back.

    ``PDBFixer.addMissingAtoms(seed=None)`` rebuilds a missing sidechain
    via template-overlay + a short local minimization; if that result
    clashes with a neighbor (< 0.13 nm — PDBFixer's own
    ``_findNearestDistance`` cutoff, matched by :func:`find_clashing_atoms`
    above), it falls back to UNSEEDED Langevin dynamics (300 K, up to
    2000 steps) to kick the new atoms apart. That's genuine stochastic
    MD: the escaped conformation differs run to run on the exact same
    input, and the fixed step budget doesn't always fully resolve the
    clash or land on L chirality. Papering over a bad rebuild
    downstream (reflect + re-minimize) only ever reacts to it after the
    fact — this retries the rebuild itself, at the source.

    Call after ``fixer.findMissingResidues()`` + ``fixer.findMissingAtoms()``
    (and, if applicable, ``findNonstandardResidues()`` +
    ``replaceNonstandardResidues()``) — i.e. whenever ``fixer`` would
    otherwise be ready for a plain ``fixer.addMissingAtoms()`` call.
    Mutates ``fixer.topology``/``fixer.positions`` in place, left set to
    the first clean attempt found (or the last attempt, with a printed
    warning, if none passed within ``max_attempts``).

    PDBFixer builds an entirely new ``Topology`` object on every call
    (``_addAtomsToTopology`` matches residues against
    ``self.missingAtoms``/``self.missingResidues`` by object identity,
    not value) and only reassigns ``self.topology``/``self.positions``
    at the very end — the pre-call objects are never mutated in place —
    so snapshotting them once and restoring before each attempt is a
    correct, cheap way to retry from scratch without reconstructing
    PDBFixer or re-running ``findMissingAtoms()``.

    Covers BOTH ``fixer.missingAtoms`` (heavy atoms added to an existing
    residue) AND ``fixer.missingResidues`` (whole residues inserted via
    SEQRES-driven gap-filling) — ``addMissingAtoms()`` runs the identical
    template-overlay + clash-escape MD over atoms from both in one call,
    so a clash on a gap-filled residue is exactly as real a risk as one
    on a rebuilt sidechain. Rather than replicate PDBFixer's internal
    residue-insertion-order arithmetic to predict a new residue's future
    ``(chain, resid, icode)`` ahead of time, this identifies "new" atoms
    by DIFFERENCE: any atom in the post-call topology whose identity key
    didn't exist in the pre-call snapshot is new, whether it came from a
    ``missingAtoms`` heavy-atom addition or a ``missingResidues`` whole
    new residue.
    """
    if not fixer.missingAtoms and not fixer.missingResidues:
        # Nothing to rebuild that carries this specific clash-escape risk.
        fixer.addMissingAtoms(seed=1)
        return

    snap_topology = fixer.topology
    snap_positions = fixer.positions
    existing_atom_keys = {
        (res.chain.id, res.id, res.insertionCode, atom.name)
        for res in snap_topology.residues() for atom in res.atoms()
    }

    for attempt in range(1, max_attempts + 1):
        fixer.topology = snap_topology
        fixer.positions = snap_positions
        fixer.addMissingAtoms(seed=attempt)

        new_atoms = [
            a for a in fixer.topology.atoms()
            if (a.residue.chain.id, a.residue.id, a.residue.insertionCode, a.name)
            not in existing_atom_keys
        ]
        new_atom_indices = [a.index for a in new_atoms]
        new_residue_keys = {
            (a.residue.chain.id, a.residue.id, a.residue.insertionCode) for a in new_atoms
        }

        d_residues = {
            (c, r) for (c, r, _name, _t) in find_d_residues(fixer.topology, fixer.positions)
        }
        chirality_bad = any((c, r) in d_residues for (c, r, _ic) in new_residue_keys)

        clashing = find_clashing_atoms(fixer.topology, fixer.positions, new_atom_indices)

        if not chirality_bad and not clashing:
            if verbose and attempt > 1:
                print(f"  {log_prefix}addMissingAtoms: clean rebuild on attempt "
                      f"{attempt} (seed={attempt})")
            return

        reasons = []
        if chirality_bad:
            reasons.append("D-Cα chirality")
        if clashing:
            reasons.append(f"{len(clashing)} clashing new atom(s)")
        if attempt < max_attempts:
            if verbose:
                print(f"  {log_prefix}addMissingAtoms attempt {attempt} "
                      f"(seed={attempt}): {', '.join(reasons)} — retrying with "
                      f"a new seed")
        else:
            print(f"  WARNING: {log_prefix}addMissingAtoms attempt {attempt} "
                  f"(seed={attempt}): {', '.join(reasons)} — exhausted "
                  f"{max_attempts} seeded rebuild attempts, keeping the last "
                  f"one. Downstream chirality reflect + minimize restraints "
                  f"remain as a last-resort safety net.")


# CYS-family residue names — used by the disulfide-bond helpers below.
# Kept alongside them since all three genuinely share one concern (SG-SG
# bond identity across a topology rebuild), unlike the general chirality/
# clash helpers above.
CYS_FAMILY_RESNAMES = frozenset({"CYS", "CYX", "CYM"})


def collect_ss_pairs(topology: Any) -> set[tuple[tuple[str, str], tuple[str, str]]]:
    """Return the set of (chain, resid) pairs currently SG-SG bonded.

    Used to snapshot the true disulfide pairing before a PDBFixer rebuild
    or a fresh ``PDBFile`` load — both re-derive disulfides by a pure
    distance-cutoff scan (``Topology.createDisulfideBonds()``) with no 1:1
    matching, so tightly packed CYS clusters can pick up spurious extra
    SG-SG bonds. Restore the true pairing afterward via
    :func:`drop_spurious_inter_aa_bonds`'s ``valid_ss_pairs`` argument.
    """
    pairs = set()
    for b in topology.bonds():
        if (b[0].name == "SG" and b[1].name == "SG"
                and b[0].residue.name in CYS_FAMILY_RESNAMES
                and b[1].residue.name in CYS_FAMILY_RESNAMES):
            key = tuple(sorted([
                (b[0].residue.chain.id, b[0].residue.id),
                (b[1].residue.chain.id, b[1].residue.id),
            ]))
            pairs.add(key)
    return pairs


def drop_spurious_inter_aa_bonds(
    topology: Any,
    verbose: bool = False,
    valid_ss_pairs: set | None = None,
    positions: Any = None,
) -> int:
    """Remove spurious bonds that break FF template matching.

    Three classes of spurious bond, the first two stemming from
    over-eager CONECT inference (OpenBabel guessing bonds by distance
    on real X-ray coordinates):

    1. **Inter-residue non-peptide bonds** between two standard protein
       residues. Mirrors the filter in
       :func:`dvbfixer.pdbutils.inference._apply_filter`; the canonical
       C(prev) - N(next) peptide bond is kept, everything else dropped.
    2. **Hydrogens with more than one heavy-atom partner**. Every H can
       bond to exactly one atom; extra bonds violate valence and create
       the "1 C-O bond too many" template error.
    3. **Extra SG-SG bonds beyond the true disulfide pairing.** OpenMM's
       own ``PDBFile.__init__`` unconditionally calls
       ``Topology.createDisulfideBonds()`` on every load (and
       ``PDBFixer.addMissingAtoms()`` does the same internally on
       rebuild) — a pure distance-cutoff scan with no 1:1 matching, so
       tightly packed CYS clusters (e.g. two chain copies whose
       N-termini sit close together) get every pairwise SG-SG contact
       within cutoff, and template matching then fails with "N S
       atom(s) too many" (a disulfide SG must have exactly one external
       S bond). Resolved two ways: if ``valid_ss_pairs`` is given (see
       :func:`collect_ss_pairs`, snapshotted before a rebuild/load
       that's known to reintroduce this), only pairs in that set
       survive. Otherwise, if ``positions`` is given, resolve by greedy
       nearest-distance 1:1 matching.

    Returns the number of bonds dropped.
    """
    from dvbfixer.ffutils import PROTEIN_RESIDUES

    # Pass 1: classify inter-residue bonds. SG-SG candidates are
    # collected separately and resolved after the main pass (Pass 1.5)
    # since disambiguating them may need to compare across all of them,
    # not just decide bond-by-bond.
    all_bonds = list(topology.bonds())
    survivors = []
    dropped = []
    ss_candidates = []
    for bond in all_bonds:
        b1 = bond[0]
        b2 = bond[1]
        # Rule 2: H with a partner keeps its FIRST partner only.
        # (Deferred to Pass 2 for correct counting.)
        if b1.residue is b2.residue:
            survivors.append(bond)
            continue
        if (b1.residue.name not in PROTEIN_RESIDUES
                or b2.residue.name not in PROTEIN_RESIDUES):
            survivors.append(bond)
            continue
        # Canonical peptide bond between adjacent residues: C of prev
        # residue to N of next. Keep.
        if {b1.name, b2.name} == {"C", "N"}:
            survivors.append(bond)
            continue
        # Disulfide bond: SG-SG between two CYS-family residues (CYS/
        # CYX/CYM). text_rename_variants_to_parent rewrote CYX → CYS
        # before load, so both endpoints are CYS-named here even for
        # user-annotated CYX pairs. Resolved in Pass 1.5 below.
        if (b1.name == "SG" and b2.name == "SG"
                and b1.residue.name in CYS_FAMILY_RESNAMES
                and b2.residue.name in CYS_FAMILY_RESNAMES):
            ss_candidates.append(bond)
            continue
        dropped.append(bond)

    # Pass 1.5: resolve SG-SG candidates to at most one partner per atom.
    if len(ss_candidates) < 2:
        survivors.extend(ss_candidates)
    elif valid_ss_pairs is not None:
        for bond in ss_candidates:
            b1, b2 = bond[0], bond[1]
            key = tuple(sorted([
                (b1.residue.chain.id, b1.residue.id),
                (b2.residue.chain.id, b2.residue.id),
            ]))
            if key in valid_ss_pairs:
                survivors.append(bond)
            else:
                dropped.append(bond)
    elif positions is not None:
        from openmm.unit import nanometer as _nm

        def _pos(atom):
            return positions[atom.index].value_in_unit(_nm)

        def _dist2(bond):
            p1, p2 = _pos(bond[0]), _pos(bond[1])
            return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

        locked = set()
        for bond in sorted(ss_candidates, key=_dist2):
            b1, b2 = bond[0], bond[1]
            if b1.index in locked or b2.index in locked:
                dropped.append(bond)
                continue
            locked.add(b1.index)
            locked.add(b2.index)
            survivors.append(bond)
    else:
        # No way to disambiguate — keep everything (old behavior).
        survivors.extend(ss_candidates)

    # Pass 2: hydrogen valence check. An H atom can bond to exactly
    # one heavy atom; extra bonds break residue-template matching
    # ("1 C-O bond too many" and similar errors from over-eager
    # CONECT inference). Keep the FIRST heavy-atom partner seen for
    # each H; drop any subsequent partners.
    h_seen: dict[int, bool] = {}
    survivors_after_h = []
    for bond in survivors:
        b1, b2 = bond[0], bond[1]
        h_atom = None
        if b1.element is not None and b1.element.symbol == "H":
            h_atom = b1
        elif b2.element is not None and b2.element.symbol == "H":
            h_atom = b2
        if h_atom is not None:
            if h_seen.get(h_atom.index):
                dropped.append(bond)
                continue
            h_seen[h_atom.index] = True
        survivors_after_h.append(bond)

    if not dropped:
        return 0

    topology._bonds = survivors_after_h  # noqa: SLF001

    if verbose:
        for bond in dropped:
            b1, b2 = bond[0], bond[1]
            print(f"  [geom] dropped spurious bond "
                  f"{b1.residue.chain.id}/{b1.residue.name}{b1.residue.id}:{b1.name} "
                  f"- {b2.residue.chain.id}/{b2.residue.name}{b2.residue.id}:{b2.name}")
    return len(dropped)


# NOTE: An earlier revision added an ``improper_chirality_restraint``
# (``CustomTorsionForce`` on CA-N-C-CB with θ₀ = +34°, k = 1000 kJ/mol/rad²)
# to actively bias residues toward L during OpenMM ``minimizeEnergy``.
# It was removed because:
#   1. Field consensus (VMD/NAMD chirality plugin, pdb4amber, CHARMM-GUI,
#      AlphaFold amber-relax) is *reflect → plain minimize → verify* —
#      nobody drives D→L via a stiff runtime force.
#   2. L and D have identical FF energies in isolation; packing already
#      picks L in a folded protein. A stiff bias fights physics on the
#      rare cases where minimize genuinely prefers D (which is the wrong
#      moment to force convergence — the input needs fixing).
#   3. The Cartesian gradient of a dihedral has two singular denominators
#      (1/|b1×b2|², 1/|b2×b3|²). Any degeneracy in either produces NaN
#      forces, tripping ``minimizeEnergy``'s startup guard even on
#      inputs that would otherwise minimize fine.
# Reflection (``fix_ca_chirality``) plus the iterative post-minimize
# check-reflect-re-minimize loop in ``minimize/pipeline.py`` covers all
# realistic cases without the numerical fragility.
