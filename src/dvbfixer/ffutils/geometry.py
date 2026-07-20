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


def repair_misplaced_hydrogens(
    topology: Any,
    positions: Any,
    verbose: bool = False,
) -> int:
    """Detect and repair hydrogens that ``addHydrogens`` placed badly.

    A hydrogen is considered misplaced when EITHER:

    - its distance to its bonded heavy-atom parent exceeds
      ``_PARENT_DIST_FACTOR`` * expected covalent bond length, OR
    - it sits within ``_COINCIDENT_TOL_NM`` of any other atom in the
      same residue.

    Each misplaced H is re-placed at ``bond_length`` from its parent,
    in the direction opposite the parent's other heavy neighbor
    (linear-anti). This is a chemically coarse but topologically sane
    placement — downstream minimize relaxes it to proper sp3 geometry.

    Args:
        topology: OpenMM ``Topology``.
        positions: Modeller / Simulation positions (a list-like of
            ``Vec3`` in nanometers). Mutated in place.
        verbose: Print each repair.

    Returns:
        Number of hydrogens repaired.
    """
    from openmm import Vec3
    from openmm.unit import nanometer

    def _pos(atom: Any) -> tuple[float, float, float]:
        p = positions[atom.index].value_in_unit(nanometer)
        return float(p[0]), float(p[1]), float(p[2])

    def _dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
        return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2

    repairs = 0
    residues = list(topology.residues())
    for res in residues:
        atoms = list(res.atoms())
        for h_atom in atoms:
            if h_atom.element is None or h_atom.element.symbol != "H":
                continue

            parent = _bonded_heavy_parent(topology, h_atom)
            if parent is None:
                # Unbonded H — nothing we can do sensibly. Skip.
                continue

            bond_len = _BOND_LEN_NM.get(parent.element.symbol, _DEFAULT_BOND_LEN_NM)
            parent_pos = _pos(parent)
            h_pos = _pos(h_atom)
            d2_parent = _dist2(h_pos, parent_pos)

            # Detection: too far from parent OR coincident with any other
            # atom in the same residue.
            misplaced = d2_parent > (bond_len * _PARENT_DIST_FACTOR) ** 2
            coincident_with: Any | None = None
            if not misplaced:
                for other in atoms:
                    if other.index == h_atom.index:
                        continue
                    d2 = _dist2(h_pos, _pos(other))
                    if d2 < _COINCIDENT_TOL_NM * _COINCIDENT_TOL_NM:
                        coincident_with = other
                        misplaced = True
                        break

            if not misplaced:
                continue

            # Repair: place H at bond_len from parent, anti to the
            # parent's other heavy neighbor.
            neighbor = _parents_other_heavy_neighbor(topology, parent, h_atom)
            if neighbor is None:
                if verbose:
                    print(f"  [geom] {res.chain.id}/{res.name}{res.id}/{h_atom.name}: "
                          f"misplaced but parent {parent.name} has no other heavy "
                          f"neighbor — leaving as-is")
                continue

            nx, ny, nz = _pos(neighbor)
            dx = parent_pos[0] - nx
            dy = parent_pos[1] - ny
            dz = parent_pos[2] - nz
            norm = (dx * dx + dy * dy + dz * dz) ** 0.5
            if norm < 1e-9:
                if verbose:
                    print(f"  [geom] {res.chain.id}/{res.name}{res.id}/{h_atom.name}: "
                          f"parent and neighbor coincident — leaving as-is")
                continue

            scale = bond_len / norm
            new_h = (
                parent_pos[0] + dx * scale,
                parent_pos[1] + dy * scale,
                parent_pos[2] + dz * scale,
            )
            positions[h_atom.index] = Vec3(*new_h) * nanometer
            repairs += 1
            if verbose:
                reason = (
                    f"coincident with {coincident_with.name}"
                    if coincident_with is not None
                    else f"parent {parent.name} distance {d2_parent ** 0.5 * 10:.2f} Å > "
                         f"{bond_len * _PARENT_DIST_FACTOR * 10:.2f} Å"
                )
                print(f"  [geom] repaired {res.chain.id}/{res.name}{res.id}/{h_atom.name} "
                      f"({reason}) → linear-anti from {neighbor.name} at "
                      f"{bond_len * 10:.2f} Å from {parent.name}")

    return repairs
