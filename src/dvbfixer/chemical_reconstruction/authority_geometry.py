"""Read-only length-unit views of explicit authority, not a chemistry inference."""
from __future__ import annotations

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import ReconstructionRequest, Unit
from .measurement import Point


def geometry_view(request: ReconstructionRequest) -> tuple[dict[ExactAtomRef, Point], dict[frozenset[ExactAtomRef], float]]:
    positions: dict[ExactAtomRef, Point] = {}
    bonds: dict[frozenset[ExactAtomRef], float] = {}
    for component in request.components:
        state = component.resolved
        if state is None:
            continue
        factor = 10. if state.ideal_unit == Unit.NANOMETER else 1.
        for atom in state.graph.atoms:
            x, y, z = atom.ideal_position
            positions[ExactAtomRef(component.identity, atom.name)] = x*factor, y*factor, z*factor
        for bond in state.graph.bonds:
            pair = frozenset((ExactAtomRef(component.identity, bond.atom1), ExactAtomRef(component.identity, bond.atom2)))
            bonds[pair] = bond.ideal_length*factor
    return positions, bonds
