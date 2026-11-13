"""Private local intervention policy; strict original proposer is unchanged."""
from __future__ import annotations

from dataclasses import dataclass, replace

from scripts.hydrogen_initialization_guard import (
    GuardRefusal,
    Input,
    Result,
    norm,
    normalized_volume,
    sub,
    validate,
)
from scripts.hydrogen_initialization_guard import (
    propose as strict_propose,
)

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

REVISION = 'private-compatible-noop-v1'


@dataclass(frozen=True)
class Warning:
    atom: ExactAtomRef
    residue_name: str
    radius_nm: float
    heavy_bond_lengths_nm: tuple[float, ...]
    reason: str = 'compatible_noop_outside_distance_admission_not_geometry_certified'


@dataclass(frozen=True)
class Proposal:
    result: Result
    original_new_hydrogens: frozenset[ExactAtomRef]
    intervention_hydrogens: frozenset[ExactAtomRef]
    warnings: tuple[Warning, ...]


def propose(value: Input) -> Proposal:
    adjacent = validate(value)
    supported, unsupported, interventions, warnings = [], [], [], []
    for h, atom in enumerate(value.atoms):
        if atom.identity not in value.new_hydrogens:
            continue
        parent = adjacent[h][0]
        center = value.atoms[parent]
        heavy = [i for i in adjacent[parent] if value.atoms[i].element not in ('H', 'D')]
        hs = [i for i in adjacent[parent] if value.atoms[i].element in ('H', 'D')]
        inventory = CANONICAL_HEAVY_ATOMS.get(VARIANT_TO_PARENT.get(center.residue_name, center.residue_name), ())
        if (center.element != 'C' or len(heavy) != 3 or hs != [h]
                or center.identity.name not in inventory
                or any(value.atoms[i].identity.component != center.identity.component for i in [parent, h]+heavy)
                or any(value.atoms[i].identity.name not in inventory for i in heavy)):
            unsupported.append(atom.identity)
            continue
        supported.append(atom.identity)
        heavy.sort(key=lambda i: value.atoms[i].identity)
        p = value.positions_nm[parent]
        neighbors = tuple(value.positions_nm[i] for i in heavy)
        hv = normalized_volume(neighbors, p)
        fv = normalized_volume(neighbors, value.positions_nm[h])
        if abs(hv) <= 1e-4 or abs(fv) <= 1e-4:
            raise GuardRefusal('degenerate supported center under local policy')
        if hv*fv < 0:
            interventions.append(atom.identity)
            continue
        radius = norm(sub(value.positions_nm[h], p))
        distances = tuple(norm(sub(q, p)) for q in neighbors)
        if not 0.05 <= radius <= 0.15 or any(not 0.10 <= d <= 0.20 for d in distances):
            warnings.append(Warning(atom.identity, atom.residue_name, radius, distances))
    mask = frozenset(interventions)
    # Complete atom/coordinate/graph inventory remains present for contact checks.
    strict = strict_propose(replace(value, new_hydrogens=mask))
    if frozenset(strict.moved) != mask:
        raise GuardRefusal('strict proposer did not resolve every intervention')
    result = Result(strict.positions_nm, strict.moved, tuple(supported), tuple(unsupported))
    return Proposal(result, value.new_hydrogens, mask, tuple(warnings))
