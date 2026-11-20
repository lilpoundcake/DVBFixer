"""Graph-cut fragment torsions and independent all-heavy-path trace measurement."""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.rotamers import periodic_degrees

from .ccd_authority import _molecule
from .contracts import (
    ContractError,
    Observation,
    Operation,
    ReconstructionRequest,
    _identity,
    _text,
    _tuple,
    digest_record,
)
from .measurement import torsion_degrees
from .peptide_state import LocalPeptideAuthority
from .proposals import Candidate


class ConformerError(ContractError):
    """Unsupported graph/motion or invalid transport; no silent fallback."""


@dataclass(frozen=True)
class FragmentTorsion:
    anchors: tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef, ExactAtomRef]
    moving: tuple[ExactAtomRef, ...]
    root_distance: int

    def __post_init__(self) -> None:
        _tuple(self.anchors, ExactAtomRef, "torsion anchors", maximum=4)
        _tuple(self.moving, ExactAtomRef, "moving fragment", maximum=4096)
        if (len(self.anchors) != 4 or len(set(self.anchors)) != 4 or not self.moving
                or len(set(self.moving)) != len(self.moving) or self.moving != tuple(sorted(self.moving))
                or set(self.anchors[1:3]) & set(self.moving) or self.anchors[3] not in self.moving
                or type(self.root_distance) is not int or not 0 <= self.root_distance <= 4096):
            raise ConformerError("invalid-bounded-fragment-torsion")
        for a in (*self.anchors, *self.moving):
            _identity(a)
            if a.component != self.anchors[0].component:
                raise ConformerError("fragment-torsion-cross-component")


@dataclass(frozen=True)
class ExcludedAxis:
    atoms: tuple[ExactAtomRef, ExactAtomRef]
    reason: str

    def __post_init__(self) -> None:
        _tuple(self.atoms, ExactAtomRef, "excluded axis", maximum=2)
        if len(self.atoms) != 2 or self.atoms[0] == self.atoms[1]:
            raise ConformerError("invalid-excluded-axis")
        for a in self.atoms:
            _identity(a)
        _text(self.reason, "excluded-axis reason")


@dataclass(frozen=True)
class TorsionTopology:
    roots: tuple[ExactAtomRef, ...]
    rotors: tuple[FragmentTorsion, ...]
    excluded: tuple[ExcludedAxis, ...]

    def __post_init__(self) -> None:
        _tuple(self.roots, ExactAtomRef, "torsion roots", maximum=64)
        _tuple(self.rotors, FragmentTorsion, "torsion axes", maximum=64)
        _tuple(self.excluded, ExcludedAxis, "excluded axes", maximum=8192)
        if (len(self.roots) < 3 or len(set(self.roots)) != len(self.roots) or self.roots != tuple(sorted(self.roots))
                or len({r.anchors[1:3] for r in self.rotors}) != len(self.rotors)):
            raise ConformerError("invalid-bounded-torsion-topology")
        for a in self.roots:
            _identity(a)

    @property
    def digest(self) -> str:
        return digest_record(("graph-cut-heavy-torsion-topology-v1", self))


def _heavy_graph(request: ReconstructionRequest) -> dict[ExactAtomRef, set[ExactAtomRef]]:
    if len(request.targets) != 1:
        raise ConformerError("conformer-needs-one-target")
    target = next(c for c in request.components if c.identity in request.targets)
    if target.resolved is None:
        raise ConformerError("conformer-needs-resolved-target")
    heavy = {ExactAtomRef(target.identity, a.name) for a in target.resolved.graph.atoms if a.element not in {"H", "D"}}
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {a: set() for a in heavy}
    for bond in target.resolved.graph.bonds:
        a, b = ExactAtomRef(target.identity, bond.atom1), ExactAtomRef(target.identity, bond.atom2)
        if a in heavy and b in heavy:
            adjacent[a].add(b)
            adjacent[b].add(a)
    return adjacent


def _reachable(adjacent: dict[ExactAtomRef, set[ExactAtomRef]], start: ExactAtomRef,
                cut: frozenset[ExactAtomRef]) -> set[ExactAtomRef]:
    found, pending = {start}, [start]
    while pending:
        a = pending.pop()
        for b in adjacent[a]:
            if frozenset((a, b)) != cut and b not in found:
                found.add(b)
                pending.append(b)
    return found


def compile_fragment_torsions(request: ReconstructionRequest, authority: LocalPeptideAuthority,
                              roots: tuple[ExactAtomRef, ...]) -> TorsionTopology:
    target = next((c for c in request.components if c.identity == authority.identity), None)
    if (request.operation != Operation.HEAVY or request.targets != {authority.identity} or target is None
            or target.resolved != authority.state or not request.ownership.generated
            or request.ownership.temporary or request.ownership.mobile):
        raise ConformerError("conformer-needs-bound-fixed-heavy-request")
    graph = authority.state.graph
    adjacent = _heavy_graph(request)
    observed = {o.identity: o for o in request.observations if o.identity.component == authority.identity}
    generated = request.ownership.generated
    if (set(observed) | generated != adjacent.keys() or set(observed) & generated
            or any(o.occupancy != 1 or o.element in {"H", "D"} for o in observed.values())):
        raise ConformerError("conformer-needs-complete-target-heavy-inventory")
    if (not isinstance(roots, tuple) or not 3 <= len(roots) <= 64 or len(set(roots)) != len(roots)
            or roots != tuple(sorted(roots)) or not set(roots) <= observed.keys()):
        raise ConformerError("conformer-needs-exact-observed-roots")
    if _reachable(adjacent, roots[0], frozenset()) != set(adjacent):
        raise ConformerError("conformer-needs-connected-heavy-graph")
    distances = {a: 0 for a in roots}
    pending = list(roots)
    while pending:
        a = pending.pop(0)
        for b in sorted(adjacent[a]):
            if b not in distances:
                distances[b] = distances[a]+1
                pending.append(b)
    flat = authority.witness.flat_fragments
    molecule = _molecule(flat, {a.name: a.ideal_position for a in flat.atoms}, extra_bonds=authority.witness.extra_bonds)
    mapping = dict(authority.witness.node_map)
    indices = {mapping[a.name]: i for i, a in enumerate(flat.atoms)}
    ideals = {ExactAtomRef(authority.identity, a.name): a.ideal_position for a in graph.atoms}
    rotors, excluded = [], []
    for bond in sorted(graph.bonds, key=lambda b: tuple(sorted((b.atom1, b.atom2)))):
        a, b = sorted((ExactAtomRef(authority.identity, bond.atom1), ExactAtomRef(authority.identity, bond.atom2)))
        if a not in adjacent or b not in adjacent:
            continue
        pair = frozenset((a, b))
        native = molecule.GetBondBetweenAtoms(indices[a], indices[b])
        hybrid = {str(molecule.GetAtomWithIdx(indices[x]).GetHybridization()) for x in (a, b)}
        reason = None
        if bond.order != 1 or bond.aromatic:
            reason = "not-nonaromatic-single-bond"
        elif native.GetIsConjugated() or not hybrid <= {"SP2", "SP3"} or "SP3" not in hybrid:
            reason = "conjugated-or-unqualified-hybridization"
        upstream = _reachable(adjacent, a, pair)
        if reason is None and b in upstream:
            reason = "cyclic-axis"
        if reason is not None:
            excluded.append(ExcludedAxis((a, b), reason))
            continue
        if set(roots) <= upstream:
            u, v = a, b
        elif set(roots) <= set(adjacent)-upstream:
            u, v, upstream = b, a, set(adjacent)-upstream
        else:
            excluded.append(ExcludedAxis((a, b), "roots-span-axis"))
            continue
        downstream = set(adjacent)-upstream
        moving = downstream-{v}
        if not moving:
            reason = "no-downstream-heavy-motion"
        elif not moving <= generated:
            reason = "observed-downstream-anchor"
        elif any({c.link.component_atom, c.link.other_atom} & (downstream & generated) for c in request.connections):
            reason = "generated-external-link"
        if reason is not None:
            excluded.append(ExcludedAxis((a, b), reason))
            continue
        frame = None
        for first in sorted(adjacent[u]-{v}):
            for last in sorted(adjacent[v] & moving):
                try:
                    torsion_degrees((ideals[first], ideals[u], ideals[v], ideals[last]))
                except ValueError:
                    continue
                frame = (first, u, v, last)
                break
            if frame is not None:
                break
        if frame is None:
            excluded.append(ExcludedAxis((a, b), "no-nondegenerate-source-frame"))
            continue
        rotors.append(FragmentTorsion(frame, tuple(sorted(moving)), distances[u]))
    return TorsionTopology(roots, tuple(sorted(rotors, key=lambda r: (r.root_distance, r.anchors[1:3]))), tuple(excluded))


def rotate_fragments(request: ReconstructionRequest, authority: LocalPeptideAuthority, baseline: Candidate,
                     topology: TorsionTopology, offsets: tuple[float, ...]) -> tuple[Observation, ...]:
    """RDKit proper dihedral owner; reject any hidden observation movement."""
    from rdkit.Chem import rdMolTransforms

    if len(offsets) != len(topology.rotors):
        raise ConformerError("conformer-offset-coverage")
    if not any(offsets):
        return baseline.additions
    points = {o.identity: o.angstrom for o in request.observations}
    points.update({o.identity: o.position for o in baseline.additions})
    flat = authority.witness.flat_fragments
    mapping = dict(authority.witness.node_map)
    positions = {a.name: points.get(mapping[a.name], a.ideal_position) for a in flat.atoms}
    molecule = _molecule(flat, positions, extra_bonds=authority.witness.extra_bonds)
    indices = {mapping[a.name]: i for i, a in enumerate(flat.atoms)}
    conf = molecule.GetConformer()
    target_heavy = _heavy_graph(request)
    observed = {o.identity for o in request.observations}

    def point(atom: ExactAtomRef) -> tuple[float, float, float]:
        p = conf.GetAtomPosition(indices[atom])
        return float(p.x), float(p.y), float(p.z)

    goals = [periodic_degrees(float(rdMolTransforms.GetDihedralDeg(conf, *(indices[a] for a in r.anchors)))+offset)
             for r, offset in zip(topology.rotors, offsets, strict=True)]
    for rotor, offset, goal in zip(topology.rotors, offsets, goals, strict=True):
        if offset == 0:
            continue
        before = {a: point(a) for a in target_heavy}
        rdMolTransforms.SetDihedralDeg(conf, *(indices[a] for a in rotor.anchors), goal)
        for a in target_heavy:
            after = point(a)
            if a in observed and after != before[a]:
                raise ConformerError("conformer-rotator-moved-observation")
            if a not in rotor.moving and math.dist(after, before[a]) > 1e-10:
                raise ConformerError("conformer-rotator-moved-outside-fragment")
    return tuple(replace(o, position=point(o.identity)) for o in baseline.additions)


@dataclass(frozen=True)
class TorsionTraceAudit:
    passed: bool
    findings: tuple[str, ...]
    checked_paths: int
    maximum_error_degrees: float


def audit_torsion_offsets(request: ReconstructionRequest, baseline: Candidate, candidate: Candidate,
                          topology: TorsionTopology, offsets: tuple[float, ...], *, rounded: bool) -> TorsionTraceAudit:
    """Measure ALL heavy four-paths; never call RDKit SetDihedral or a constructor."""
    if len(offsets) != len(topology.rotors):
        return TorsionTraceAudit(False, ("torsion-offset-coverage",), 0, 0.)
    observed = {o.identity: o.angstrom for o in request.observations}
    base = observed | {o.identity: o.position for o in baseline.additions}
    points = observed | {o.identity: o.position for o in candidate.additions}
    adjacent = _heavy_graph(request)
    delta = {frozenset(r.anchors[1:3]): value for r, value in zip(topology.rotors, offsets, strict=True)}
    paths: set[tuple[ExactAtomRef, ExactAtomRef, ExactAtomRef, ExactAtomRef]] = set()
    for b in adjacent:
        for c in adjacent[b]:
            for a in adjacent[b]-{c}:
                for d in adjacent[c]-{b, a}:
                    path = (a, b, c, d)
                    if set(path) & request.ownership.generated:
                        paths.add(min(path, (d, c, b, a)))
    errors, maximum = set(), 0.
    for path in sorted(paths):
        a, b, c, d = path
        try:
            original = torsion_degrees((base[a], base[b], base[c], base[d]))
            current = torsion_degrees((points[a], points[b], points[c], points[d]))
            error = abs(periodic_degrees(current-original-delta.get(frozenset(path[1:3]), 0.)))
        except (KeyError, ValueError):
            errors.add("torsion-unmeasurable-path")
            continue
        maximum = max(maximum, error)
        if error > (.25 if rounded else 1e-6):
            errors.add("torsion-offset-or-unsampled-bond-conflict")
    return TorsionTraceAudit(not errors, tuple(sorted(errors)), len(paths), maximum)
