"""Explicitly typed real tertiary peptide-N geometry; not atom-type inference."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from dvbfixer.domain.structure_identity import ExactAtomRef

from .amide_authority import FFBONDED_SHA256, _select
from .contracts import (
    Constraint,
    ConstraintKind,
    LinkKind,
    ReconstructionRequest,
    Unit,
    _identity,
    _sha,
    digest_record,
)
from .declared_peptide import DeclaredPeptideAuthority

PORT_ROWS_DIGEST = "a6da05b22f952873169cdd1603c596d3b990f56f46ebd15070aedcb028eaafa6"


@dataclass(frozen=True)
class PlanarPortAuthority:
    source: str
    source_sha256: str
    rows: tuple[tuple[str, ...], ...]

    def __post_init__(self) -> None:
        if (not self.source or self.source_sha256 != FFBONDED_SHA256
                or digest_record(self.rows) != PORT_ROWS_DIGEST):
            raise ValueError("exact-pinned-planar-port-rows-required")

    @property
    def digest(self) -> str:
        return digest_record(("ff19sb-tertiary-amide-port-v1", self))


def load_planar_port_authority(path: Path) -> PlanarPortAuthority:
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != FFBONDED_SHA256:
        raise ValueError("planar-port-parameter-source-pin")
    text = data.decode("ascii")
    rows = (_select(text, "bondtypes", ("C", "N")), _select(text, "angletypes", ("C", "N", "CT")),
            _select(text, "angletypes", ("C", "N", "XC")), _select(text, "angletypes", ("CT", "N", "XC")))
    return PlanarPortAuthority(str(path), FFBONDED_SHA256, rows)


@dataclass(frozen=True)
class PlanarPortDeclaration:
    incoming: ExactAtomRef
    nitrogen: ExactAtomRef
    alpha: ExactAtomRef
    distal: ExactAtomRef
    atom_types: tuple[str, str, str, str]
    peptide_authority_digest: str
    parameter_authority_digest: str

    def __post_init__(self) -> None:
        for atom in (self.incoming, self.nitrogen, self.alpha, self.distal):
            _identity(atom)
        _sha(self.peptide_authority_digest)
        _sha(self.parameter_authority_digest)
        if (not isinstance(self.atom_types, tuple) or len(self.atom_types) != 4
                or any(type(t) is not str or not t for t in self.atom_types)):
            raise ValueError("explicit-four-port-atom-types-required")

    @property
    def digest(self) -> str:
        return digest_record(("explicit-real-tertiary-peptide-port-v1", self))


def planar_port_constraints(original: ReconstructionRequest, authority: DeclaredPeptideAuthority,
                             declaration: PlanarPortDeclaration, parameters: PlanarPortAuthority) -> tuple[Constraint, ...]:
    """Use real fixed port positions and declared exact types; never witness XYZ."""
    incoming, n, ca, cd = declaration.incoming, declaration.nitrogen, declaration.alpha, declaration.distal
    if (declaration.peptide_authority_digest != authority.digest or declaration.parameter_authority_digest != parameters.digest
            or declaration.atom_types != ("C", "N", "XC", "CT")
            or (incoming.name, n.name, ca.name, cd.name) != ("C", "N", "CA", "CD")
            or {n.component, ca.component, cd.component} != {authority.identity}
            or incoming.component == authority.identity or original.targets != {authority.identity}):
        raise ValueError("planar-port-exact-directed-type-and-authority-declaration")
    component = next((c for c in original.components if c.identity == authority.identity), None)
    if component is None or component.resolved != authority.state:
        raise ValueError("planar-port-target-state-binding")
    graph = authority.state.graph
    atoms = {a.name: a for a in graph.atoms}
    neighbours = {b.atom2 if b.atom1 == "N" else b.atom1 for b in graph.bonds if "N" in {b.atom1, b.atom2}}
    if (neighbours != {"CA", "CD"} or atoms["N"].element != "N" or atoms["N"].formal_charge != 0
            or atoms["N"].aromatic or any(atoms[name].element != "C" for name in neighbours)
            or any(b.order != 1 or b.aromatic for b in graph.bonds if "N" in {b.atom1, b.atom2})):
        raise ValueError("planar-port-needs-neutral-two-carbon-no-H-nitrogen")
    observed = {o.identity: o for o in original.observations}
    if (not {incoming, n, ca} <= observed.keys() or cd not in original.ownership.generated
            or any(observed[a].occupancy != 1 or observed[a].element != element
                   for a, element in ((incoming, "C"), (n, "N"), (ca, "C")))
            or not any(c.kind == LinkKind.COVALENT and c.link.order == 1
                       and {c.link.component_atom, c.link.other_atom} == {incoming, n} for c in original.connections)):
        raise ValueError("planar-port-needs-real-fixed-incoming-link-and-generated-distal")
    label = "real-planar-peptide-port:" + declaration.digest
    length = 10 * float(parameters.rows[0][-2])
    angles = tuple(float(row[-2]) for row in parameters.rows[1:])
    return (Constraint(ConstraintKind.DISTANCE, (incoming, n), Unit.ANGSTROM, length-.10, length+.10, label=label),
            *(Constraint(ConstraintKind.ANGLE, refs, Unit.DEGREE, value-12., value+12., label=label)
              for refs, value in zip(((incoming, n, cd), (incoming, n, ca), (cd, n, ca)), angles, strict=True)),
            Constraint(ConstraintKind.PLANARITY, (n, incoming, ca, cd), Unit.ANGSTROM, upper=.05, label=label))
