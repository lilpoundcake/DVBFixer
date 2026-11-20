"""Pinned family-universal carbonyl type evidence, never charge/state assignment."""
from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.ffutils.variants import AMBER_VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS
from dvbfixer.rtp_parser import RTPResidue, parse_rtp

from .ccd_authority import CcdGeometryAuthority
from .contracts import (
    ChemistryKind,
    LinkKind,
    Observation,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    _identity,
    _number,
    _sha,
    _text,
    _tuple,
    canonical_graph_digest,
    digest_record,
)

RTP_SHA256 = "ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419"
VARIANT_OWNER_DIGEST = "f92fc4f12e8997e14efbd407703e1df34460dbbaa78f34e1ac137891ad678806"
CATALOG_PAYLOAD_DIGEST = "2f34f225fadd62b052ad1d8d74421af05cc94e7adace041c01fd9d8ba808b940"


def variant_owner_digest() -> str:
    return digest_record(tuple(sorted(AMBER_VARIANT_TO_PARENT.items())))


@dataclass(frozen=True)
class RtpTypeState:
    name: str
    atoms: tuple[tuple[str, str, float, int], ...]
    bonds: tuple[tuple[str, str], ...]
    impropers: tuple[tuple[str, ...], ...]

    def __post_init__(self) -> None:
        _text(self.name, "RTP state")
        _tuple(self.atoms, tuple, "RTP atom rows", maximum=256)
        _tuple(self.bonds, tuple, "RTP bond rows", maximum=512)
        _tuple(self.impropers, tuple, "RTP improper rows", maximum=256)
        if (not self.atoms or self.atoms != tuple(sorted(self.atoms))
                or len({row[0] for row in self.atoms if row}) != len(self.atoms)):
            raise ValueError("exact-unique-sorted-RTP-atoms-required")
        for row in self.atoms:
            if (len(row) != 4 or type(row[0]) is not str or type(row[1]) is not str or not row[0] or not row[1]
                    or not _number(row[2]) or type(row[3]) is not int):
                raise ValueError("malformed-RTP-atom-row")
        if self.bonds != tuple(sorted(set(self.bonds))) or any(len(row) != 2 or row != tuple(sorted(row)) for row in self.bonds):
            raise ValueError("exact-unique-normalized-RTP-bonds-required")
        if self.impropers != tuple(sorted(self.impropers)) or any(len(row) != 4 for row in self.impropers):
            raise ValueError("exact-RTP-improper-rows-required")


@dataclass(frozen=True)
class IncomingTypeCatalog:
    source: str
    source_sha256: str
    owner_digest: str
    states: tuple[RtpTypeState, ...]

    def __post_init__(self) -> None:
        _text(self.source, "RTP source")
        _tuple(self.states, RtpTypeState, "complete RTP catalog", maximum=64)
        payload = tuple((s.name, s.atoms, s.bonds, s.impropers) for s in self.states)
        if (self.source_sha256 != RTP_SHA256 or self.owner_digest != VARIANT_OWNER_DIGEST
                or variant_owner_digest() != VARIANT_OWNER_DIGEST
                or tuple(s.name for s in self.states) != tuple(sorted({s.name for s in self.states}))
                or digest_record(payload) != CATALOG_PAYLOAD_DIGEST):
            raise ValueError("pinned-complete-RTP-catalog-or-variant-owner-mismatch")


def load_incoming_type_catalog(path: Path) -> IncomingTypeCatalog:
    with path.open("rb") as handle:
        data = handle.read(4*1024*1024+1)
    if len(data) > 4*1024*1024 or hashlib.sha256(data).hexdigest() != RTP_SHA256:
        raise ValueError("incoming-type-RTP-source-pin")
    # Typed boundary around the existing legacy parser, not a second RTP parser.
    parse = cast(Callable[[Path], tuple[object, dict[str, RTPResidue]]], parse_rtp)
    _, entries = parse(path)
    names = sorted((set(CANONICAL_HEAVY_ATOMS) | set(AMBER_VARIANT_TO_PARENT)) & entries.keys())
    states = tuple(RtpTypeState(name,
        tuple(sorted((str(r[0]), str(r[1]), float(r[2]), int(r[3])) for r in entries[name].atoms)),
        tuple(sorted((min(str(r[0]), str(r[1])), max(str(r[0]), str(r[1]))) for r in entries[name].bonds)),
        tuple(sorted(tuple(str(x) for x in row) for row in entries[name].impropers))) for name in names)
    if hashlib.sha256(path.read_bytes()).hexdigest() != RTP_SHA256:
        raise ValueError("incoming-type-RTP-changed-during-read")
    return IncomingTypeCatalog(str(path), RTP_SHA256, variant_owner_digest(), states)


@dataclass(frozen=True)
class IncomingSourceBinding:
    original_source_sha256: str
    prepared_source_sha256: str
    preparation_digest: str
    request_digest: str
    observations_digest: str
    probe_digest: str
    incoming: ExactAtomRef
    recipient: ExactAtomRef
    residue_name: str
    core: tuple[Observation, ...]
    annotation_pairs: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]

    def __post_init__(self) -> None:
        for value in (self.original_source_sha256, self.prepared_source_sha256, self.preparation_digest,
                      self.request_digest, self.observations_digest, self.probe_digest):
            _sha(value)
        _identity(self.incoming)
        _identity(self.recipient)
        _text(self.residue_name, "literal source residue label")
        _tuple(self.core, Observation, "observed carbonyl core", maximum=4)
        _tuple(self.annotation_pairs, tuple, "carbonyl source annotations", maximum=16)
        if (len(self.core) != 4 or len({o.identity for o in self.core}) != 4
                or self.core != tuple(sorted(self.core, key=lambda o: o.identity))):
            raise ValueError("exact-four-source-core-records-required")
        if self.annotation_pairs != tuple(sorted(set(self.annotation_pairs))):
            raise ValueError("exact-unique-source-annotation-pairs-required")
        for pair in self.annotation_pairs:
            if len(pair) != 2 or pair != tuple(sorted(pair)) or pair[0] == pair[1]:
                raise ValueError("invalid-source-annotation-pair")
            for atom in pair:
                _identity(atom)


def candidate_family(label: str, catalog: IncomingTypeCatalog) -> tuple[str, tuple[str, ...]]:
    catalog.__post_init__()
    states = {s.name for s in catalog.states}
    if label in AMBER_VARIANT_TO_PARENT:
        if label not in states:
            raise ValueError("explicit-incoming-state-template-missing")
        return AMBER_VARIANT_TO_PARENT[label], (label,)
    if label not in CANONICAL_HEAVY_ATOMS:
        raise ValueError("incoming-label-namespace-or-family-unqualified")
    variants = {name for name, parent in AMBER_VARIANT_TO_PARENT.items() if parent == label}
    if not variants <= states:
        raise ValueError("incomplete-incoming-variant-family")
    family = variants | ({label} if label in states else set())
    if not family:
        raise ValueError("incoming-family-has-no-template")
    return label, tuple(sorted(family))


def invariant_carbonyl_type(states: tuple[RtpTypeState, ...]) -> str:
    """Local topology/type calculation; never compare or assign partial charges."""
    types = set()
    for state in states:
        atoms = {row[0]: row for row in state.atoms}
        neighbours = {b if a == "C" else a for a, b in state.bonds if "C" in {a, b}}
        if not {"CA", "C", "O"} <= atoms.keys() or neighbours != {"CA", "O"} or ("CA", "+N", "C", "O") not in state.impropers:
            raise ValueError("incompatible-incoming-carbonyl-RTP-topology")
        types.add(atoms["C"][1])
    if not states or len(types) != 1 or next(iter(types)) != "C":
        raise ValueError("incoming-carbonyl-type-not-family-invariant-C")
    return "C"


def _ccd_carbonyl(ccd: CcdGeometryAuthority, parent: str) -> None:
    canonical_graph_digest(ccd.state.graph)
    if ccd.component.component_id != parent or "PEPTIDE" not in ccd.component.component_type.upper():
        raise ValueError("incoming-canonical-CCD-parent-mismatch")
    atoms = {a.name: a for a in ccd.state.graph.atoms}
    required = {"C": "C", "CA": "C", "O": "O", "OXT": "O", "HXT": "H"}
    if any(name not in atoms or atoms[name].element != element for name, element in required.items()):
        raise ValueError("incoming-CCD-carbonyl-footprint-missing")
    if any(atoms[name].aromatic or atoms[name].formal_charge != 0 for name in ("C", "O", "CA")):
        raise ValueError("incoming-CCD-retained-core-not-neutral-nonaromatic")
    leaving = set(ccd.leaving_atoms)
    if not {"OXT", "HXT"} <= leaving or leaving & {"C", "O", "CA"}:
        raise ValueError("incoming-CCD-leaving-branch-unproven")
    incident = {b.atom2 if b.atom1 == "C" else b.atom1: (b.order, b.aromatic)
                for b in ccd.state.graph.bonds if "C" in {b.atom1, b.atom2}}
    if incident != {"CA": (1., False), "O": (2., False), "OXT": (1., False)}:
        raise ValueError("incoming-CCD-carbonyl-bond-orders-unproven")
    if not any({b.atom1, b.atom2} == {"OXT", "HXT"} and b.order == 1 and not b.aromatic for b in ccd.state.graph.bonds):
        raise ValueError("incoming-CCD-leaving-H-bond-unproven")


def _proof_values(request: ReconstructionRequest, source: IncomingSourceBinding, catalog: IncomingTypeCatalog,
                   ccd: CcdGeometryAuthority) -> tuple[str, tuple[str, ...], str]:
    if (type(request) is not PreparedContextRequest or source.request_digest != request.digest
            or source.prepared_source_sha256 != request.source_sha256
            or source.preparation_digest != request.context_policy.preparation_digest
            or source.observations_digest != digest_record(request.observations)):
        raise ValueError("incoming-proof-source-request-binding")
    settings = dict(request.origin.settings) if request.origin else {}
    if settings.get("original_source_sha256") != source.original_source_sha256:
        raise ValueError("incoming-proof-original-source-binding")
    c, n = source.incoming, source.recipient
    ca, o = ExactAtomRef(c.component, "CA"), ExactAtomRef(c.component, "O")
    core = {c, ca, o, n}
    observed = {p.identity: p for p in request.observations}
    donor = next((component for component in request.components if component.identity == c.component), None)
    active = request.ownership.generated | request.ownership.mobile | request.ownership.temporary
    if (c.name != "C" or n.name != "N" or c.component == n.component or n.component not in request.targets
            or c.component in request.targets or donor is None or donor.kind != ChemistryKind.UNKNOWN
            or donor.states or donor.selected_state is not None or not core <= observed.keys() or core & active):
        raise ValueError("incoming-proof-needs-opaque-fixed-donor-and-directed-target-port")
    if any(observed[a].occupancy != 1 or observed[a].element != element or observed[a].unit != Unit.ANGSTROM
           for a, element in ((c, "C"), (ca, "C"), (o, "O"), (n, "N"))):
        raise ValueError("incoming-proof-core-element-or-occupancy")
    if digest_record(source.core) != digest_record(tuple(observed[a] for a in sorted(core))):
        raise ValueError("incoming-proof-core-observation-transport")
    if any(a.component == c.component and a.name in {"OXT", "HXT"} for a in observed):
        raise ValueError("incoming-proof-terminal-leaving-atoms-present")
    links = tuple(link for link in request.connections if c in {link.link.component_atom, link.link.other_atom})
    if not links or any(link.kind != LinkKind.COVALENT or link.link.order != 1
                        or {link.link.component_atom, link.link.other_atom} != {c, n} for link in links):
        raise ValueError("incoming-proof-carbonyl-external-link-unqualified")
    allowed = {frozenset((c, other)) for other in (ca, o, n)}
    if any(c not in pair or frozenset(pair) not in allowed for pair in source.annotation_pairs):
        raise ValueError("incoming-proof-contradictory-carbonyl-attachment")
    parent, names = candidate_family(source.residue_name, catalog)
    _ccd_carbonyl(ccd, parent)
    by_name = {s.name: s for s in catalog.states}
    common = invariant_carbonyl_type(tuple(by_name[name] for name in names))
    return parent, names, common


@dataclass(frozen=True)
class IncomingCarbonylTypeProof:
    request: ReconstructionRequest
    source: IncomingSourceBinding
    catalog: IncomingTypeCatalog
    ccd: CcdGeometryAuthority
    canonical_parent: str
    candidate_states: tuple[str, ...]
    atom_type: str
    selected_state: Literal[None] = None
    assigned_charge: Literal[None] = None
    scope: str = "local carbonyl-C type only; donor state/charge/H/sidechain and full polymer unassigned"

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        expected = _proof_values(self.request, self.source, self.catalog, self.ccd)
        if (expected != (self.canonical_parent, self.candidate_states, self.atom_type)
                or self.selected_state is not None or self.assigned_charge is not None
                or self.scope != "local carbonyl-C type only; donor state/charge/H/sidechain and full polymer unassigned"):
            raise ValueError("incoming-type-proof-family-type-or-scope-binding")

    @property
    def digest(self) -> str:
        return digest_record(("state-invariant-incoming-carbonyl-type-v1", self))


def prove_incoming_type(request: ReconstructionRequest, source: IncomingSourceBinding, catalog: IncomingTypeCatalog,
                         ccd: CcdGeometryAuthority) -> IncomingCarbonylTypeProof:
    parent, states, atom_type = _proof_values(request, source, catalog, ccd)
    return IncomingCarbonylTypeProof(request, source, catalog, ccd, parent, states, atom_type)
