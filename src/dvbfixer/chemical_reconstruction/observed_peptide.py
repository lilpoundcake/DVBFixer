"""Full observed-PDB transport for conditional local heavy probes; no output writer."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import PdbAtom, _component_ref, _parse_atom_lines
from dvbfixer.pdbutils.sequence_evidence import observed_peptide_continuity
from dvbfixer.protein_completion.links import parse_declared_links

from .contracts import (
    ChemistryKind,
    Component,
    Connection,
    Constraint,
    ConstraintKind,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from .peptide_state import LocalPeptideAuthority, PeptideStateError
from .plan_compiler import compile_chain_plan
from .proposals import AuditPolicy, GeometryPlan, Placement, PlacementKind


def _parse_source(source: bytes) -> tuple[list[str], list[PdbAtom]]:
    if not isinstance(source, bytes) or not 0 < len(source) <= 16*1024*1024:
        raise PeptideStateError("bounded nonempty original PDB bytes required")
    lines = source.decode("ascii").splitlines(keepends=True)
    if any(not line.endswith("\n") or "\r" in line for line in lines):
        raise PeptideStateError("LF source records required")
    for line in lines:
        if line[:6].strip() in {"ATOM", "HETATM"} and not line.startswith(("ATOM  ", "HETATM")):
            raise PeptideStateError("malformed structural prefix")
    atoms = _parse_atom_lines(lines)
    if not atoms:
        raise PeptideStateError("no source atoms")
    return lines, atoms


def select_interior_target(source: bytes, residue_name: str, *, chain: str, number: str | None = None) -> ComponentInstanceRef:
    """File-order identity selection only, never geometric success selection."""
    _, atoms = _parse_source(source)
    groups: dict[ComponentInstanceRef, str] = {}
    for atom in atoms:
        if atom.record == "ATOM  " and atom.chain_id == chain:
            groups.setdefault(_component_ref(atom), atom.residue_name)
    ordered = list(groups)
    for identity in ordered[1:-1]:
        if groups[identity] == residue_name and (number is None or identity.sequence_number == number):
            return identity
    raise PeptideStateError("declared interior target not found")


@dataclass(frozen=True)
class ObservedPeptideProbe:
    request: ReconstructionRequest
    source_atom_count: int
    withheld: frozenset[ExactAtomRef]
    real_peptide_partners: tuple[ExactAtomRef, ExactAtomRef]
    raw_link_record_count: int
    unassessed_internal_connections: tuple[tuple[ExactAtomRef, ExactAtomRef], ...]


def prepare_observed_probe(source: bytes, authority: LocalPeptideAuthority, *, withheld: frozenset[str],
                           operation: Operation = Operation.HEAVY,
                           minimum_volume: float = .05, planar_tolerance_angstrom: float = .05) -> ObservedPeptideProbe:
    if operation != Operation.HEAVY:
        raise PeptideStateError("amide-hydrogen-geometry-authority-unqualified")
    lines, atoms = _parse_source(source)
    target = authority.identity
    selected = [a for a in atoms if _component_ref(a) == target]
    if not selected or any(a.record != "ATOM  " or a.residue_name != authority.parent.component.component_id for a in selected):
        raise PeptideStateError("exact target/source state mismatch")
    if any(a.element in {"H", "D"} for a in selected):
        raise PeptideStateError("target-measured-hydrogen-geometry-unqualified")
    present = {a.name for a in selected}
    if not isinstance(withheld, frozenset) or not withheld <= present or withheld & {"N", "CA", "C", "O"}:
        raise PeptideStateError("explicit present sidechain-only withhold mask required")
    withheld_refs = frozenset(ExactAtomRef(target, n) for n in withheld)
    groups = list(dict.fromkeys(_component_ref(a) for a in atoms if a.record == "ATOM  " and a.chain_id == target.chain_id))
    index = groups.index(target)
    if index == 0 or index == len(groups)-1:
        raise PeptideStateError("both real observed peptide neighbours required")
    before, after = groups[index-1], groups[index+1]
    continuity = observed_peptide_continuity(lines, target.chain_id,
                                           [(i.sequence_number, i.insertion_code) for i in (before, target, after)])
    if continuity != (True, True):
        raise PeptideStateError("observed peptide ports are broken or ambiguous")
    incoming = (ExactAtomRef(before, "C"), ExactAtomRef(target, "N"))
    outgoing = (ExactAtomRef(target, "C"), ExactAtomRef(after, "N"))
    links: dict[frozenset[ExactAtomRef], Connection] = {
        frozenset(pair): Connection(ExternalLink(*pair, 1.), LinkKind.COVALENT) for pair in (incoming, outgoing)}
    original_refs = {a.serial: ExactAtomRef(_component_ref(a), a.name) for a in atoms}
    intra = []
    target_bonds = {frozenset((b.atom1, b.atom2)) for b in authority.state.graph.bonds}
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        fields = [line[i:i+5].strip() for i in range(6, len(line.rstrip("\n")), 5)]
        serials = [int(value) for value in fields if value]
        if not serials or any(s not in original_refs for s in serials):
            raise PeptideStateError("dangling or empty source CONECT")
        a = original_refs[serials[0]]
        for serial in serials[1:]:
            b = original_refs[serial]
            if a == b:
                raise PeptideStateError("self source CONECT")
            if a.component == b.component:
                if a.component == target and frozenset((a.name, b.name)) not in target_bonds:
                    raise PeptideStateError("source target CONECT contradicts declared state")
                if a.component != target:
                    intra.append((a, b))
            else:
                links.setdefault(frozenset((a, b)), Connection(ExternalLink(a, b), LinkKind.UNKNOWN))
    for declared in parse_declared_links(lines, atoms):
        pair = frozenset((declared.left, declared.right))
        if pair in links and links[pair].kind == LinkKind.COVALENT:
            continue
        links[pair] = (Connection(ExternalLink(declared.left, declared.right, 1.), LinkKind.COVALENT)
                       if declared.record == "SSBOND" else Connection(ExternalLink(declared.left, declared.right), LinkKind.UNKNOWN))
    observations = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM,
                                     float(lines[a.line_index][54:60])) for a in atoms
                         if ExactAtomRef(_component_ref(a), a.name) not in withheld_refs)
    expected = {a.name for a in authority.state.graph.atoms if a.element not in {"H", "D"}}
    retained_target = {o.identity.name for o in observations if o.identity.component == target}
    generated = frozenset(ExactAtomRef(target, name) for name in expected-retained_target)
    ids = tuple(dict.fromkeys(o.identity.component for o in observations))
    components = tuple(Component(i, ChemistryKind.PROTEIN, (authority.state,)) if i == target
                       else Component(i, ChemistryKind.UNKNOWN) for i in ids)
    raw_records = tuple(line for line in lines if line[:6].strip() in {"CONECT", "LINK", "SSBOND"})
    constraints = (*authority.constraints(minimum_volume=minimum_volume, planar_tolerance_angstrom=planar_tolerance_angstrom),
                   *(Constraint(ConstraintKind.DISTANCE, pair, Unit.ANGSTROM, 1., 1.8) for pair in (incoming, outgoing)))
    request = ReconstructionRequest(hashlib.sha256(source).hexdigest(), Operation.HEAVY, components, observations,
                                    frozenset({target}), Ownership(generated), tuple(links.values()), constraints,
                                    origin=Origin("conditional-local-peptide-observed-v1", authority.digest,
                                    (("withheld", digest_record(withheld_refs)), ("raw_link_records", digest_record(raw_records)),
                                     ("chemistry_scope", "ALA witnesses are not actual source neighbours; full context unassessed"))))
    return ObservedPeptideProbe(request, len(atoms), withheld_refs, (incoming[0], outgoing[1]), len(raw_records), tuple(intra))


def compile_peptide_plan(probe: ObservedPeptideProbe, policy: AuditPolicy) -> GeometryPlan:
    return compile_peptide_request_plan(probe.request, policy)


def compile_peptide_request_plan(request: ReconstructionRequest, policy: AuditPolicy) -> GeometryPlan:
    """Shared deterministic sidechain plan, without inventing observed probe metadata."""
    target = next(iter(request.targets))
    cb = ExactAtomRef(target, "CB")
    initial: tuple[Placement, ...] = ()
    if cb in request.ownership.generated:
        initial = (Placement(cb, (ExactAtomRef(target, "CA"), ExactAtomRef(target, "N"), ExactAtomRef(target, "C")), PlacementKind.CONE),)
    return compile_chain_plan(request, policy, initial_placements=initial)
