"""Source-byte admission for incoming type proofs; never forwards withheld points."""
from __future__ import annotations

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref

from .contracts import digest_record
from .declared_peptide import DeclaredPeptideAuthority
from .incoming_type import IncomingSourceBinding
from .observed_peptide import ObservedPeptideProbe, _parse_source
from .source_preparation import SourcePreparation, audit_source_transport, prepare_source_probe


def bind_incoming_source(preparation: SourcePreparation, probe: ObservedPeptideProbe,
                           authority: DeclaredPeptideAuthority) -> IncomingSourceBinding:
    if not audit_source_transport(preparation, probe).passed:
        raise ValueError("incoming-source-transport-refused")
    expected = prepare_source_probe(preparation, authority, withheld=frozenset(a.name for a in probe.withheld))
    if digest_record(expected) != digest_record(probe):
        raise ValueError("incoming-source-probe-or-annotation-transport")
    incoming = probe.real_peptide_partners[0]
    recipient = ExactAtomRef(authority.identity, "N")
    _, atoms = _parse_source(preparation.prepared_text.encode("ascii"))
    donor = tuple(a for a in atoms if _component_ref(a) == incoming.component)
    labels = {a.residue_name for a in donor}
    if len(labels) != 1 or any(a.record != "ATOM  " or a.alternate_location for a in donor):
        raise ValueError("incoming-source-exact-protein-residue-label-required")
    core_refs = {incoming, ExactAtomRef(incoming.component, "O"), ExactAtomRef(incoming.component, "CA"), recipient}
    core = tuple(sorted((o for o in probe.request.observations if o.identity in core_refs), key=lambda o: o.identity))
    annotations = {(min(a, b), max(a, b)) for a, b in probe.unassessed_internal_connections if incoming in {a, b}}
    annotations |= {(min(c.link.component_atom, c.link.other_atom), max(c.link.component_atom, c.link.other_atom)) for c in probe.request.connections
                    if incoming in {c.link.component_atom, c.link.other_atom}}
    return IncomingSourceBinding(preparation.original_sha256, preparation.prepared_sha256, preparation.digest,
        probe.request.digest, digest_record(probe.request.observations), digest_record(probe), incoming, recipient,
        next(iter(labels)), core, tuple(sorted(annotations)))


def audit_incoming_source(preparation: SourcePreparation, probe: ObservedPeptideProbe,
                           authority: DeclaredPeptideAuthority, source: IncomingSourceBinding) -> tuple[str, ...]:
    try:
        expected = bind_incoming_source(preparation, probe, authority)
    except ValueError as exc:
        return ("incoming-source-audit-refused:" + str(exc),)
    return () if digest_record(expected) == digest_record(source) else ("incoming-source-byte-slice-binding",)
