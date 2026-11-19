"""Source-bound metadata preparation; never changes structural observations."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace

from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref
from dvbfixer.pdbutils.io import update_conect

from .contracts import (
    Observation,
    OpaqueContextPolicy,
    Origin,
    PreparedContextRequest,
    Unit,
    digest_record,
)
from .observed_peptide import ObservedPeptideProbe, _parse_source, prepare_observed_probe
from .peptide_state import LocalPeptideAuthority, PeptideStateError


@dataclass(frozen=True)
class RemovedMetadata:
    line_number: int
    original_line: str
    reason: str


def _project(original: str) -> tuple[str, tuple[RemovedMetadata, ...]]:
    lines, atoms = _parse_source(original.encode("ascii"))
    if any(a.alternate_location for a in atoms):
        raise PeptideStateError("explicit-altloc-selection-required")
    if len({a.model for a in atoms}) != 1:
        raise PeptideStateError("explicit-model-selection-required")
    serial_map = {a.serial: a.serial for a in atoms}
    output, removed = [], []
    for index, line in enumerate(lines, 1):
        if line[:6].strip() != "CONECT":
            output.append(line)
            continue
        if not line.startswith("CONECT"):
            raise PeptideStateError("malformed-source-CONECT-prefix")
        tail = line[6:].rstrip("\n")
        fields = [tail[i:i+5] for i in range(0, len(tail), 5)]
        if any(field.strip() and (len(field) != 5 or not field.strip().isascii() or not field.strip().isdigit()) for field in fields):
            raise PeptideStateError("malformed-source-CONECT-field")
        serials = tuple(int(field) for field in fields if field.strip())
        if serials and serials[0] in serials[1:]:
            raise PeptideStateError("self-source-CONECT")
        # Only the existing whole-record decision is used; valid bytes stay exact.
        dangling = update_conect(line, serial_map) is None
        if dangling or len(serials) < 2:
            removed.append(RemovedMetadata(index, line, "dangling-CONECT" if dangling else "empty-CONECT"))
        else:
            output.append(line)
    return "".join(output), tuple(removed)


@dataclass(frozen=True)
class SourcePreparation:
    original_text: str
    prepared_text: str
    removed_metadata: tuple[RemovedMetadata, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.original_text, str) or not isinstance(self.prepared_text, str):
            raise PeptideStateError("immutable-source-text-required")
        prepared, removed = _project(self.original_text)
        if self.prepared_text != prepared or self.removed_metadata != removed:
            raise PeptideStateError("source-preparation-projection-binding-mismatch")

    @property
    def original_sha256(self) -> str:
        return hashlib.sha256(self.original_text.encode("ascii")).hexdigest()

    @property
    def prepared_sha256(self) -> str:
        return hashlib.sha256(self.prepared_text.encode("ascii")).hexdigest()

    @property
    def digest(self) -> str:
        return digest_record(("observed-source-metadata-preparation-v1", self))


def prepare_source(source: bytes) -> SourcePreparation:
    # Validate/bound bytes before decoding, even for a no-op projection.
    _parse_source(source)
    original = source.decode("ascii")
    prepared, removed = _project(original)
    return SourcePreparation(original, prepared, removed)


def prepare_source_probe(preparation: SourcePreparation, authority: LocalPeptideAuthority, *,
                         withheld: frozenset[str]) -> ObservedPeptideProbe:
    source = preparation.prepared_text.encode("ascii")
    lines, atoms = _parse_source(source)
    # Check the complete target BEFORE withholding; a mask cannot hide disorder.
    for atom in atoms:
        if _component_ref(atom) == authority.identity and float(lines[atom.line_index][54:60]) != 1:
            raise PeptideStateError("target-partial-occupancy-unqualified")
    probe = prepare_observed_probe(source, authority, withheld=withheld)
    r = probe.request
    fractional = frozenset(o.identity for o in r.observations if o.occupancy != 1)
    policy = OpaqueContextPolicy(preparation.prepared_sha256, preparation.digest, fractional, digest_record(r.observations))
    origin = Origin("prepared-observed-peptide-v1", preparation.digest,
                    (("original_source_sha256", preparation.original_sha256),
                     ("prepared_source_sha256", preparation.prepared_sha256),
                     ("original_adapter_origin", digest_record(r.origin)),
                     ("context_policy", "retained-opaque-fractional-v1-conservative-clashes")))
    request = PreparedContextRequest(r.source_sha256, r.operation, r.components, r.observations,
                                      r.targets, r.ownership, r.connections, r.constraints, r.ph, origin,
                                      context_policy=policy)
    return replace(probe, request=request)


@dataclass(frozen=True)
class SourceTransportAudit:
    passed: bool
    findings: tuple[str, ...]


def audit_source_transport(preparation: SourcePreparation, probe: ObservedPeptideProbe) -> SourceTransportAudit:
    """Reparse source observations independently of the adapter and constructor."""
    request = probe.request
    errors = set()
    if not isinstance(request, PreparedContextRequest):
        return SourceTransportAudit(False, ("prepared-context-request-required",))
    if (request.source_sha256 != preparation.prepared_sha256
            or request.context_policy.preparation_digest != preparation.digest):
        errors.add("source-preparation-binding")
    lines, atoms = _parse_source(preparation.prepared_text.encode("ascii"))
    identities = {ExactAtomRef(_component_ref(a), a.name) for a in atoms}
    if not probe.withheld <= identities or any(a.component not in request.targets for a in probe.withheld):
        errors.add("source-withhold-identity")
    if any(a.name in {"N", "CA", "C", "O"} for a in probe.withheld):
        errors.add("source-withhold-backbone-unqualified")
    expected = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position, Unit.ANGSTROM,
                                 float(lines[a.line_index][54:60])) for a in atoms
                     if ExactAtomRef(_component_ref(a), a.name) not in probe.withheld)
    if digest_record(request.observations) != digest_record(expected):
        errors.add("source-observation-transport")
    if probe.source_atom_count != len(atoms):
        errors.add("source-atom-count")
    if any(_component_ref(a) in request.targets and float(lines[a.line_index][54:60]) != 1 for a in atoms):
        errors.add("target-partial-occupancy-unqualified")
    settings = dict(request.origin.settings) if request.origin is not None else {}
    if (settings.get("original_source_sha256") != preparation.original_sha256
            or settings.get("prepared_source_sha256") != preparation.prepared_sha256):
        errors.add("source-origin-binding")
    return SourceTransportAudit(not errors, tuple(sorted(errors)))
