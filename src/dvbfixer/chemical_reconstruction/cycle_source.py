"""Source-name transport for private cycle audits, not chemical state inference."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Literal

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.pdb_component_reconstruction import _component_ref

from .contracts import ChemistryKind, _sha, _text, _tuple, digest_record
from .cycle_contract import CycleConstructionProblem
from .observed_peptide import _parse_source
from .template_comparison import masked_probe_source, validate_masked_source


def _names(masked: bytes) -> tuple[tuple[ComponentInstanceRef, str], ...]:
    _, atoms = _parse_source(masked)
    names: dict[ComponentInstanceRef, str] = {}
    for atom in atoms:
        component = _component_ref(atom)
        if component in names and names[component] != atom.residue_name:
            raise ValueError("cycle-source-conflicting-residue-names")
        names[component] = atom.residue_name
    return tuple(sorted(names.items()))


@dataclass(frozen=True)
class CycleSourceMetadata:
    problem_digest: str
    request_digest: str
    source_sha256: str
    observations_digest: str
    masked_sha256: str
    removed_records_digest: str
    masked_pdb: str
    residue_names: tuple[tuple[ComponentInstanceRef, str], ...]
    revision: Literal["cycle-source-metadata-v1"] = "cycle-source-metadata-v1"

    def __post_init__(self) -> None:
        for pin in (self.problem_digest, self.request_digest, self.source_sha256, self.observations_digest,
                    self.masked_sha256, self.removed_records_digest):
            _sha(pin)
        _tuple(self.residue_names, tuple, "cycle source residue names", maximum=100000)
        if not isinstance(self.masked_pdb, str) or not 0 < len(self.masked_pdb) <= 16*1024*1024:
            raise ValueError("cycle-source-masked-byte-budget")
        masked = self.masked_pdb.encode("ascii")
        if (self.revision != "cycle-source-metadata-v1" or hashlib.sha256(masked).hexdigest() != self.masked_sha256
                or digest_record(_names(masked)) != digest_record(self.residue_names)):
            raise ValueError("cycle-source-metadata-byte-name-binding")
        for _, name in self.residue_names:
            _text(name, "source residue name")

    @property
    def digest(self) -> str:
        return digest_record(("cycle-source-metadata-v1", self))

    @property
    def source_coverage_verified(self) -> Literal[False]:
        return False


def bind_cycle_source(problem: CycleConstructionProblem, prepared: bytes) -> CycleSourceMetadata:
    """I/O admission adapter; never sends complete withheld records to construction.

    The original/prepared bytes must be reopened independently for saved readback.
    A self-consistent masked-name DTO alone cannot authenticate its source label.
    Only a digest of removed records, NOT their native coordinates, is retained.
    """
    if not isinstance(problem, CycleConstructionProblem) or type(prepared) is not bytes:
        raise ValueError("cycle-source-typed-problem-and-bytes-required")
    pin = hashlib.sha256(prepared).hexdigest()
    if pin != problem.source.source_sha256:
        raise ValueError("cycle-source-original-byte-binding")
    masked, removed = masked_probe_source(prepared, problem.source)
    names = _names(masked)
    lookup = dict(names)
    for component in problem.source.components:
        if component.identity in problem.source.targets and component.kind == ChemistryKind.PROTEIN and component.resolved is not None:
            source_name = lookup[component.identity]
            declared_name = component.resolved.graph.component_id
            if VARIANT_TO_PARENT.get(source_name, source_name) != VARIANT_TO_PARENT.get(declared_name, declared_name):
                raise ValueError("cycle-source-target-parent-name-mismatch")
    return CycleSourceMetadata(problem.digest, problem.source.digest, pin, problem.coverage.observations_digest,
        hashlib.sha256(masked).hexdigest(), digest_record(removed), masked.decode("ascii"), names)


def audit_cycle_source(problem: CycleConstructionProblem, prepared: bytes, metadata: CycleSourceMetadata) -> None:
    """Verify actual source bytes/names/masks, not a self-asserted digest cache."""
    if not isinstance(metadata, CycleSourceMetadata):
        raise ValueError("cycle-source-typed-metadata-required")
    validate_masked_source(metadata.masked_pdb.encode("ascii"), problem.source)
    expected = bind_cycle_source(problem, prepared)
    if digest_record(expected) != digest_record(metadata):
        raise ValueError("cycle-source-independent-readback-binding")
