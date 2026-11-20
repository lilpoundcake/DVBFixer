"""Private numeric H input with explicit generated-heavy, never source, provenance.

This contract proves transport and masks, not parent geometry or selection. The
staged service must independently audit those before accepting a final envelope.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from dvbfixer.domain.structure_identity import ExactAtomRef

from .contracts import (
    ContractError,
    Operation,
    Origin,
    Ownership,
    PreparedContextRequest,
    ReconstructionRequest,
    Unit,
    _sha,
    digest_record,
    validate_fractional_context,
)
from .proposals import Candidate, Strategy


def conditioning_snapshot_digest(original: ReconstructionRequest, checkpoint: Candidate) -> str:
    return digest_record(("rounded-heavy-conditioning-snapshot-v1", original.digest, checkpoint))


def conditioning_fractional_atoms(original: ReconstructionRequest) -> frozenset[ExactAtomRef]:
    if isinstance(original, PreparedContextRequest):
        return original.context_policy.fractional_atoms
    return frozenset()


def validate_heavy_conditioning(request: ReconstructionRequest, original: ReconstructionRequest, checkpoint: Candidate, *,
                                 strategy: Strategy, adapter: str,
                                 snapshot_digest: Callable[[ReconstructionRequest, Candidate], str]) -> None:
    """Pure shared transport validation, not parent/scientific acceptance."""
    # Exact types deny nested conditioning and foreign subclasses with hidden roles.
    if type(original) not in {ReconstructionRequest, PreparedContextRequest}:
        raise ContractError("conditioning-needs-original-source-request")
    if (original.operation != Operation.HEAVY or not original.ownership.generated
            or original.ownership.temporary or original.ownership.mobile
            or not isinstance(checkpoint, Candidate) or checkpoint.request_digest != original.digest
            or checkpoint.strategy != strategy):
        raise ContractError("conditioning-heavy-checkpoint-binding")
    if (request.components != original.components or request.targets != original.targets
            or request.connections != original.connections or request.ph != original.ph
            or request.constraints[:len(original.constraints)] != original.constraints):
        raise ContractError("conditioning-original-chemistry-or-constraints-changed")
    if (len(checkpoint.additions) != len(original.ownership.generated)
            or {o.identity for o in checkpoint.additions} != original.ownership.generated
            or original.ownership.generated & {o.identity for o in original.observations}):
        raise ContractError("conditioning-exact-heavy-mask-required")
    graphs = {ExactAtomRef(c.identity, a.name): a for c in original.components if c.resolved is not None
              for a in c.resolved.graph.atoms}
    if any(a.component not in original.targets or a not in graphs or graphs[a].element in {"H", "D"}
           for a in original.ownership.generated):
        raise ContractError("conditioning-needs-generated-canonical-heavy-elements")
    for observation in checkpoint.additions:
        if (observation.element != graphs[observation.identity].element or observation.unit != Unit.ANGSTROM
                or observation.occupancy != 1 or any(len(f"{x:8.3f}") != 8 or x != float(f"{x:8.3f}")
                                                    for x in observation.position)):
            raise ContractError("conditioning-needs-rounded-unit-heavy-checkpoint")
    if digest_record(request.observations) != digest_record(original.observations + checkpoint.additions):
        raise ContractError("conditioning-observation-transport-mismatch")
    snapshot = snapshot_digest(original, checkpoint)
    if (request.source_sha256 != snapshot or request.origin is None
            or request.origin.adapter != adapter or request.origin.request_digest != snapshot):
        raise ContractError("conditioning-snapshot-origin-binding")
    if request.origin.settings:
        keys = tuple(k for k, _ in request.origin.settings)
        declared = ("peptide_state_declaration", "declared_peptide_authority")
        complete = ("complete_H_policy", "complete_H_authority")
        if keys not in {declared, complete, declared + complete}:
            raise ContractError("conditioning-unknown-declarative-origin-settings")
        for key, value in request.origin.settings:
            if key == "complete_H_policy":
                if value not in {"complete-H-groups-v2-explicit-sp2-N-thiol", "complete-H-groups-v1-fixed-heavy"}:
                    raise ContractError("conditioning-unknown-complete-H-policy")
            else:
                _sha(value)
    target_observed = {o.identity for o in request.observations if o.identity.component in request.targets}
    heavy = {a for a, node in graphs.items() if a.component in request.targets and node.element not in {"H", "D"}}
    hydrogens = frozenset(a for a, node in graphs.items() if a.component in request.targets and node.element == "H")
    if target_observed != heavy:
        raise ContractError("conditioning-needs-complete-heavy-only-target")
    if request.ownership.mobile or request.ownership.temporary:
        raise ContractError("conditioning-cannot-move-or-own-temporary-atoms")
    if request.operation == Operation.HEAVY:
        if request.ownership.generated:
            raise ContractError("conditioning-heavy-mode-is-fixed-planning-only")
    elif request.operation != Operation.HYDROGENS or request.ownership.generated != hydrogens or not hydrogens:
        raise ContractError("conditioning-needs-exact-complete-H-mask")
    validate_fractional_context(request, conditioning_fractional_atoms(original))


@dataclass(frozen=True, kw_only=True)
class HeavyConditionedRequest(ReconstructionRequest):
    """Original observations plus immutable generated-heavy numeric conditioning."""
    original: ReconstructionRequest
    checkpoint: Candidate

    def __post_init__(self) -> None:
        super().__post_init__()
        validate_heavy_conditioning(self, self.original, self.checkpoint, strategy=Strategy.TORSION,
                                     adapter="generated-heavy-conditioning-v1", snapshot_digest=conditioning_snapshot_digest)

    @property
    def allowed_fractional_atoms(self) -> frozenset[ExactAtomRef]:
        return conditioning_fractional_atoms(self.original)

    @property
    def digest(self) -> str:
        return digest_record(("chemical-heavy-conditioned-request-v1", self))


def condition_on_heavy(original: ReconstructionRequest, checkpoint: Candidate) -> HeavyConditionedRequest:
    """Transport only; no coordinate generation, selection or scientific acceptance."""
    snapshot = conditioning_snapshot_digest(original, checkpoint)
    return HeavyConditionedRequest(
        snapshot, Operation.HEAVY, original.components, original.observations + checkpoint.additions,
        original.targets, Ownership(frozenset()), original.connections, original.constraints, original.ph,
        Origin("generated-heavy-conditioning-v1", snapshot), original=original, checkpoint=checkpoint,
    )
