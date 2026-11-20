"""One immutable rounded-heavy checkpoint followed by one complete-H proposal."""
from __future__ import annotations

from typing import Literal

from dvbfixer.domain.structure_identity import ExactAtomRef

from .amide_authority import AmideMetricAuthority
from .conformer_search import ConformerSearchResult
from .construction import propose
from .contracts import ReconstructionRequest, digest_record
from .declared_peptide import DeclaredPeptideAuthority
from .staged_audit import compile_staged_view, measure_staged
from .staged_contracts import StagedEnvelope, StagedResult, parent_digest


def complete_staged(original: ReconstructionRequest, parent: ConformerSearchResult,
                     authority: DeclaredPeptideAuthority, incoming: ExactAtomRef,
                     parameters: AmideMetricAuthority, *, backbone_types: tuple[str, str, str] | None) -> StagedResult:
    """No search/holdout/phase/optimizer interface, and no H-dependent reselection."""
    binding = parent_digest(parent)
    try:
        view = compile_staged_view(original, parent, authority, incoming, parameters, backbone_types)
    except ValueError as exc:
        # Refusal contains no accepted output; full supplied parent remains caller evidence.
        status: Literal["parent-refused", "H-refused"] = (
            "parent-refused" if str(exc).startswith("staged-parent") or parent.status != "selected" else "H-refused")
        return StagedResult(original.digest, binding, status, str(exc))
    proposal = propose(view.request, view.plan)
    assert parent.selected_index is not None
    checkpoint = parent.attempts[parent.selected_index].rounded
    assert checkpoint is not None
    raw = (StagedEnvelope(original.digest, binding, view.request.digest, digest_record(proposal.raw),
                          "raw-H-after-rounded-heavy", checkpoint.additions, proposal.raw.additions)
           if proposal.raw is not None else None)
    rounded = (StagedEnvelope(original.digest, binding, view.request.digest, digest_record(proposal.rounded),
                              "rounded-H-after-rounded-heavy", checkpoint.additions, proposal.rounded.additions)
               if proposal.rounded is not None else None)
    audit = measure_staged(original, parent, authority, incoming, parameters, backbone_types, view, proposal, raw, rounded)
    return StagedResult(original.digest, binding, "staged-geometry-passed" if audit.passed else "H-rejected",
                        "conditional-declared-checks-only" if audit.passed else "retained-H-stage-failure-no-heavy-reselection",
                        view, proposal, raw, rounded, audit)
