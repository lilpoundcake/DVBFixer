"""One complete-H proposal after an immutable typed projection checkpoint."""
from __future__ import annotations

from .amide_authority import AmideMetricAuthority
from .construction import propose
from .contracts import digest_record
from .projected_hydrogen_audit import compile_projected_h_view, measure_projected_hydrogens
from .projected_hydrogen_contracts import ProjectedHeavyParent, ProjectedHydrogenResult
from .staged_contracts import StagedEnvelope


def complete_projected_hydrogens(parent: ProjectedHeavyParent, parameters: AmideMetricAuthority, *,
                                  backbone_types: tuple[str, str, str] | None) -> ProjectedHydrogenResult:
    try:
        view = compile_projected_h_view(parent, parameters, backbone_types)
    except ValueError as exc:
        return ProjectedHydrogenResult(parent.original.digest, parent.request.digest, parent.digest,
            "parent-refused" if str(exc).startswith("projected-parent-audit:") else "H-refused", str(exc))
    proposal = propose(view.request, view.plan)
    checkpoint = parent.projection.rounded
    assert checkpoint is not None
    raw = (StagedEnvelope(parent.original.digest, parent.digest, view.request.digest, digest_record(proposal.raw),
                          "raw-H-after-rounded-heavy", checkpoint.additions, proposal.raw.additions)
           if proposal.raw is not None else None)
    rounded = (StagedEnvelope(parent.original.digest, parent.digest, view.request.digest, digest_record(proposal.rounded),
                              "rounded-H-after-rounded-heavy", checkpoint.additions, proposal.rounded.additions)
               if proposal.rounded is not None else None)
    audit = measure_projected_hydrogens(parent, parameters, backbone_types, view, proposal, raw, rounded)
    return ProjectedHydrogenResult(parent.original.digest, parent.request.digest, parent.digest,
        "projected-H-geometry-passed" if audit.passed else "H-rejected",
        "projected-heavy-H-conditional-declared-checks-only" if audit.passed else "retained-projected-H-failure-no-heavy-change",
        view, proposal, raw, rounded, audit)
