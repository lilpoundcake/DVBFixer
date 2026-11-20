"""Constructor/solver-free projected-parent binding and shared combined-H audit."""
from __future__ import annotations

from .amide_authority import AmideMetricAuthority
from .contracts import digest_record
from .declared_peptide import compile_declared_hydrogens
from .hydrogen_completion import CompleteHydrogenView
from .observed_peptide import compile_peptide_request_plan
from .planar_port import planar_port_constraints
from .projected_conditioning import condition_on_projected
from .projected_hydrogen_contracts import ProjectedHeavyParent, ProjectedHydrogenResult
from .projection_audit import audit_projection_result
from .projection_contracts import append_projection_constraints
from .proposals import ProposalResult
from .staged_audit import measure_bound_hydrogens
from .staged_contracts import FIXED_POLICY, HYDROGEN_POLICY, StagedAudit, StagedEnvelope


def compile_projected_h_view(parent: ProjectedHeavyParent, parameters: AmideMetricAuthority,
                              backbone_types: tuple[str, str, str] | None) -> CompleteHydrogenView:
    try:
        p = parent.projection
        if p.plan is None or p.rounded is None or p.audit is None or not p.audit.passed or p.status != "projected-geometry-passed":
            raise ValueError("passing-complete-rounded-projection-required")
        if p.plan.original.digest != parent.original.digest:
            raise ValueError("projection-original-source-request-binding")
        extra = planar_port_constraints(parent.original, parent.peptide_authority, parent.port_declaration, parent.port_parameters)
        if append_projection_constraints(parent.original, extra).digest != parent.request.digest:
            raise ValueError("actual-port-declaration-request-binding")
        if p.plan.baseline_plan != compile_peptide_request_plan(parent.original, FIXED_POLICY):
            raise ValueError("original-peptide-baseline-plan-binding")
        findings = audit_projection_result(parent.original, parent.request, p.plan.baseline_plan,
                                           p.plan.baseline, p.plan.anchors, p)
        if findings:
            raise ValueError(",".join(findings))
        base = condition_on_projected(parent.request, p)
    except ValueError as exc:
        raise ValueError("projected-parent-audit:" + str(exc)) from exc
    return compile_declared_hydrogens(base, parent.peptide_authority, parent.port_declaration.incoming, parameters,
                                      backbone_types=backbone_types, fixed_policy=FIXED_POLICY, generated_policy=HYDROGEN_POLICY)


def measure_projected_hydrogens(parent: ProjectedHeavyParent, parameters: AmideMetricAuthority,
                                 backbone_types: tuple[str, str, str] | None, view: CompleteHydrogenView,
                                 proposal: ProposalResult, raw: StagedEnvelope | None,
                                 rounded: StagedEnvelope | None) -> StagedAudit:
    try:
        expected = compile_projected_h_view(parent, parameters, backbone_types)
    except ValueError as exc:
        return StagedAudit(False, ("projected-H-recompile-refused:" + str(exc),))
    if digest_record(view) != digest_record(expected):
        return StagedAudit(False, ("projected-H-declaration-binding",))
    assert parent.projection.rounded is not None
    return measure_bound_hydrogens(parent.original, parent.digest, parent.projection.rounded,
                                    parent.peptide_authority, view, proposal, raw, rounded)


def audit_projected_hydrogen_result(parent: ProjectedHeavyParent, parameters: AmideMetricAuthority,
                                     backbone_types: tuple[str, str, str] | None,
                                     result: ProjectedHydrogenResult) -> tuple[str, ...]:
    errors = set()
    if (result.original_request_digest != parent.original.digest or result.port_request_digest != parent.request.digest
            or result.parent_digest != parent.digest):
        errors.add("projected-H-original-port-parent-binding")
    try:
        compile_projected_h_view(parent, parameters, backbone_types)
    except ValueError as exc:
        status = "parent-refused" if str(exc).startswith("projected-parent-audit:") else "H-refused"
        if (result.status != status or result.reason != str(exc)
                or any(x is not None for x in (result.view, result.hydrogen_result, result.raw, result.rounded, result.audit))):
            errors.add("projected-H-refusal-record-mismatch")
        return tuple(sorted(errors))
    if result.view is None or result.hydrogen_result is None:
        return tuple(sorted(errors | {"projected-H-missing-evidence"}))
    audit = measure_projected_hydrogens(parent, parameters, backbone_types, result.view,
                                        result.hydrogen_result, result.raw, result.rounded)
    if digest_record(audit) != digest_record(result.audit):
        errors.add("projected-H-saved-audit-mismatch")
    status = "projected-H-geometry-passed" if audit.passed else "H-rejected"
    reason = "projected-heavy-H-conditional-declared-checks-only" if audit.passed else "retained-projected-H-failure-no-heavy-change"
    if result.status != status or result.reason != reason:
        errors.add("projected-H-status-or-reason-mismatch")
    return tuple(sorted(errors))
