"""Testable raw-based two-stage repair controller, separate from publication.

Independent assessors supply geometry and admission findings. A geometry repair
is not a publication decision. No existing runner/validator policy is changed.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

from dvbfixer.model.diffusion.anchor_policy import (
    AnchorPositionRestraint,
    NumericalRefinementError,
    build_anchor_restraints,
)
from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity


@dataclass(frozen=True)
class RepairIssue:
    category: str
    residues: tuple[ResidueIdentity, ...]
    message: str
    actionable: bool = True


@dataclass(frozen=True)
class RepairAssessment:
    issues: tuple[RepairIssue, ...] = ()
    admission_failures: tuple[str, ...] = ()

    @property
    def geometry_passed(self) -> bool:
        return not any(issue.actionable for issue in self.issues)


@dataclass(frozen=True)
class RepairAttempt:
    stage: Literal["raw", "generated-only", "generated-and-anchors"]
    input_sha256: str
    output_text: str | None
    assessment: RepairAssessment | None
    restraints: tuple[AnchorPositionRestraint, ...] = ()
    error: str | None = None
    failure_kind: Literal["numerical", "unclassified"] | None = None
    numerical_failure: NumericalFailureEvidence | None = None


@dataclass(frozen=True)
class NumericalFailureEvidence:
    """Bounded stage-1 failure provenance; eligibility is still recomputed from raw."""

    input_sha256: str
    error_type: str
    message: str
    revision: int = 1

    def __post_init__(self) -> None:
        if type(self.revision) is not int or self.revision != 1:
            raise ValueError("unsupported numerical failure evidence revision")
        if not isinstance(self.input_sha256, str) or len(self.input_sha256) != 64 or any(c not in "0123456789abcdef" for c in self.input_sha256):
            raise ValueError("numerical failure evidence needs a lowercase SHA-256")
        if not isinstance(self.error_type, str) or self.error_type not in {"NumericalRefinementError", "BoundaryNumericalError"}:
            raise ValueError("unclassified failures cannot authorize numerical recovery")
        if not isinstance(self.message, str) or not self.message or len(self.message) > 2048:
            raise ValueError("numerical failure message must be nonempty and bounded")


@dataclass(frozen=True)
class StagedRepairResult:
    selected_stage: str | None
    coordinate_text: str | None
    attempts: tuple[RepairAttempt, ...]
    selected_anchors: tuple[ResidueIdentity, ...] = ()
    # This controller intentionally never grants an observed-motion exception to
    # the production validator. Qualification belongs to the versioned contract.
    requires_anchor_acceptance_contract: bool = False


_ANCHOR_CATEGORIES = frozenset({"junction-peptide-connectivity", "bond_length", "bond_angle",
                               "non_planar_amide", "ramachandran_outlier", "clash"})


def problematic_anchors(request: DiffusionRequest, assessment: RepairAssessment) -> tuple[ResidueIdentity, ...]:
    """Only direct flanks actually named by residual actionable findings."""
    anchors = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    generated = {r for gap in request.gaps for r in gap.generated_residues}
    chosen = set()
    for issue in assessment.issues:
        if not issue.actionable or issue.category not in _ANCHOR_CATEGORIES:
            continue
        involved = set(issue.residues)
        if issue.category == "clash" and not involved & generated:
            continue
        chosen.update(involved & anchors)
    return tuple(sorted(chosen))


def run_refinement_stages(
    raw_text: str,
    source_text: str,
    request: DiffusionRequest,
    *,
    assess: Callable[[str], RepairAssessment],
    minimize: Callable[[str, tuple[AnchorPositionRestraint, ...]], str],
    relax_anchors: bool = True,
    bfactor_source: Literal["auto", "experimental", "disabled"] = "auto",
    notify: Callable[[str], None] | None = None,
) -> StagedRepairResult:
    """Default-enabled escalation; both minimizations receive the SAME raw text.

    assess must identify nongeometric admission failures separately, cover outer
    junctions and all-L, and retain all affected residue identities. The caller
    must run full independent acceptance before publishing any selected output.
    Stage 2 requires a separate anchor acceptance contract even on geometry pass.
    RuntimeError/ValueError are retained as attempt failures, but only an explicit
    NumericalRefinementError can escalate after a failure. No message guessing.
    Other programming errors, interrupts and cancellation propagate.
    """
    raw_digest = hashlib.sha256(raw_text.encode("ascii")).hexdigest()
    raw_assessment = assess(raw_text)
    attempts = [RepairAttempt("raw", raw_digest, raw_text, raw_assessment)]
    if raw_assessment.admission_failures:
        return StagedRepairResult(None, None, tuple(attempts))
    if raw_assessment.geometry_passed:
        return StagedRepairResult("raw", raw_text, tuple(attempts))

    def attempt(stage, restraints):
        try:
            text = minimize(raw_text, restraints)
            assessment = assess(text)
            result = RepairAttempt(stage, raw_digest, text, assessment, restraints)
        except (RuntimeError, ValueError) as exc:
            kind = "numerical" if isinstance(exc, NumericalRefinementError) else "unclassified"
            evidence = None
            if kind == "numerical":
                error_type = "BoundaryNumericalError" if type(exc).__name__ == "BoundaryNumericalError" else "NumericalRefinementError"
                evidence = NumericalFailureEvidence(raw_digest, error_type, str(exc)[:2048] or "Numerical minimization failure")
            result = RepairAttempt(stage, raw_digest, None, None, restraints,
                                   f"{type(exc).__name__}: {exc}"[:2048], kind, evidence)
        attempts.append(result)
        return result

    first = attempt("generated-only", ())
    if first.assessment is not None:
        if first.assessment.admission_failures:
            return StagedRepairResult(None, None, tuple(attempts))
        if first.assessment.geometry_passed:
            return StagedRepairResult(first.stage, first.output_text, tuple(attempts))
    elif first.failure_kind != "numerical":
        return StagedRepairResult(None, None, tuple(attempts))
    residual = first.assessment if first.assessment is not None else raw_assessment
    selected = problematic_anchors(request, residual)
    if not relax_anchors or not selected:
        return StagedRepairResult(None, None, tuple(attempts))
    restraints = build_anchor_restraints(source_text, selected, observed_atoms=request.fixed_atoms,
                                         bfactor_source=bfactor_source)
    if notify is not None:
        labels = ", ".join(f"{r.chain}:{r.residue_number}{r.insertion_code}" for r in selected)
        notify(f"WARNING: Generated-only repair did not resolve junction geometry. "
               f"Retrying from raw with restrained observed anchors: {labels}.")
    second = attempt("generated-and-anchors", restraints)
    if second.assessment is not None and not second.assessment.admission_failures and second.assessment.geometry_passed:
        return StagedRepairResult(second.stage, second.output_text, tuple(attempts), selected, True)
    return StagedRepairResult(None, None, tuple(attempts), selected, True)
