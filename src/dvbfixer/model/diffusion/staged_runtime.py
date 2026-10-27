"""Scientific adapters for geometry-gated runner execution and independent audit."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

from dvbfixer.model.diffusion.anchor_acceptance import (
    _assessment,
    _observed_auxiliary,
    _observed_lines,
    _structure,
    audit_staged_anchor_repair,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionContractError,
    DiffusionRequest,
    ValidationSummary,
)
from dvbfixer.model.diffusion.refinement_evidence import StagedRefinementEvidence
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy, read_refinement_policy
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment

DRIFT_GATES = frozenset({"fixed-heavy-atom-rmsd", "fixed-heavy-atom-max-displacement"})


def assessment_adapter(request: DiffusionRequest, source_text: str):
    """Assess expanded geometry and exact observed records under an explicit mask."""
    source = _structure(source_text, "source.pdb")
    source_lines = _observed_lines(source_text)
    observed = source.atoms.keys() - set(request.generated_atoms)

    def assess(text, allowed=frozenset()):
        candidate = _structure(text, "candidate.pdb")
        assessment, _ = _assessment(request, source, candidate)
        failures = [gate for gate in assessment.admission_failures if not allowed or gate not in DRIFT_GATES]
        lines = _observed_lines(text)
        for atom in observed:
            original, current = source_lines.get(atom), lines.get(atom)
            if current is None or original is None:
                failures.append("staged-observed-record-missing")
            elif atom not in allowed and original != current:
                failures.append("staged-unselected-observed-record-changed")
            elif atom in allowed and (original[:30] != current[:30] or original[54:] != current[54:]):
                failures.append("staged-observed-metadata-changed")
        if _observed_auxiliary(source_text, observed) != _observed_auxiliary(text, observed):
            failures.append("staged-observed-auxiliary-changed")
        return RepairAssessment(assessment.issues, tuple(dict.fromkeys(failures)))

    return assess


def validate_staged_evidence(
    request: DiffusionRequest,
    evidence: StagedRefinementEvidence,
    *,
    source_text: str,
    selected: ArtifactReference,
    load_text: Callable[[ArtifactReference], str],
    seed: int | None = None,
) -> ValidationSummary:
    """Recompute decisions from verified artifacts, never runner masks or scores."""
    evidence.validate_request(request, seed=seed)
    policy = read_refinement_policy(request.backend_options)
    assert isinstance(policy, RefinementPolicy)
    if selected != evidence.selected:
        raise DiffusionContractError("selected candidate does not match staged evidence")
    raw_text = load_text(evidence.raw)
    candidate_text = load_text(evidence.selected)
    first_text = load_text(evidence.stage1) if evidence.stage1 is not None else None
    assess = assessment_adapter(request, source_text)
    raw = assess(raw_text)
    if evidence.selected_stage == "generated-and-anchors":
        audit = audit_staged_anchor_repair(
            request, source_text=source_text, raw_text=raw_text, stage1_text=first_text,
            candidate_text=candidate_text, stage1_input_sha256=evidence.stage1_input_sha256,
            stage2_input_sha256=evidence.stage2_input_sha256, limits=policy.movement_limits,
            stage1_failure=evidence.stage1_failure,
        )
        # This versioned public adapter, not the research auditor alone, grants
        # the request-bound exception. Historical strict results stay in metrics.
        return replace(audit.summary, warnings=tuple(
            warning for warning in audit.summary.warnings
            if warning != "Research anchor audit only; not qualified for production publication."
        ))
    final = assess(candidate_text)
    failures = [*raw.admission_failures, *final.admission_failures]
    if not final.geometry_passed:
        failures.append("staged-expanded-geometry")
    if evidence.selected_stage == "generated-only" and raw.geometry_passed:
        failures.append("staged-unnecessary-raw-clean")
    return ValidationSummary(not failures, tuple(dict.fromkeys(failures)))
