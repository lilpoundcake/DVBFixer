"""Constructor-free combined lineage/geometry/CIP and explicit pair accounting."""
from __future__ import annotations

from dataclasses import replace

from dvbfixer.domain.structure_identity import ExactAtomRef

from .amide_authority import AmideMetricAuthority
from .candidate_audit import audit_candidate
from .conformer_search import ConformerSearchResult, audit_search_result
from .construction import proposal_admission
from .contracts import LinkKind, ReconstructionRequest, digest_record
from .declared_peptide import DeclaredPeptideAuthority, compile_declared_hydrogens
from .heavy_conditioning import condition_on_heavy
from .hydrogen_completion import CompleteHydrogenView, audit_hydrogen_cip
from .observed_peptide import compile_peptide_request_plan
from .proposals import Candidate, ProposalResult
from .staged_contracts import (
    FIXED_POLICY,
    HYDROGEN_POLICY,
    PairCoverage,
    StagedAudit,
    StagedEnvelope,
    StagedResult,
    parent_digest,
)


def compile_staged_view(original: ReconstructionRequest, parent: ConformerSearchResult,
                        authority: DeclaredPeptideAuthority, incoming: ExactAtomRef,
                        parameters: AmideMetricAuthority, backbone_types: tuple[str, str, str] | None
                        ) -> CompleteHydrogenView:
    """Recompile declarations ONLY; never construct/rotate coordinates."""
    try:
        errors = audit_search_result(original, authority, parent)
        if errors:
            raise ValueError(",".join(errors))
        if (parent.plan is None or parent.status != "selected" or type(parent.selected_index) is not int
                or not 0 <= parent.selected_index < len(parent.attempts)
                or parent.plan.base_plan != compile_peptide_request_plan(original, FIXED_POLICY)):
            raise ValueError("needs-authority-bound-selected-parent")
        baseline, plan = parent.baseline, parent.plan.base_plan
        if (baseline.plan_digest != plan.digest or baseline.admission != proposal_admission(original, plan)
                or baseline.raw is None or baseline.rounded is None
                or baseline.raw_audit != audit_candidate(original, plan, baseline.raw)
                or baseline.rounded_audit != audit_candidate(original, plan, baseline.rounded)):
            raise ValueError("baseline-metadata-or-audit-binding")
        selected = parent.attempts[parent.selected_index]
        if not selected.selectable or selected.rounded is None:
            raise ValueError("needs-passing-rounded-heavy-checkpoint")
        base = condition_on_heavy(original, selected.rounded)
    except ValueError as exc:
        raise ValueError("staged-parent-audit:" + str(exc)) from exc
    return compile_declared_hydrogens(base, authority, incoming, parameters, backbone_types=backbone_types,
                                      fixed_policy=FIXED_POLICY, generated_policy=HYDROGEN_POLICY)


def pair_coverage(original: ReconstructionRequest, view: CompleteHydrogenView) -> PairCoverage:
    """Enumerate the union once, partitioned by actual source/heavy/H ownership."""
    source = {o.identity for o in original.observations}
    heavy, hydrogens = original.ownership.generated, view.request.ownership.generated
    if source & (heavy | hydrogens) or heavy & hydrogens:
        raise ValueError("staged-pair-ownership-overlap")
    if {o.identity for o in view.request.observations} != source | heavy:
        raise ValueError("staged-pair-context-inventory")
    adjacent: dict[ExactAtomRef, set[ExactAtomRef]] = {}
    for component in original.components:
        if component.resolved is not None:
            for bond in component.resolved.graph.bonds:
                a, b = ExactAtomRef(component.identity, bond.atom1), ExactAtomRef(component.identity, bond.atom2)
                adjacent.setdefault(a, set()).add(b)
                adjacent.setdefault(b, set()).add(a)
    for connection in original.connections:
        if connection.kind == LinkKind.COVALENT:
            a, b = connection.link.component_atom, connection.link.other_atom
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    pairs: dict[tuple[ExactAtomRef, ExactAtomRef], str] = {}
    for a in sorted(heavy | hydrogens):
        near = adjacent.get(a, set())
        excluded = {a} | near | {c for b in near for c in adjacent.get(b, set())}
        for b in sorted(source | heavy | hydrogens):
            if b in excluded:
                continue
            pair = min(a, b), max(a, b)
            if a in hydrogens or b in hydrogens:
                kind = "HH" if a in hydrogens and b in hydrogens else "Hh" if a in heavy or b in heavy else "Hs"
            else:
                kind = "hh" if a in heavy and b in heavy else "hs"
            pairs[pair] = kind
            if len(pairs) > 1_000_000:
                raise ValueError("staged-pair-accounting-budget")
    kinds = tuple(pairs.values())
    return PairCoverage(kinds.count("hs"), kinds.count("hh"), kinds.count("Hs"), kinds.count("Hh"), kinds.count("HH"),
                        digest_record(tuple((a, b, kind) for (a, b), kind in sorted(pairs.items()))))


def _round(candidate: Candidate) -> Candidate:
    additions = []
    for atom in candidate.additions:
        fields = tuple(f"{x:8.3f}" for x in atom.position)
        if any(len(field) != 8 for field in fields):
            raise ValueError("staged-H-PDB-coordinate-range")
        additions.append(replace(atom, position=(float(fields[0]), float(fields[1]), float(fields[2]))))
    return replace(candidate, additions=tuple(additions))


def measure_staged(original: ReconstructionRequest, parent: ConformerSearchResult,
                    authority: DeclaredPeptideAuthority, incoming: ExactAtomRef,
                    parameters: AmideMetricAuthority, backbone_types: tuple[str, str, str] | None,
                    view: CompleteHydrogenView, proposal: ProposalResult,
                    raw: StagedEnvelope | None, rounded: StagedEnvelope | None) -> StagedAudit:
    """Rebind every declarative input and measure saved candidates, never propose."""
    try:
        expected = compile_staged_view(original, parent, authority, incoming, parameters, backbone_types)
    except ValueError as exc:
        return StagedAudit(False, ("staged-recompile-refused:" + str(exc),))
    if digest_record(view) != digest_record(expected):
        return StagedAudit(False, ("staged-H-declaration-binding",))
    errors = set()
    if proposal.plan_digest != view.plan.digest:
        errors.add("staged-H-result-plan-binding")
    if proposal.admission != proposal_admission(view.request, view.plan):
        errors.add("staged-H-admission-metadata")
    raw_geometry = audit_candidate(view.request, view.plan, proposal.raw) if proposal.raw is not None else None
    rounded_geometry = audit_candidate(view.request, view.plan, proposal.rounded) if proposal.rounded is not None else None
    raw_cip = audit_hydrogen_cip(view, authority, proposal.raw) if proposal.raw is not None else None
    rounded_cip = audit_hydrogen_cip(view, authority, proposal.rounded) if proposal.rounded is not None else None
    if (raw_geometry != proposal.raw_audit or rounded_geometry != proposal.rounded_audit):
        errors.add("staged-H-saved-geometry-audit-mismatch")
    if not raw_geometry or not raw_geometry.passed or not rounded_geometry or not rounded_geometry.passed:
        errors.add("staged-H-raw-or-rounded-geometry")
    if not raw_cip or not raw_cip.passed or not rounded_cip or not rounded_cip.passed:
        errors.add("staged-H-raw-or-rounded-CIP")
    if proposal.status != "declared-geometry-passed":
        errors.add("staged-H-proposal-not-passed")
    try:
        if proposal.raw is None or proposal.rounded is None or digest_record(_round(proposal.raw)) != digest_record(proposal.rounded):
            errors.add("staged-H-rounding-transport")
    except ValueError:
        errors.add("staged-H-rounding-transport")
    assert parent.selected_index is not None
    checkpoint = parent.attempts[parent.selected_index].rounded
    assert checkpoint is not None
    for candidate, envelope, phase in ((proposal.raw, raw, "raw-H-after-rounded-heavy"),
                                        (proposal.rounded, rounded, "rounded-H-after-rounded-heavy")):
        if candidate is None:
            if envelope is not None:
                errors.add("staged-envelope-without-H-candidate")
            continue
        if envelope is None:
            errors.add("staged-missing-envelope")
            continue
        if (envelope.original_request_digest != original.digest or envelope.parent_digest != parent_digest(parent)
                or envelope.conditioning_digest != view.request.digest
                or envelope.hydrogen_candidate_digest != digest_record(candidate) or envelope.phase != phase
                or digest_record(envelope.generated_heavy) != digest_record(checkpoint.additions)
                or digest_record(envelope.generated_hydrogens) != digest_record(candidate.additions)):
            errors.add("staged-envelope-lineage-or-role-transport")
    try:
        if (proposal.raw is None or proposal.rounded is None
                or {a.identity for a in proposal.raw.additions} != view.request.ownership.generated
                or {a.identity for a in proposal.rounded.additions} != view.request.ownership.generated):
            raise ValueError("staged-pair-H-inventory-incomplete")
        coverage = pair_coverage(original, view)
    except ValueError as exc:
        errors.add(str(exc))
        coverage = None
    return StagedAudit(not errors, tuple(sorted(errors)), raw_geometry, rounded_geometry, raw_cip, rounded_cip, coverage)


def audit_staged_result(original: ReconstructionRequest, parent: ConformerSearchResult,
                         authority: DeclaredPeptideAuthority, incoming: ExactAtomRef,
                         parameters: AmideMetricAuthority, backbone_types: tuple[str, str, str] | None,
                         result: StagedResult) -> tuple[str, ...]:
    """Readback integrity, distinct from scientific acceptance of the retained result."""
    errors = set()
    if result.original_request_digest != original.digest or result.parent_digest != parent_digest(parent):
        errors.add("staged-result-parent-or-source-binding")
    try:
        expected = compile_staged_view(original, parent, authority, incoming, parameters, backbone_types)
    except ValueError as exc:
        status = "parent-refused" if str(exc).startswith("staged-parent") or parent.status != "selected" else "H-refused"
        if (result.status != status or result.reason != str(exc)
                or any(x is not None for x in (result.view, result.hydrogen_result, result.raw, result.rounded, result.audit))):
            errors.add("staged-refusal-record-mismatch")
        return tuple(sorted(errors))
    if result.view is None or result.hydrogen_result is None:
        return tuple(sorted(errors | {"staged-missing-H-evidence"}))
    if digest_record(result.view) != digest_record(expected):
        errors.add("staged-H-declaration-binding")
    audit = measure_staged(original, parent, authority, incoming, parameters, backbone_types,
                            result.view, result.hydrogen_result, result.raw, result.rounded)
    if audit != result.audit:
        errors.add("staged-saved-combined-audit-mismatch")
    status = "staged-geometry-passed" if audit.passed else "H-rejected"
    if result.status != status:
        errors.add("staged-result-status-mismatch")
    reason = "conditional-declared-checks-only" if audit.passed else "retained-H-stage-failure-no-heavy-reselection"
    if result.reason != reason:
        errors.add("staged-result-reason-mismatch")
    return tuple(sorted(errors))
