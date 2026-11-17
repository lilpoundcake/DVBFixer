"""Post-assembly bounded-anchor repair; donors and frozen requests never move.

The assembled projection is the raw restart frame. Its generated-only local
attempt supplies residual evidence, not a new starting frame or a new source.
Existing staged anchor acceptance and movement limits remain authoritative.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import (
    _findings,
    _structure,
    audit_staged_anchor_repair,
)
from dvbfixer.model.diffusion.anchor_policy import build_anchor_restraints, check_anchor_movements
from dvbfixer.model.diffusion.contract import DiffusionContractError
from dvbfixer.model.diffusion.gap_portfolio import _coordinates, digest, plan_gap_repairs
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, problematic_anchors
from dvbfixer.model.diffusion.runner_refinement import (
    _protein_only_refinement_text,
    rewrite_generated_coordinates,
)
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

REVISION = "post-portfolio-bounded-anchors-v1"


def anchor_repair_scope(request, source_text, selected_text, stage1_text):
    """Recompute eligible direct anchors and issue-local generated ownership."""
    if digest(source_text) != request.normalized_pdb.sha256:
        raise DiffusionContractError("post-portfolio source binding differs")
    assess = assessment_adapter(request, source_text)
    raw, first = assess(selected_text), assess(stage1_text)
    if raw.admission_failures or first.admission_failures:
        raise DiffusionContractError("post-portfolio inputs have admission failures")
    if raw.geometry_passed or first.geometry_passed:
        raise DiffusionContractError("post-portfolio anchor repair is unnecessary")
    first_scope = plan_gap_repairs(request, raw)
    first_residues = {r for i in first_scope.gap_indices for r in request.gaps[i].generated_residues}
    first_atoms = {a for a in request.generated_atoms
        if any((a.chain, a.residue_number, a.insertion_code) ==
               (r.chain, r.residue_number, r.insertion_code) for r in first_residues)}
    first_xyz, _ = _coordinates(stage1_text)
    if not first_atoms <= first_xyz.keys() or stage1_text != rewrite_generated_coordinates(
            selected_text, {a: first_xyz[a] for a in first_atoms}):
        raise DiffusionContractError("post-portfolio stage-1 changed observations or clean gaps")
    baseline = {(issue.category, issue.residues) for issue in _findings(_structure(source_text, "source.pdb"), request)}
    eligible = RepairAssessment(tuple(issue for issue in first.issues
        if (issue.category, issue.residues) not in baseline))
    anchors = problematic_anchors(request, eligible)
    if not anchors or set(problematic_anchors(request, first)) - set(anchors):
        raise DiffusionContractError("post-portfolio has no exclusively reconstruction-justified anchors")
    indices = tuple(sorted(set(plan_gap_repairs(request, raw).gap_indices) |
                           set(plan_gap_repairs(request, first).gap_indices)))
    residues = tuple(dict.fromkeys(r for i in indices for r in request.gaps[i].generated_residues))
    identities = {(r.chain, r.residue_number, r.insertion_code) for r in residues}
    atoms = tuple(a for a in request.generated_atoms
        if (a.chain, a.residue_number, a.insertion_code) in identities)
    if not atoms:
        raise DiffusionContractError("post-portfolio anchor repair has no affected generated gaps")
    return indices, residues, atoms, anchors


def audit_portfolio_anchor_repair(request, source_text, selected_text, stage1_text, candidate_text):
    """Independent custody/mobility audit; never trusts runner masks or receipts.

    Whole-gap donor lineage and stage-1 local custody are checked by the outer
    portfolio auditor before this function. Here every nonmobile record, including
    headers, links, clean generated gaps and anchor metadata, must remain exact.
    """
    indices, residues, atoms, anchors = anchor_repair_scope(request, source_text, selected_text, stage1_text)
    restraints = build_anchor_restraints(source_text, anchors, observed_atoms=request.fixed_atoms, bfactor_source="auto")
    xyz, rows = _coordinates(candidate_text)
    _, original = _coordinates(selected_text)
    allowed = set(atoms) | {r.atom for r in restraints}
    if set(rows) != set(original) or not allowed <= xyz.keys():
        raise DiffusionContractError("post-portfolio candidate changed atom inventory")
    if candidate_text != rewrite_generated_coordinates(selected_text, {a: xyz[a] for a in allowed}):
        raise DiffusionContractError("post-portfolio candidate changed nonmobile records or metadata")
    audit = audit_staged_anchor_repair(request, source_text=source_text, raw_text=selected_text,
        stage1_text=stage1_text, candidate_text=candidate_text,
        stage1_input_sha256=digest(selected_text), stage2_input_sha256=digest(selected_text),
        limits=RefinementPolicy().movement_limits)
    if audit.authorized_anchors != anchors:
        raise DiffusionContractError("post-portfolio independent anchor authorization differs")
    structure = _structure(candidate_text, "post-portfolio.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    return audit


@dataclass(frozen=True)
class PortfolioAnchorRepair:
    coordinate_text: str
    report: dict


def repair_portfolio_anchors(request, source_text, selected_text, stage1_text, *,
                             refinement_seed, max_iterations=250, restart_count=4, platform_name="CPU",
                             retain_trial=None):
    """One explicit escalation, restarted from the assembled donor projection.

    Returns rejected geometry as diagnostics too; publication always re-audits.
    Numerical/configuration failures propagate, leaving script-owned evidence.
    """
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region

    if (type(refinement_seed) is not int or refinement_seed < 0 or type(max_iterations) is not int
            or not 1 <= max_iterations <= 2000 or type(restart_count) is not int or not 0 <= restart_count <= 8):
        raise DiffusionContractError("invalid post-portfolio refinement budget")
    indices, residues, atoms, anchors = anchor_repair_scope(request, source_text, selected_text, stage1_text)
    assess = assessment_adapter(request, source_text)
    restraints = build_anchor_restraints(source_text, anchors, observed_atoms=request.fixed_atoms, bfactor_source="auto")
    allowed = frozenset(r.atom for r in restraints)
    limits = RefinementPolicy().movement_limits
    trials = []

    def score(xyz):
        projected = rewrite_generated_coordinates(selected_text, xyz)
        assessment = assess(projected, allowed)
        movement = check_anchor_movements(restraints, _coordinates(projected)[0], limits)
        value = (len(assessment.admission_failures), sum(not m.within_limits for m in movement),
                 sum(i.actionable for i in assessment.issues))
        trials.append({"score": value, "assessment": asdict(assessment), "movements": [asdict(m) for m in movement],
                       "coordinate_sha256": digest(projected)})
        if retain_trial is not None:
            retain_trial(len(trials) - 1, projected, trials[-1])
        return value

    result = refine_generated_region(
        _protein_only_refinement_text(selected_text) if request.heterogen_contexts else selected_text,
        generated_residues=residues, generated_atoms=atoms, observed_anchor_restraints=restraints,
        max_iterations=max_iterations, restart_count=restart_count, random_seed=refinement_seed,
        perturbation_angstrom=0.25, platform_name=platform_name, geometry_preconditioning=True,
        candidate_score=score, anchor_movement_limits=limits)
    text = rewrite_generated_coordinates(selected_text, result.coordinates_angstrom)
    audit = audit_portfolio_anchor_repair(request, source_text, selected_text, stage1_text, text)
    report = {"operation": REVISION, "input_sha256": digest(selected_text), "stage1_sha256": digest(stage1_text),
        "output_sha256": digest(text), "request_sha256": digest(request.to_json()), "source_sha256": digest(source_text),
        "gap_indices": indices, "mobile_generated_atoms": [asdict(a) for a in atoms],
        "authorized_anchors": [asdict(r) for r in anchors], "restraints": [asdict(r) for r in restraints],
        "movement_limits": asdict(limits), "refinement_seed": refinement_seed, "max_iterations": max_iterations,
        "restart_count": restart_count, "platform": result.platform, "audit": asdict(audit),
        "trials": trials, "search_diagnostics": result.search_diagnostics,
        "publication_qualified": False, "new_sampling": False}
    return PortfolioAnchorRepair(text, report)
