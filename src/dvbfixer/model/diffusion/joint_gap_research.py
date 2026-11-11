"""Private, source-only union transport for one jointly sampled group of gaps.

Not a public profile/admission change. In particular no independent runner is
relabelled joint and no other generated prediction becomes conditioning.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    BackendOption,
    ContextCoordinate,
    DiffusionContractError,
    DiffusionRequest,
    GapKind,
    ProteinContextPlan,
    ProteinContextPolicy,
    ProteinContextSpan,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, ATOM37_ORDER
from dvbfixer.model.diffusion.protein_context import (
    EngineContextView,
    _split_native_segments,
    engine_context_view,
    independent_gap_requests,
    source_context_position_features,
    source_coordinates,
    validate_cc94_request,
    validate_context_source,
)

REVISION = "private-cc94-joint-gap-source-union-v1"
STRUCTURE_REVISION = "private-cc94-joint-structure-source-union-v2"
COMPLETE_REVISION = "private-cc94-complete-observed-context-v3"
# Explicitly authorized 1HE8 complete-context research shape; not public support.
COMPLETE_EXECUTION_TOKEN_LIMIT = 1131
MOTIF_COMPONENT_REVISION = "private-cc94-neighbor-motif-components-v4"
# Keep the explicitly authorized private ceiling; >512 remains experimental.
MOTIF_COMPONENT_TOKEN_LIMIT = COMPLETE_EXECUTION_TOKEN_LIMIT


@dataclass(frozen=True)
class JointGapResearchView:
    request: DiffusionRequest
    engine_view: EngineContextView
    gap_indices: tuple[int, ...]
    parent_request_sha256: str
    positions: tuple[int, ...]
    chains: tuple[int, ...]


def plan_complete_context(parent, source: bytes, gap_indices: tuple[int, ...], *, seed: int):
    """Plan ALL observed protein atoms plus the selected unknown gaps.

    No observed body or partner chain can become output-only. Unselected unknown
    gaps remain absent, never predictions used as conditioning. Planning can
    describe an over-budget input; complete_context_admission separately and
    mandatorily controls execution. Old v1/v2 plans remain historical, unchanged.
    """
    if (not 1 <= len(gap_indices) <= 16 or tuple(sorted(set(gap_indices))) != gap_indices
            or any(type(i) is not int or not 0 <= i < len(parent.gaps) for i in gap_indices)
            or type(seed) is not int or not 0 <= seed < 2**32):
        raise DiffusionContractError("complete context needs bounded distinct gaps and one explicit seed")
    validate_cc94_request(parent)
    if parent.heterogen_contexts or parent.retained_explicit_links:
        raise DiffusionContractError("complete-context pilot supports canonical protein without explicit covalent context")
    points = source_coordinates(parent, source)
    if set(points) != set(parent.fixed_atoms):
        raise DiffusionContractError("complete context source does not cover every observed atom")
    if any(a.atom_name not in ATOM37_ORDER for a in parent.fixed_atoms):
        raise DiffusionContractError("an observed atom cannot be represented in atom37; refusing to drop it")
    expected_names = {r: _CANONICAL_ONE_TO_THREE[target.sequence[i]]
        for target in parent.target_sequences for placement in parent.sequence_placements
        if target.chain == placement.chain
        for i, r in zip(placement.observed_target_indices, placement.observed_residues, strict=True)}
    observed_residues = set()
    observed_atoms = set()
    for line in source.decode().splitlines():
        if not line.startswith("ATOM  "):
            continue
        r = ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())
        if r not in expected_names or line[17:20].strip() != expected_names[r]:
            raise DiffusionContractError("complete source chemistry/target placement differs")
        observed_residues.add(r)
        observed_atoms.add(AtomIdentity(r.chain, r.residue_number, r.insertion_code, line[12:16].strip()))
    if observed_residues != set(expected_names):
        raise DiffusionContractError("complete source omits an observed target residue")
    if observed_atoms != set(parent.fixed_atoms):
        raise DiffusionContractError("complete source contains untracked observed atoms")
    gaps = tuple(parent.gaps[i] for i in gap_indices)
    generated_residues = {r for gap in gaps for r in gap.generated_residues}
    generated = tuple(a for a in parent.generated_atoms
        if ResidueIdentity(a.chain, a.residue_number, a.insertion_code) in generated_residues)
    spans = []
    for placement in parent.sequence_placements:
        indices = set(placement.observed_target_indices)
        indices.update(i for g in gaps if g.chain == placement.chain
                       for i in range(g.target_interval.start, g.target_interval.stop))
        for index in sorted(indices):
            if spans and spans[-1].chain == placement.chain and spans[-1].stop == index:
                spans[-1] = ProteinContextSpan(placement.chain, spans[-1].start, index + 1)
            else:
                spans.append(ProteinContextSpan(placement.chain, index, index + 1))
    spans = _split_native_segments(parent, source, points, spans)
    if not spans or len(spans) > 26:
        raise DiffusionContractError("complete context exceeds the native entity budget; no local fallback")
    tokens = sum(s.stop - s.start for s in spans)
    junctions = {r for gap in gaps for r in (gap.generated_residues[0], gap.generated_residues[-1])}
    plan = ProteinContextPlan(
        policy=ProteinContextPolicy(name=COMPLETE_REVISION, target_flanks=None,
            partner_mode="full", residue_budget=tokens),  # requirement, NOT execution permission
        source_sha256=parent.normalized_pdb.sha256, spans=tuple(spans),
        generated=tuple(a for a in generated if ResidueIdentity(a.chain, a.residue_number, a.insertion_code) not in junctions),
        movable_junction=tuple(a for a in generated if ResidueIdentity(a.chain, a.residue_number, a.insertion_code) in junctions),
        fixed_conditioning=tuple(ContextCoordinate(a, tuple(float(x) for x in points[a])) for a in parent.fixed_atoms),
        preserved_only=())
    options = tuple(o for o in parent.backend_options if o.name != "cc89_gap_execution") + (
        BackendOption("cc89_gap_execution", COMPLETE_REVISION),)
    group = replace(parent, gaps=gaps, generated_atoms=generated, protein_context=plan, gap_contexts=None,
        seeds=(seed,), candidate_count=1, backend_options=options,
        regenerated_template_ownership=(tuple(o for o in parent.regenerated_template_ownership
            if o.residue in generated_residues) if parent.regenerated_template_ownership is not None else None))
    view = engine_context_view(group)
    positions, chains = source_context_position_features(group)
    if {s for _, s in view.atom_to_source} != set(parent.fixed_atoms) | set(generated):
        raise DiffusionContractError("complete engine mapping lost source or generated atoms")
    return JointGapResearchView(group, view, gap_indices,
        hashlib.sha256(parent.to_json().encode()).hexdigest(), positions, chains)


def complete_context_admission(view: JointGapResearchView):
    """Fail closed, including when a complete plan exceeds the allowed budget."""
    plan = view.request.protein_context
    reasons = []
    if plan is None or plan.policy.name != COMPLETE_REVISION:
        reasons.append("complete-context-revision-required")
    else:
        if plan.preserved_only or {c.atom for c in plan.fixed_conditioning} != set(view.request.fixed_atoms):
            reasons.append("observed-context-incomplete")
        represented = {s for _, s in view.engine_view.atom_to_source}
        if not set(view.request.fixed_atoms) <= represented:
            reasons.append("engine-observed-context-incomplete")
    if len(view.positions) > COMPLETE_EXECUTION_TOKEN_LIMIT:
        reasons.append("complete-context-exceeds-explicit-research-token-budget")
    return {"supported": not reasons, "reasons": reasons, "required_tokens": len(view.positions),
        "execution_token_limit": COMPLETE_EXECUTION_TOKEN_LIMIT,
        "public_token_limit": 512, "public_profile_enabled": False,
        "required_fixed_atoms": len(view.request.fixed_atoms),
        "conditioning_fixed_atoms": len(plan.fixed_conditioning) if plan else 0,
        "output_only_observed_atoms": len(plan.preserved_only) if plan else 0}


def build_joint_gap_view(parent, source: bytes, gap_indices: tuple[int, ...], *, seed: int):
    """Union independently audited source contexts without any axis truncation.

    This pilot supports 2-4 internal gaps in one physical chain and <=512 tokens.
    Context is exactly the union of source-selected spans; no predicted reference
    coordinates are used. A private policy name prevents public reinterpretation.
    """
    return _build_joint_view(parent, source, gap_indices, seed=seed, whole_structure=False)


def build_joint_structure_view(parent, source: bytes, *, seed: int):
    """Explicit full-structure research extension, <=1024 tokens/16 gaps.

    512 is the public audited budget and training fixed_size, NOT relaxed here.
    This distinct private operation uses the pinned model's variable-length
    sampling axis; callers must exercise live memory guards and native probes.
    No context is removed to fit, and the checkpoint/config remain unchanged.
    """
    return _build_joint_view(parent, source, tuple(range(len(parent.gaps))), seed=seed, whole_structure=True)


def plan_motif_component(parent, source: bytes, gap_indices: tuple[int, ...], *, seed: int):
    """Jointly sample exactly one source-neighbor component using local motifs.

    No complete-complex conditioning and no independent fallback inside a
    component. Different components never condition on each other's predictions.
    Planning may exceed the budget; motif_component_admission controls execution.
    """
    from dvbfixer.model.diffusion.gap_neighbors import source_gap_components

    components, _ = source_gap_components(parent, source.decode())
    if gap_indices not in components:
        raise DiffusionContractError("motif selection must equal one complete neighboring-gap component")
    return _build_joint_view(parent, source, gap_indices, seed=seed,
        whole_structure=False, motif_component=True)


def motif_component_admission(view: JointGapResearchView):
    plan = view.request.protein_context
    reasons = []
    if plan is None or plan.policy.name != MOTIF_COMPONENT_REVISION:
        reasons.append("motif-component-revision-required")
    if len(view.positions) > MOTIF_COMPONENT_TOKEN_LIMIT:
        reasons.append("neighbor-component-exceeds-motif-token-budget-no-splitting")
    return {"supported": not reasons, "reasons": reasons,
        "required_tokens": len(view.positions), "execution_token_limit": MOTIF_COMPONENT_TOKEN_LIMIT,
        "public_token_limit": 512, "public_profile_enabled": False,
        "experimental_over_public_budget": len(view.positions) > 512,
        "conditioning_fixed_atoms": len(plan.fixed_conditioning) if plan else 0,
        "output_only_observed_atoms": len(plan.preserved_only) if plan else 0}


def _build_joint_view(parent, source, gap_indices, *, seed, whole_structure, motif_component=False):
    revision = (MOTIF_COMPONENT_REVISION if motif_component else
                STRUCTURE_REVISION if whole_structure else REVISION)
    budget = (sum(len(t.sequence) for t in parent.target_sequences) if motif_component else
              1024 if whole_structure else 512)
    maximum_gaps = 64 if motif_component else 16 if whole_structure else 4
    if (not (1 if motif_component else 2) <= len(gap_indices) <= maximum_gaps or tuple(sorted(set(gap_indices))) != gap_indices
            or any(type(i) is not int or not 0 <= i < len(parent.gaps) for i in gap_indices)
            or type(seed) is not int or seed < 0):
        raise DiffusionContractError("joint research requires bounded ordered distinct gaps and one explicit seed")
    if hashlib.sha256(source).hexdigest() != parent.normalized_pdb.sha256:
        raise DiffusionContractError("joint research source digest differs")
    validate_cc94_request(parent)
    if parent.gap_contexts is None or parent.heterogen_contexts:
        raise DiffusionContractError("joint research needs frozen protein-only source contexts")
    gaps = tuple(parent.gaps[i] for i in gap_indices)
    if not motif_component and (len({g.chain for g in gaps}) != 1 or
            (not whole_structure and any(g.gap_kind is not GapKind.INTERNAL for g in gaps))):
        raise DiffusionContractError("joint research pilot requires internal gaps on one physical chain")
    singles = independent_gap_requests(parent)
    plans = []
    for i in gap_indices:
        validate_context_source(singles[i], source)
        plans.append(parent.gap_contexts[i])
    if any(plan.policy != plans[0].policy for plan in plans):
        raise DiffusionContractError("joint research source-context policies differ")
    chain_order = tuple(target.chain for target in parent.target_sequences)
    indices = {(s.chain, i) for p in plans for s in p.spans for i in range(s.start, s.stop)}
    if motif_component:
        # Other unknown gaps are neither fixed motifs nor unowned dummy tokens.
        indices.difference_update((g.chain, i) for j, g in enumerate(parent.gaps) if j not in gap_indices
            for i in range(g.target_interval.start, g.target_interval.stop))
        budget = len(indices)  # required shape, not execution permission
    spans = []
    for chain in chain_order:
        for index in sorted(i for c, i in indices if c == chain):
            if spans and spans[-1].chain == chain and spans[-1].stop == index:
                spans[-1] = ProteinContextSpan(chain, spans[-1].start, index + 1)
            else:
                spans.append(ProteinContextSpan(chain, index, index + 1))
    if len(indices) > budget or len(spans) > 26:
        raise DiffusionContractError("joint source union exceeds the token/entity budget; no independent fallback")
    spans = _split_native_segments(parent, source, source_coordinates(parent, source), spans)
    if len(spans) > 26:
        raise DiffusionContractError("joint source union exceeds native entity budget")
    generated_residues = {r for gap in gaps for r in gap.generated_residues}
    generated = tuple(a for a in parent.generated_atoms
        if ResidueIdentity(a.chain, a.residue_number, a.insertion_code) in generated_residues)
    mobile = {a for p in plans for a in p.movable_junction}
    fixed = {c.atom: c for p in plans for c in p.fixed_conditioning}
    if not set(fixed) <= set(parent.fixed_atoms) or set(fixed) & set(parent.generated_atoms):
        raise DiffusionContractError("joint conditioning must contain original observed atoms only")
    plan = ProteinContextPlan(
        policy=replace(plans[0].policy, name=revision, residue_budget=budget), source_sha256=parent.normalized_pdb.sha256,
        spans=tuple(spans), generated=tuple(a for a in generated if a not in mobile),
        movable_junction=tuple(a for a in generated if a in mobile),
        fixed_conditioning=tuple(fixed[a] for a in parent.fixed_atoms if a in fixed),
        preserved_only=tuple(a for a in parent.fixed_atoms if a not in fixed))
    options = tuple(o for o in parent.backend_options if o.name != "cc89_gap_execution") + (
        BackendOption("cc89_gap_execution", revision),)
    group = replace(parent, gaps=gaps, generated_atoms=generated, gap_contexts=None, protein_context=plan,
        seeds=(seed,), candidate_count=1, backend_options=options,
        regenerated_template_ownership=(tuple(o for o in parent.regenerated_template_ownership
            if o.residue in generated_residues) if parent.regenerated_template_ownership is not None else None))
    view = engine_context_view(group)
    positions, chains = source_context_position_features(group)
    if len(positions) != sum(len(t.sequence) for t in view.request.target_sequences):
        raise DiffusionContractError("joint engine axis differs from physical source axis")
    return JointGapResearchView(group, view, gap_indices,
        hashlib.sha256(parent.to_json().encode()).hexdigest(), positions, chains)
