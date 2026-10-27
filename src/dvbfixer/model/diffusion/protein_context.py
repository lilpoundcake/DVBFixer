"""Deterministic, source-only protein context selection and engine span mapping."""

from __future__ import annotations

import hashlib
import json
import string
from dataclasses import dataclass, replace

import numpy as np

from dvbfixer.domain.structure_identity import allocate_chain_ids
from dvbfixer.model.diffusion.cc94 import PROFILES as CC94_PROFILES
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    BackendOption,
    ContextCoordinate,
    DiffusionRequest,
    GapKind,
    ProteinContextEvidence,
    ProteinContextPlan,
    ProteinContextPolicy,
    ProteinContextSpan,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
    TemplateOwnership,
    atom_identity_digest,
    protein_context_digest,
)
from dvbfixer.model.diffusion.scope import (
    _CANONICAL_ONE_TO_THREE,
    CANONICAL_HEAVY_ATOMS,
    canonical_heavy_atom_identities,
)


class ProteinContextError(ValueError):
    """Required context is incomplete, incompatible, or exceeds the engine budget."""


CUDA_CC89_PROFILE = "protpardelle-1c-cuda"
MPS_CC89_PROFILE = "protpardelle-1c-mps-local-v1"
CUDA_CC89_CONTEXT = ProteinContextPolicy(
    name="atom37-terminal-preservation-v1", target_flanks=8,
    same_chain_spatial=True, partner_mode="local",
)
ANCHOR_SEGMENT_OPTION = BackendOption("cc89_context_geometry", "anchor-shell-segment-union-v1")
ANCHOR_SEGMENT_CONTEXT = replace(CUDA_CC89_CONTEXT, name=ANCHOR_SEGMENT_OPTION.value)
TERMINAL_SHELL_OPTION = BackendOption("cc89_terminal_context", "eight-observed-residues-20a-union-v1")
TERMINAL_SHELL_CONTEXT = replace(CUDA_CC89_CONTEXT, name="anchor-shell-segment-terminal20-union-v1")


def _cc89_context_policy(request: DiffusionRequest) -> ProteinContextPolicy:
    option = next((item for item in request.backend_options if item.name == ANCHOR_SEGMENT_OPTION.name), None)
    terminal = next((item for item in request.backend_options if item.name == TERMINAL_SHELL_OPTION.name), None)
    if terminal is not None:
        if terminal != TERMINAL_SHELL_OPTION or option != ANCHOR_SEGMENT_OPTION:
            raise ProteinContextError("terminal context expansion requires explicit matching source geometry options")
        return TERMINAL_SHELL_CONTEXT
    if option is None:
        return CUDA_CC89_CONTEXT
    if option != ANCHOR_SEGMENT_OPTION:
        raise ProteinContextError("unsupported cc89 context geometry")
    return ANCHOR_SEGMENT_CONTEXT


CUDA_CC89_OPTIONS = (
    BackendOption("cc89_chain_phase", "neutral-v1"),
    BackendOption("cc89_position_features", "source-physical-v1"),
)
INDEPENDENT_GAPS_OPTION = BackendOption("cc89_gap_execution", "independent-local-v1")


def independent_gap_requests(request: DiffusionRequest) -> tuple[DiffusionRequest, ...]:
    """Immutable per-gap views: other missing regions never become conditioning."""
    contexts = request.gap_contexts or (None,) * len(request.gaps)
    result = []
    for gap, plan in zip(request.gaps, contexts, strict=True):
        residues = set(gap.generated_residues)
        result.append(replace(
            request, gaps=(gap,), gap_contexts=None, protein_context=plan,
            generated_atoms=tuple(atom for atom in request.generated_atoms if _residue(atom) in residues),
            regenerated_template_ownership=(tuple(owner for owner in request.regenerated_template_ownership if owner.residue in residues)
                if request.regenerated_template_ownership is not None else None),
        ))
    return tuple(result)


def is_cuda_cc89(request: DiffusionRequest) -> bool:
    return any(option.name == "profile" and option.value == CUDA_CC89_PROFILE for option in request.backend_options)


def is_mps_cc89(request: DiffusionRequest) -> bool:
    return any(option.name == "profile" and option.value == MPS_CC89_PROFILE for option in request.backend_options)


def is_corrected_cc89(request: DiffusionRequest) -> bool:
    """Compatibility name for the shared source-only context family, not weights.

    Select numerical interventions with is_cuda_cc89/is_mps_cc89/is_cc94 instead.
    """
    return is_cuda_cc89(request) or is_mps_cc89(request) or is_cc94(request)


def is_cc94(request: DiffusionRequest) -> bool:
    return any(option.name == "profile" and option.value in CC94_PROFILES for option in request.backend_options)


def _position_options(request: DiffusionRequest) -> tuple[BackendOption, ...]:
    if is_cc94(request):
        return (BackendOption("cc89_chain_phase", "native"), CUDA_CC89_OPTIONS[1])
    return CUDA_CC89_OPTIONS


def prepare_cc94_request(request: DiffusionRequest, pdb: bytes, context_mode: str) -> DiffusionRequest:
    """Reuse source-only local selection with native relative position encoding."""
    if not is_cc94(request):
        raise ProteinContextError("cc94 request must declare its versioned profile")
    return prepare_mps_cc89_request(request, pdb, context_mode)


def validate_cc94_request(request: DiffusionRequest) -> None:
    if not is_cc94(request):
        raise ProteinContextError("cc94 request must declare its versioned profile")
    profile = next(o.value for o in request.backend_options if o.name == "profile")
    _validate_corrected_cc89_request(request, profile, "cc94")


def prepare_cuda_cc89_request(request: DiffusionRequest, pdb: bytes, context_mode: str) -> DiffusionRequest:
    """Apply the public CUDA profile's declared local context, never a fallback."""
    return _prepare_corrected_cc89_request(request, pdb, context_mode, CUDA_CC89_PROFILE, "CUDA cc89")


def prepare_mps_cc89_request(request: DiffusionRequest, pdb: bytes, context_mode: str) -> DiffusionRequest:
    """Prepare the versioned Apple profile without relabeling a CUDA request."""
    independent = (len(request.gaps) > 1 or request.gaps[0].gap_kind is not GapKind.INTERNAL
        or sum(len(target.sequence) for target in request.target_sequences) > 512
        or INDEPENDENT_GAPS_OPTION in request.backend_options)
    if independent:
        if len(request.gaps) > 64:
            raise ProteinContextError("independent local execution supports at most 64 gaps per candidate")
        if context_mode not in {"auto", "local"}:
            raise ProteinContextError("independent local gaps require --diffusion-context auto|local")
        if request.protein_context is not None:
            if len(request.gaps) != 1 or request.gap_contexts is not None:
                raise ProteinContextError("independent gaps cannot reinterpret a joint context plan")
            # A frozen per-gap view is already materialized. Verify it exactly;
            # do not wrap it again in a parent execution plan or change its hash.
            (validate_cc94_request if is_cc94(request) else validate_mps_cc89_request)(request)
            validate_context_source(request, pdb)
            return request
        options = list(request.backend_options)
        for option in (*_position_options(request), INDEPENDENT_GAPS_OPTION):
            existing = next((item for item in options if item.name == option.name), None)
            if existing is not None and existing != option:
                raise ProteinContextError("independent gap execution policy conflicts with the request")
            if existing is None:
                options.append(option)
        request = replace(request, backend_options=tuple(options))
        if request.gap_contexts is None:
            plans = []
            for single in independent_gap_requests(request):
                plan = select_protein_context(single, pdb, _cc89_context_policy(request)).protein_context
                assert plan is not None
                plans.append(plan)
            request = replace(request, gap_contexts=tuple(plans))
        (validate_cc94_request if is_cc94(request) else validate_mps_cc89_request)(request)
        validate_context_source(request, pdb)
        return request
    profile = next(o.value for o in request.backend_options if o.name == "profile") if is_cc94(request) else MPS_CC89_PROFILE
    return _prepare_corrected_cc89_request(request, pdb, context_mode, profile, "cc94" if is_cc94(request) else "MPS cc89")


def _prepare_corrected_cc89_request(
    request: DiffusionRequest, pdb: bytes, context_mode: str, profile: str, label: str,
) -> DiffusionRequest:
    if not any(option.name == "profile" and option.value == profile for option in request.backend_options):
        raise ProteinContextError(f"{label} request must declare its profile")
    if context_mode not in {"auto", "local"}:
        raise ProteinContextError(f"{profile} requires --diffusion-context auto|local")
    if len(request.gaps) != 1 or request.gaps[0].gap_kind is not GapKind.INTERNAL:
        raise ProteinContextError(f"{label} supports one internal gap per request")
    if request.protein_context is None:
        request = select_protein_context(request, pdb, _cc89_context_policy(request))
    options = list(request.backend_options)
    for required in _position_options(request):
        existing = next((option for option in options if option.name == required.name), None)
        if existing is not None and existing != required:
            raise ProteinContextError(f"{label} numerical policy conflicts with the frozen profile")
        if existing is None:
            options.append(required)
    request = replace(request, backend_options=tuple(options))
    _validate_corrected_cc89_request(request, profile, label)
    validate_context_source(request, pdb)
    return request


def validate_cuda_cc89_request(request: DiffusionRequest) -> None:
    """Validate mandatory profile limits, independent of length evidence."""
    _validate_corrected_cc89_request(request, CUDA_CC89_PROFILE, "CUDA cc89")


def validate_mps_cc89_request(request: DiffusionRequest) -> None:
    _validate_corrected_cc89_request(request, MPS_CC89_PROFILE, "MPS cc89")


def _validate_corrected_cc89_request(request: DiffusionRequest, profile: str, label: str) -> None:
    if not any(option.name == "profile" and option.value == profile for option in request.backend_options):
        raise ProteinContextError(f"{label} request must declare its profile")
    independent = (profile == MPS_CC89_PROFILE or profile in CC94_PROFILES) and INDEPENDENT_GAPS_OPTION in request.backend_options
    if request.gap_contexts is not None:
        if not independent:
            raise ProteinContextError("independent gap contexts require their versioned execution policy")
        for single in independent_gap_requests(request):
            _validate_corrected_cc89_request(single, profile, label)
        return
    if len(request.gaps) != 1 or (not independent and request.gaps[0].gap_kind is not GapKind.INTERNAL):
        raise ProteinContextError(f"{label} supports one internal gap per request")
    # Length is an evidence boundary, reported by scope admission as a warning.
    # Keep topology, chemistry, context and total-resource limits mandatory.
    gap = request.gaps[0]
    sequence = next(target.sequence for target in request.target_sequences if target.chain == gap.chain)
    if not independent and "C" in sequence[gap.target_interval.start:gap.target_interval.stop]:
        raise ProteinContextError(f"{label} gap support is limited to CYS-free windows")
    if request.heterogen_contexts:
        raise ProteinContextError(f"{label} supports canonical protein inputs only")
    if not independent and sum(len(target.sequence) for target in request.target_sequences) > 512:
        raise ProteinContextError(f"{label} whole-source budget is 512 protein residues")
    plan = request.protein_context
    if (plan is None or plan.policy != _cc89_context_policy(request) or plan.trajectory_queries is not None
            or plan.retained_conditioning or plan.previous_plan_sha256 is not None):
        raise ProteinContextError(f"{label} requires the frozen source-only local context policy")
    if any(option not in request.backend_options for option in _position_options(request)):
        raise ProteinContextError(f"{label} requires its declared chain phase and physical source positions")
    source_context_position_features(request)
    if sum(span.stop - span.start for span in plan.spans) > 512:
        raise ProteinContextError("each independent local engine axis must fit 512 residues")


def source_coordinates(request: DiffusionRequest, pdb: bytes) -> dict[AtomIdentity, np.ndarray]:
    if hashlib.sha256(pdb).hexdigest() != request.normalized_pdb.sha256:
        raise ProteinContextError("context source digest mismatch")
    try:
        text = pdb.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ProteinContextError("context source is not valid UTF-8") from exc
    fixed = set(request.fixed_atoms)
    points = {}
    for line in text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        try:
            atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        except (ValueError, IndexError) as exc:
            raise ProteinContextError("malformed context source atom identity") from exc
        if atom not in fixed:
            continue
        if atom in points:
            raise ProteinContextError("duplicate context source atom identity")
        try:
            coordinate = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError as exc:
            raise ProteinContextError("malformed context source coordinate") from exc
        if not np.isfinite(coordinate).all():
            raise ProteinContextError("non-finite context source coordinate")
        points[atom] = coordinate
    if set(points) != set(request.fixed_atoms):
        raise ProteinContextError("protein context requires every fixed atom in the normalized protein source")
    return points


def _residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


def select_protein_context(
    request: DiffusionRequest, pdb: bytes, policy: ProteinContextPolicy, *,
    trajectory_queries: tuple[ContextCoordinate, ...] | None = None,
    retained_conditioning: tuple[AtomIdentity, ...] | None = None,
    previous_plan_sha256: str | None = None,
    adaptation_source: ArtifactReference | None = None,
) -> DiffusionRequest:
    """Select complete spans from observed anchor atoms; never inspect withheld positions."""
    targets = {target.chain: target for target in request.target_sequences}
    placements = {placement.chain: placement for placement in request.sequence_placements}
    gap_chains = {gap.chain for gap in request.gaps}
    if not gap_chains or request.heterogen_contexts:
        raise ProteinContextError("protein context requires gap-bearing chains and canonical protein context")
    points = source_coordinates(request, pdb)
    atoms_by_residue: dict[ResidueIdentity, list[AtomIdentity]] = {}
    for atom in request.fixed_atoms:
        atoms_by_residue.setdefault(_residue(atom), []).append(atom)
    anchors = {anchor for gap in request.gaps for anchor in (gap.left_anchor, gap.right_anchor) if anchor is not None}
    anchor_points = np.asarray([points[atom] for residue in sorted(anchors) for atom in atoms_by_residue.get(residue, [])])
    if len(anchor_points) < 3:
        raise ProteinContextError("context requires complete observed anchor geometry")
    anchor_segment = None
    has_segment_union = policy.name in {ANCHOR_SEGMENT_CONTEXT.name, TERMINAL_SHELL_CONTEXT.name}
    if has_segment_union:
        if len(request.gaps) != 1:
            raise ProteinContextError("anchor-segment context requires one gap per sampling view")
    if has_segment_union and request.gaps[0].gap_kind is GapKind.INTERNAL:
        gap = request.gaps[0]
        assert gap.left_anchor is not None and gap.right_anchor is not None
        endpoints = []
        for anchor in (gap.left_anchor, gap.right_anchor):
            point = points.get(AtomIdentity(anchor.chain, anchor.residue_number, anchor.insertion_code, "CA"))
            if point is None:
                raise ProteinContextError("anchor-segment context requires both observed CA atoms")
            endpoints.append(point)
        anchor_segment = (endpoints[0], endpoints[1])
    terminal_points = None
    if policy.name == TERMINAL_SHELL_CONTEXT.name and request.gaps[0].gap_kind is not GapKind.INTERNAL:
        gap = request.gaps[0]
        observed = placements[gap.chain].observed_residues
        anchor = gap.right_anchor if gap.gap_kind is GapKind.N_TERMINAL else gap.left_anchor
        assert anchor is not None
        anchor_index = observed.index(anchor)
        flank = (observed[anchor_index:anchor_index + 8] if gap.gap_kind is GapKind.N_TERMINAL
                 else observed[max(0, anchor_index - 7):anchor_index + 1])
        # Only measured atoms from the terminal observed window, never from
        # generated tails or another invocation's predictions.
        terminal_points = np.asarray([points[atom] for residue in flank for atom in atoms_by_residue.get(residue, [])])
    if trajectory_queries is not None:
        expected = {atom for atom in request.generated_atoms if atom.atom_name in {"N", "CA", "C", "O"}}
        if {query.atom for query in trajectory_queries} != expected:
            raise ProteinContextError("adaptive selection requires the complete generated backbone, not reference positions")
        anchor_points = np.concatenate((anchor_points, np.asarray([query.coordinate_angstrom for query in trajectory_queries])))
    selected: set[tuple[str, int]] = set()
    for gap in request.gaps:
        start = 0 if policy.target_flanks is None else max(0, gap.target_interval.start - policy.target_flanks)
        stop = len(targets[gap.chain].sequence) if policy.target_flanks is None else min(len(targets[gap.chain].sequence), gap.target_interval.stop + policy.target_flanks)
        selected.update((gap.chain, index) for index in range(start, stop))
    if INDEPENDENT_GAPS_OPTION in request.backend_options:
        available = {(placement.chain, index) for placement in request.sequence_placements
            for index in placement.observed_target_indices}
        available.update((gap.chain, index) for gap in request.gaps
            for index in range(gap.target_interval.start, gap.target_interval.stop))
        selected.intersection_update(available)
    if retained_conditioning is not None:
        if not set(retained_conditioning) <= set(request.fixed_atoms):
            raise ProteinContextError("adaptive retained identities must belong to the immutable fixed source")
        retained_residues = {_residue(atom) for atom in retained_conditioning}
        selected.update((placement.chain, index) for placement in request.sequence_placements
                        for index, residue in zip(placement.observed_target_indices, placement.observed_residues) if residue in retained_residues)
    priorities = []
    for chain, placement in placements.items():
        if chain not in gap_chains and policy.partner_mode == "none":
            continue
        for index, residue in zip(placement.observed_target_indices, placement.observed_residues):
            if (chain, index) in selected:
                continue
            if chain in gap_chains and not policy.same_chain_spatial:
                continue
            residue_atoms = atoms_by_residue.get(residue, [])
            if not residue_atoms:
                continue
            residue_points = np.asarray([points[atom] for atom in residue_atoms])
            distance = float(np.min(np.linalg.norm(residue_points[:, None, :] - anchor_points[None, :, :], axis=2)))
            if anchor_segment is not None:
                # Union with the existing atom-based anchor shells, not a
                # replacement. Clamp projection to the finite CA–CA segment;
                # a residue beyond an endpoint is not near an infinite line.
                distance = min(distance, float(np.min(_point_segment_distances(residue_points, *anchor_segment))))
            full_partner = chain not in gap_chains and policy.partner_mode == "full"
            terminal_contact = (terminal_points is not None and
                float(np.min(np.linalg.norm(residue_points[:, None, :] - terminal_points[None, :, :], axis=2))) <= 20.0)
            if not full_partner and distance > policy.secondary_radius_angstrom and not terminal_contact:
                continue
            rank = 0 if distance <= policy.direct_radius_angstrom else 1
            priorities.append((rank, distance, chain, index))
    # Include every required shell contact and complete neighboring sequence
    # windows. Priority is deterministic; budget failure never trims contacts.
    for _rank, _distance, chain, index in sorted(priorities):
        if policy.partner_mode == "full" and chain not in gap_chains:
            selected.update((chain, i) for i in range(len(targets[chain].sequence)))
        else:
            observed = set(placements[chain].observed_target_indices)
            left = index
            right = index + 1
            for _ in range(policy.sequence_padding):
                if left - 1 in observed:
                    left -= 1
                if right in observed:
                    right += 1
            selected.update((chain, i) for i in range(left, right))
    spans: list[ProteinContextSpan] = []
    # Retain declared target order for gap-bearing entities and deterministic
    # partner order. Single-chain legacy plans retain their exact span ordering.
    for chain in (*(target.chain for target in request.target_sequences if target.chain in gap_chains), *sorted(set(targets) - gap_chains)):
        indices = sorted(index for c, index in selected if c == chain)
        for index in indices:
            if spans and spans[-1].chain == chain and spans[-1].stop == index:
                spans[-1] = ProteinContextSpan(chain, spans[-1].start, index + 1)
            else:
                spans.append(ProteinContextSpan(chain, index, index + 1))
    if policy.name == "intact-native-segments-v1":
        spans = _split_native_segments(request, pdb, points, spans)
    if sum(span.stop - span.start for span in spans) > policy.residue_budget:
        raise ProteinContextError("required protein context exceeds residue budget; no neighbors were silently dropped")
    selected_residues = set()
    for chain, placement in placements.items():
        for index, residue in zip(placement.observed_target_indices, placement.observed_residues):
            if (chain, index) not in selected:
                continue
            expected = set(canonical_heavy_atom_identities(residue, _CANONICAL_ONE_TO_THREE[targets[chain].sequence[index]]))
            observed_atoms = set(atoms_by_residue.get(residue, []))
            if not expected <= observed_atoms:
                backbone = {AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, name)
                    for name in ("N", "CA", "C", "O")}
                if INDEPENDENT_GAPS_OPTION not in request.backend_options or not backbone <= observed_atoms:
                    raise ProteinContextError(f"selected context residue is heavy-atom-incomplete: {residue}")
                # Sparse atom37 motif masks represent only measured sidechain
                # atoms. Never invent conditioning coordinates for absent atoms.
            selected_residues.add(residue)
    preserved_terminal_atoms: set[AtomIdentity] = set()
    if policy.name in {"atom37-terminal-preservation-v1", ANCHOR_SEGMENT_CONTEXT.name, TERMINAL_SHELL_CONTEXT.name}:
        terminals = {placement.observed_residues[-1] for placement in request.sequence_placements
                     if placement.observed_target_indices and placement.observed_target_indices[-1] == placement.target_length - 1}
        preserved_terminal_atoms = {atom for atom in request.fixed_atoms if atom.atom_name == "OXT"}
        if any(_residue(atom) not in terminals for atom in preserved_terminal_atoms):
            raise ProteinContextError("atom37 terminal preservation requires OXT on an observed C terminus")
    selected_fixed = tuple(atom for atom in request.fixed_atoms
                           if _residue(atom) in selected_residues and atom not in preserved_terminal_atoms)
    junction = {residue for gap in request.gaps for residue in (gap.generated_residues[0], gap.generated_residues[-1])}
    plan = ProteinContextPlan(
        policy, request.normalized_pdb.sha256, tuple(spans),
        tuple(atom for atom in request.generated_atoms if _residue(atom) not in junction),
        tuple(atom for atom in request.generated_atoms if _residue(atom) in junction),
        tuple(ContextCoordinate(atom, (float(points[atom][0]), float(points[atom][1]), float(points[atom][2]))) for atom in selected_fixed),
        tuple(atom for atom in request.fixed_atoms if atom not in selected_fixed),
        trajectory_queries, retained_conditioning, previous_plan_sha256, adaptation_source,
    )
    return replace(request, protein_context=plan)


def _point_segment_distances(points: np.ndarray, start: np.ndarray, stop: np.ndarray) -> np.ndarray:
    """Euclidean distances in Å to a finite segment, including degenerate ends."""
    vector = stop - start
    squared_length = float(np.dot(vector, vector))
    if squared_length == 0:
        return np.linalg.norm(points - start, axis=-1)
    fraction = np.clip((points - start) @ vector / squared_length, 0.0, 1.0)
    nearest = start + fraction[:, None] * vector
    return np.linalg.norm(points - nearest, axis=-1)


def _split_native_segments(
    request: DiffusionRequest,
    pdb: bytes,
    points: dict[AtomIdentity, np.ndarray],
    spans: list[ProteinContextSpan],
) -> list[ProteinContextSpan]:
    """Keep observed native breaks as separate engine entities, not deletions.

    The intact-input research policy uses an observed-residue sequence axis.
    Consecutive axis indices do not imply a peptide bond across a native break.
    Artificial gaps are distinguished by their intervening generated indices.
    """
    indices = {residue: (placement.chain, index) for placement in request.sequence_placements
               for index, residue in zip(placement.observed_target_indices, placement.observed_residues)}
    boundaries: set[tuple[str, int]] = set()
    for placement in request.sequence_placements:
        observed = tuple(zip(placement.observed_target_indices, placement.observed_residues))
        for (left_index, left), (right_index, right) in zip(observed, observed[1:]):
            if right_index != left_index + 1:
                continue
            carbon = points.get(AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C"))
            nitrogen = points.get(AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N"))
            if carbon is None or nitrogen is None or not 1.0 <= float(np.linalg.norm(carbon - nitrogen)) <= 1.8:
                boundaries.add((placement.chain, right_index))
    previous: ResidueIdentity | None = None
    terminated = False
    for line in pdb.decode().splitlines():
        if line.startswith("TER"):
            terminated = True
        elif line.startswith("ATOM  "):
            residue = ResidueIdentity(line[21], line[22:26].strip(), line[26].strip())
            if residue != previous:
                if terminated and previous is not None and previous in indices and residue in indices:
                    left_chain, left_index = indices[previous]
                    right_chain, right_index = indices[residue]
                    if left_chain == right_chain and right_index == left_index + 1:
                        boundaries.add((right_chain, right_index))
                previous = residue
                terminated = False
    separated: list[ProteinContextSpan] = []
    for span in spans:
        cuts = [span.start, *(index for chain, index in sorted(boundaries) if chain == span.chain and span.start < index < span.stop), span.stop]
        separated.extend(ProteinContextSpan(span.chain, left, right) for left, right in zip(cuts, cuts[1:]))
    return separated


@dataclass(frozen=True)
class SamplingEntity:
    engine_chain: str
    target: TargetSequence
    placement: SequencePlacement
    residues: tuple[ResidueIdentity, ...]


@dataclass(frozen=True)
class EngineContextView:
    """Private engine identities plus an exact inverse to the scientific request."""

    request: DiffusionRequest
    atom_to_source: tuple[tuple[AtomIdentity, AtomIdentity], ...]

    def __post_init__(self) -> None:
        engine = [pair[0] for pair in self.atom_to_source]
        source = [pair[1] for pair in self.atom_to_source]
        if len(set(engine)) != len(engine) or len(set(source)) != len(source):
            raise ProteinContextError("engine context mapping must be bijective")
        if set(engine) != set(self.request.fixed_atoms) | set(self.request.generated_atoms):
            raise ProteinContextError("engine context mapping must cover every represented atom")


def sampling_entities(request: DiffusionRequest) -> tuple[SamplingEntity, ...]:
    targets = {target.chain: target for target in request.target_sequences}
    placements = {placement.chain: placement for placement in request.sequence_placements}
    spans = request.protein_context.spans if request.protein_context else tuple(ProteinContextSpan(target.chain, 0, len(target.sequence)) for target in request.target_sequences)
    reserved = {span.chain for span in spans}
    extra_count = len(spans) - len(reserved - {" "})
    extra_labels = iter(allocate_chain_ids(extra_count, alphabet=tuple(string.ascii_uppercase + string.ascii_lowercase + string.digits), reserved=reserved))
    used = set()
    entities = []
    for span in spans:
        engine_chain = span.chain if span.chain not in used and span.chain != " " else next(extra_labels)
        used.add(engine_chain)
        placement = placements[span.chain]
        mapped = dict(zip(placement.observed_target_indices, placement.observed_residues))
        for gap in request.gaps:
            if gap.chain == span.chain:
                mapped.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues))
        if any(index not in mapped for index in range(span.start, span.stop)):
            raise ProteinContextError("sampling span contains an unmapped sequence position")
        observed = [(index - span.start, residue) for index, residue in zip(placement.observed_target_indices, placement.observed_residues) if span.start <= index < span.stop]
        entities.append(SamplingEntity(
            engine_chain, TargetSequence(span.chain, targets[span.chain].sequence[span.start:span.stop]),
            SequencePlacement(span.chain, span.stop - span.start, tuple(index for index, _residue_id in observed), tuple(residue for _index, residue in observed)),
            tuple(mapped[index] for index in range(span.start, span.stop)),
        ))
    return tuple(entities)


def source_context_position_features(
    request: DiffusionRequest, *, chain_gap: int = 200,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Position a packed context axis using physical source chains, not proxies.

    Motif serialization may split one chain into several entities. Their source
    sequence separations must survive packing. Different physical chains retain
    an explicit position break; omitted residues never become peptide neighbors.
    """
    plan = request.protein_context
    if plan is None:
        raise ProteinContextError("source position features require an explicit context plan")
    if type(chain_gap) is not int or chain_gap < 2:
        raise ProteinContextError("source chains require a non-peptide position break")
    order = tuple(dict.fromkeys(span.chain for span in plan.spans))
    seen: set[str] = set()
    previous_chain: str | None = None
    ends: dict[str, int] = {}
    for span in plan.spans:
        if span.chain != previous_chain and span.chain in seen:
            raise ProteinContextError("source context spans must be grouped by physical chain")
        if span.chain in ends and span.start < ends[span.chain]:
            raise ProteinContextError("source context spans must retain sequence order")
        seen.add(span.chain)
        ends[span.chain] = span.stop
        previous_chain = span.chain
    origins: dict[str, int] = {}
    previous_end: int | None = None
    for chain in order:
        spans = tuple(span for span in plan.spans if span.chain == chain)
        first = min(span.start for span in spans)
        last = max(span.stop for span in spans)
        start = 1 if previous_end is None else previous_end + chain_gap
        origins[chain] = start - first
        previous_end = origins[chain] + last - 1
    chain_numbers = {chain: index for index, chain in enumerate(order)}
    indices = tuple(origins[span.chain] + index for span in plan.spans for index in range(span.start, span.stop))
    chains = tuple(chain_numbers[span.chain] for span in plan.spans for _ in range(span.start, span.stop))
    return indices, chains


def engine_context_view(request: DiffusionRequest) -> EngineContextView:
    """Represent discontiguous spans as separate real chains, never concatenations.

    Only this internal view uses proxy chain IDs. Scientific source identities,
    including insertion codes and case-distinct chains, round trip through the map.
    """
    if request.protein_context is None:
        raise ProteinContextError("engine context view requires an explicit plan")
    entities = sampling_entities(request)
    selected = {item.atom for item in request.protein_context.fixed_conditioning}
    residue_map = {}
    position_map = {}
    for entity, span in zip(entities, request.protein_context.spans, strict=True):
        for index, residue in enumerate(entity.residues):
            residue_map[residue] = ResidueIdentity(entity.engine_chain, residue.residue_number, residue.insertion_code)
            position_map[(span.chain, span.start + index)] = (entity.engine_chain, index)

    def atom_map(atom: AtomIdentity) -> AtomIdentity:
        residue = residue_map[_residue(atom)]
        return AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom.atom_name)

    gaps = []
    for gap in request.gaps:
        chain, start = position_map[(gap.chain, gap.target_interval.start)]
        if any(position_map[(gap.chain, index)] != (chain, start + index - gap.target_interval.start)
               for index in range(gap.target_interval.start, gap.target_interval.stop)):
            raise ProteinContextError("one generated gap must fit completely within one engine entity")
        gaps.append(replace(gap, chain=chain, target_interval=type(gap.target_interval)(start, start + len(gap.generated_residues)),
                            left_anchor=residue_map[gap.left_anchor] if gap.left_anchor else None,
                            right_anchor=residue_map[gap.right_anchor] if gap.right_anchor else None,
                            generated_residues=tuple(residue_map[residue] for residue in gap.generated_residues),
                            movable_junction_residues=tuple(residue_map[residue] for residue in gap.movable_junction_residues)))

    def owner_map(owner: TemplateOwnership) -> TemplateOwnership:
        chain, index = position_map[(owner.chain, owner.target_index)]
        return replace(owner, chain=chain, target_index=index, residue=residue_map[owner.residue])

    represented = tuple(atom for atom in request.fixed_atoms if atom in selected) + request.generated_atoms
    internal = replace(request, protein_context=None,
        target_sequences=tuple(TargetSequence(entity.engine_chain, entity.target.sequence) for entity in entities),
        sequence_placements=tuple(SequencePlacement(entity.engine_chain, entity.placement.target_length,
            entity.placement.observed_target_indices, tuple(residue_map[residue] for residue in entity.placement.observed_residues)) for entity in entities),
        gaps=tuple(gaps), fixed_atoms=tuple(atom_map(atom) for atom in request.fixed_atoms if atom in selected),
        generated_atoms=tuple(atom_map(atom) for atom in request.generated_atoms), retained_explicit_links=(),
        template_ownership=tuple(owner_map(owner) for owner in request.template_ownership if owner.residue in {_residue(atom) for atom in selected}),
        regenerated_template_ownership=tuple(owner_map(owner) for owner in request.regenerated_template_ownership)
            if request.regenerated_template_ownership is not None else None)
    return EngineContextView(internal, tuple((atom_map(atom), atom) for atom in represented))


def selected_coordinate_digest(plan: ProteinContextPlan) -> str:
    payload = [[item.atom.chain, item.atom.residue_number, item.atom.insertion_code, item.atom.atom_name, *item.coordinate_angstrom] for item in plan.fixed_conditioning]
    return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def context_evidence(request: DiffusionRequest, atom_axis: tuple[AtomIdentity, ...]) -> ProteinContextEvidence:
    plan = request.protein_context
    if plan is None:
        raise ProteinContextError("context evidence requires an explicit selection plan")
    selected = tuple(item.atom for item in plan.fixed_conditioning)
    target_chains = {gap.chain for gap in request.gaps}
    primary_residues = {
        residue for entity in sampling_entities(request)
        if set(entity.residues) & {_residue(atom) for atom in request.generated_atoms}
        for residue in entity.residues
    }
    alignment = tuple(atom for atom in selected if _residue(atom) in primary_residues)
    capability = "multichain-conditioned" if len(target_chains) > 1 or any(atom.chain not in target_chains for atom in selected) else "same-chain-conditioned"
    if len(target_chains) == 1 and plan.policy.target_flanks == 1 and plan.policy.partner_mode == "none" and not plan.policy.same_chain_spatial:
        capability = "gap-only"
    return ProteinContextEvidence(protein_context_digest(plan), atom_identity_digest(atom_axis), selected_coordinate_digest(plan), selected, plan.preserved_only, alignment, capability)


def validate_context_source(request: DiffusionRequest, pdb: bytes) -> None:
    if request.gap_contexts is not None:
        for single in independent_gap_requests(request):
            validate_context_source(single, pdb)
        return
    if request.protein_context is None:
        return
    points = source_coordinates(request, pdb)
    for item in request.protein_context.fixed_conditioning:
        if not np.array_equal(points[item.atom], item.coordinate_angstrom):
            raise ProteinContextError("declared conditioning coordinates do not match normalized source")
    # Recompute the complete deterministic selection, not merely its supplied hash.
    plan = request.protein_context
    expected = select_protein_context(replace(request, protein_context=None), pdb, plan.policy,
        trajectory_queries=plan.trajectory_queries, retained_conditioning=plan.retained_conditioning,
        previous_plan_sha256=plan.previous_plan_sha256, adaptation_source=plan.adaptation_source)
    if expected.protein_context != request.protein_context:
        raise ProteinContextError("protein context plan differs from deterministic source selection")


def validate_context_capability(request: DiffusionRequest, model: str) -> None:
    """Reject engine-incompatible plans before preflight or checkpoint loading."""
    if request.gap_contexts is not None:
        if not (is_mps_cc89(request) or is_cc94(request)) or model != "protpardelle":
            raise ProteinContextError("independent local gap contexts require cc94 or corrected MPS cc89")
        for single in independent_gap_requests(request):
            validate_context_capability(single, model)
        return
    if is_cuda_cc89(request):
        validate_cuda_cc89_request(request)
    if is_cc94(request):
        validate_cc94_request(request)
    if is_mps_cc89(request):
        validate_mps_cc89_request(request)
    if request.protein_context is None:
        if model == "protenix" and (len(request.target_sequences) != 1 or len(request.gaps) != 1):
            raise ProteinContextError("Protenix multiple chains or gaps require explicit --diffusion-context local|full")
        return
    if sum(span.stop - span.start for span in request.protein_context.spans) > 512:
        raise ProteinContextError("explicit protein context exceeds the audited 512-residue engine-axis budget")
    if model == "protenix":
        try:
            sampling_entities(request)
        except ValueError as exc:
            raise ProteinContextError(f"unsupported Protenix context mapping: {exc}") from exc
        return
    if model != "protpardelle":
        raise ProteinContextError("selected engine does not support explicit protein context")
    if len({gap.chain for gap in request.gaps}) != 1:
        raise ProteinContextError("Protpardelle explicit context does not support joint gap-bearing chains")
    spans = request.protein_context.spans
    if not is_corrected_cc89(request) and len({span.chain for span in spans}) != len(spans):
        raise ProteinContextError("Protpardelle explicit context does not support discontiguous source-chain entities")
    if len(spans) > 26:
        raise ProteinContextError("Protpardelle explicit context supports at most 26 engine chains")
    canonical_names = {name for names in CANONICAL_HEAVY_ATOMS.values() for name in names}
    if any(item.atom.atom_name not in canonical_names for item in request.protein_context.fixed_conditioning):
        raise ProteinContextError("Protpardelle explicit context contains atoms outside its audited canonical-heavy mask")


def select_legacy_partner_spans(
    request: DiffusionRequest,
    target_chain: str,
    coordinates: dict[AtomIdentity, np.ndarray],
    *,
    residue_budget: int,
) -> tuple[ProteinContextSpan, ...]:
    """Compatibility policy for the previously accepted cc89 local crops.

    This preserves the legacy 12 Å, one-run-per-chain, budget-trim behavior.
    Explicit context plans use the complete fail-closed selector above instead.
    Keeping the old policy here removes scientific selection from the adapter.
    """
    if residue_budget <= 0:
        return ()
    anchors = {anchor for gap in request.gaps if gap.chain == target_chain
               for anchor in (gap.left_anchor, gap.right_anchor) if anchor is not None}
    anchor_points = np.asarray([xyz for atom, xyz in coordinates.items() if _residue(atom) in anchors])
    if not len(anchor_points):
        return ()
    by_residue: dict[ResidueIdentity, list[np.ndarray]] = {}
    for atom, xyz in coordinates.items():
        by_residue.setdefault(_residue(atom), []).append(xyz)
    candidates = []
    for placement in request.sequence_placements:
        if placement.chain == target_chain:
            continue
        observed = tuple(zip(placement.observed_target_indices, placement.observed_residues))
        distances = {}
        for index, residue in observed:
            points = by_residue.get(residue, [])
            if points:
                distance = float(np.linalg.norm(np.asarray(points)[:, None] - anchor_points[None], axis=2).min())
                if distance <= 12.0:
                    distances[index] = distance
        if not distances:
            continue
        runs: list[list[tuple[int, ResidueIdentity]]] = []
        for item in observed:
            if not runs or item[0] != runs[-1][-1][0] + 1:
                runs.append([item])
            else:
                runs[-1].append(item)
        eligible = [run for run in runs if any(index in distances for index, _residue_id in run)]
        run = min(eligible, key=lambda items: (-sum(index in distances for index, _residue_id in items),
                                              min(distances.get(index, float("inf")) for index, _residue_id in items), items[0][0]))
        contacts = [index for index, _residue_id in run if index in distances]
        start = max(run[0][0], min(contacts) - 2)
        stop = min(run[-1][0] + 1, max(contacts) + 3)
        candidates.append((min(distances[index] for index in contacts), placement.chain, start, stop,
                           min(contacts, key=distances.__getitem__), run[0][0]))
    remaining = residue_budget
    spans = []
    for _distance, chain, start, stop, closest, run_start in sorted(candidates):
        if remaining <= 0:
            break
        if stop - start > remaining:
            start = max(start, closest - (remaining - 1) // 2)
            stop = min(stop, start + remaining)
            start = max(run_start, stop - remaining)
        spans.append(ProteinContextSpan(chain, start, stop))
        remaining -= stop - start
    return tuple(spans)
