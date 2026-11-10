"""Private, multi-parent gap selection; never a schema-4 sampling candidate.

Only generated-only donors in the exact same request/source frame are admitted.
Every whole-complex trial is independently assessed. Local cleanliness is only
a search heuristic and grants neither acceptance nor public publication rights.
"""
from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionContractError, DiffusionRequest
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment

PORTFOLIO_REVISION = "generated-only-gap-portfolio-v1"
PORTFOLIO_KIND = "experimental-multiparent-gap-portfolio"
MAX_PARENTS = 16


def digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class GapParent:
    """Declared lineage plus content-addressed coordinates and sampler custody.

    The caller must load the trace/evidence from verified private artifacts;
    parent labels and seeds do not replace their digests as identities.
    """

    candidate_id: str
    seed: int
    request_sha256: str
    source_sha256: str
    coordinate_sha256: str
    sampler_trace_sha256: str
    evidence_sha256: str
    pdb_text: str
    stage: str = "generated-only"

    def __post_init__(self):
        if not self.candidate_id or type(self.seed) is not int or self.seed < 0:
            raise DiffusionContractError("portfolio parent requires a candidate identity and seed")
        for value in (self.request_sha256, self.source_sha256, self.coordinate_sha256,
                      self.sampler_trace_sha256, self.evidence_sha256):
            if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise DiffusionContractError("portfolio parent requires canonical SHA-256 custody")
        if digest(self.pdb_text) != self.coordinate_sha256 or self.stage != "generated-only":
            raise DiffusionContractError("portfolio parent coordinates/stage differ from custody")


@dataclass(frozen=True)
class GapChoice:
    gap_index: int
    parent_index: int


@dataclass(frozen=True)
class GapRepairScope:
    gap_indices: tuple[int, ...]
    fixed_gap_indices: tuple[int, ...]
    unattributed_categories: tuple[str, ...]


@dataclass(frozen=True)
class PortfolioTrial:
    choices: tuple[GapChoice, ...]
    coordinate_sha256: str
    assessment: RepairAssessment

    @property
    def score(self):
        return (len(self.assessment.admission_failures),
                sum(issue.actionable for issue in self.assessment.issues))


@dataclass(frozen=True)
class PortfolioResult:
    coordinate_text: str
    selected: PortfolioTrial
    trials: tuple[PortfolioTrial, ...]
    parent_gap_issue_counts: tuple[tuple[int, ...], ...]
    repair_scope: GapRepairScope
    budget_exhausted: bool
    revision: str = PORTFOLIO_REVISION
    kind: str = PORTFOLIO_KIND
    publication_qualified: bool = False


@dataclass(frozen=True)
class PortfolioLocalRepair:
    coordinate_text: str
    scope: GapRepairScope
    mobile_atoms: tuple[AtomIdentity, ...]
    input_sha256: str
    output_sha256: str
    assessment: RepairAssessment
    refinement_seed: int
    max_iterations: int
    restart_count: int
    search_diagnostics: tuple[dict, ...]
    operation: str = "generated-only-issue-local-refinement-v1"
    publication_qualified: bool = False


def _gap_residues(request):
    return [set(gap.generated_residues) | {r for r in (gap.left_anchor, gap.right_anchor) if r is not None}
            for gap in request.gaps]


def plan_gap_repairs(request: DiffusionRequest, assessment: RepairAssessment) -> GapRepairScope:
    """Freeze clean gaps; unattributed findings conservatively require all gaps.

    Aggregate gate findings with detailed counterparts do not duplicate scope.
    Non-geometric admission failures never authorize coordinate repair.
    """
    counterparts = {
        "severe-steric-overlap": {"clash"},
        "generated-or-junction-bond-length": {"bond_length"},
        "generated-or-junction-bond-angle": {"bond_angle"},
        "generated-or-junction-amide-planarity": {"non_planar_amide"},
        "generated-or-junction-ramachandran": {"ramachandran_outlier"},
        "generated-or-junction-sidechain-chi12": {"sidechain_chi12", "sidechain_chi12_outlier"},
    }
    attributed = {issue.category for issue in assessment.issues if issue.actionable and issue.residues}
    unknown = set()
    selected = set()
    for issue in assessment.issues:
        if not issue.actionable:
            continue
        if issue.residues:
            selected.update(i for i, residues in enumerate(_gap_residues(request)) if residues & set(issue.residues))
        elif issue.category not in attributed and not counterparts.get(issue.category, set()) & attributed:
            unknown.add(issue.category)
    if unknown:
        selected.update(range(len(request.gaps)))
    return GapRepairScope(tuple(sorted(selected)),
        tuple(i for i in range(len(request.gaps)) if i not in selected), tuple(sorted(unknown)))


def _coordinates(text):
    import numpy as np

    result, records = {}, {}
    for line in text.splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in result or line[16:17].strip():
            raise DiffusionContractError("portfolio parent has duplicate identities or unselected altlocs")
        xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise DiffusionContractError("portfolio parent has nonfinite coordinates")
        result[atom], records[atom] = xyz, line
    return result, records


def search_gap_portfolio(
    request: DiffusionRequest,
    source_text: str,
    parents: tuple[GapParent, ...],
    *,
    max_evaluations: int = 64,
    beam_width: int = 4,
    retain_trial: Callable[[PortfolioTrial, str], None] | None = None,
) -> PortfolioResult:
    """Bounded deterministic beam over complete complexes, no sampling/refinement.

    The unselected suffix always comes from an explicit parent, never missing
    coordinates. All trial states have ONE donor per complete gap. Source records
    and generated atom metadata stay in the immutable base projection.
    """
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    if (not 1 <= len(parents) <= MAX_PARENTS or type(max_evaluations) is not int
            or not len(parents) <= max_evaluations <= 4096 or type(beam_width) is not int
            or not 1 <= beam_width <= 32):
        raise DiffusionContractError("invalid bounded portfolio search budget")
    if len({(p.candidate_id, p.seed) for p in parents}) != len(parents):
        raise DiffusionContractError("duplicate portfolio parent identity")
    if digest(source_text) != request.normalized_pdb.sha256:
        raise DiffusionContractError("portfolio source differs from request frame")
    if not request.gaps:
        raise DiffusionContractError("portfolio requires declared whole gaps")
    owned = [{atom for atom in request.generated_atoms
              if any((atom.chain, atom.residue_number, atom.insertion_code) ==
                     (r.chain, r.residue_number, r.insertion_code) for r in gap.generated_residues)}
             for gap in request.gaps]
    if set.union(*owned) != set(request.generated_atoms) or sum(map(len, owned)) != len(request.generated_atoms):
        raise DiffusionContractError("portfolio gap ownership is incomplete or overlapping")
    assess = assessment_adapter(request, source_text)
    axes, records, assessments = [], [], []
    for parent in parents:
        if (parent.request_sha256 != digest(request.to_json()) or parent.source_sha256 != digest(source_text)
                or parent.seed not in request.seeds):
            raise DiffusionContractError("portfolio parent request/source/seed differs")
        assessment = assess(parent.pdb_text)
        if assessment.admission_failures:
            raise DiffusionContractError("portfolio parent changed observations, identity or chemistry")
        structure = _structure(parent.pdb_text, "portfolio-parent.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        xyz, rows = _coordinates(parent.pdb_text)
        if not set(request.generated_atoms) <= xyz.keys():
            raise DiffusionContractError("portfolio parent omits generated atoms")
        if records:
            for atom in request.generated_atoms:
                # Serial numbers can differ; residue chemistry/occupancy cannot.
                a, b = records[0][atom], rows[atom]
                if a[11:30] != b[11:30] or a[54:] != b[54:]:
                    raise DiffusionContractError("portfolio parent generated metadata differs")
        axes.append(xyz)
        records.append(rows)
        assessments.append(assessment)
    counts = tuple(tuple(sum(issue.actionable and bool(set(issue.residues) & residues)
                             for issue in assessment.issues)
                         for residues in _gap_residues(request)) for assessment in assessments)
    trials, texts, seen = [], {}, {}

    def evaluate(indices):
        if indices in seen:
            return seen[indices]
        if len(trials) >= max_evaluations:
            return None
        coordinates = {atom: axes[parent_index][atom] for gap_index, parent_index in enumerate(indices)
                       for atom in owned[gap_index]}
        text = rewrite_generated_coordinates(parents[0].pdb_text, coordinates)
        trial = PortfolioTrial(tuple(GapChoice(i, p) for i, p in enumerate(indices)), digest(text), assess(text))
        trials.append(trial)
        texts[indices], seen[indices] = text, trial
        if retain_trial is not None:
            retain_trial(trial, text)
        return trial

    def key(indices):
        return (*seen[indices].score, indices)

    # Evaluate homogeneous parents as controls before any multi-parent choice.
    beam = []
    for p in range(len(parents)):
        indices = (p,) * len(request.gaps)
        evaluate(indices)
        beam.append(indices)
    beam = sorted(beam, key=key)[:beam_width]
    for gap_index in sorted(range(len(request.gaps)), key=lambda i: (-max(c[i] for c in counts), i)):
        if seen[beam[0]].score == (0, 0):
            break
        expanded = set(beam)
        for current in beam:
            for parent_index in sorted(range(len(parents)), key=lambda p: (counts[p][gap_index], p)):
                indices = current[:gap_index] + (parent_index,) + current[gap_index + 1:]
                if evaluate(indices) is not None:
                    expanded.add(indices)
        beam = sorted(expanded, key=key)[:beam_width]
        if len(trials) >= max_evaluations:
            break
    selected = min(seen, key=key)
    trial = seen[selected]
    return PortfolioResult(texts[selected], trial, tuple(trials), counts,
        plan_gap_repairs(request, trial.assessment), len(trials) >= max_evaluations)


def repair_portfolio_gaps(
    request: DiffusionRequest,
    source_text: str,
    coordinate_text: str,
    *,
    refinement_seed: int,
    max_iterations: int = 250,
    restart_count: int = 4,
    platform_name: str = "CPU",
) -> PortfolioLocalRepair:
    """One bounded repair; all clean generated gaps and observations stay fixed.

    The refinement seed is an optimizer parameter, never a new sampling seed.
    Any residual/new whole-complex defect remains an independently failed result.
    """
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    if (digest(source_text) != request.normalized_pdb.sha256 or type(refinement_seed) is not int
            or refinement_seed < 0 or type(max_iterations) is not int or not 1 <= max_iterations <= 2000
            or type(restart_count) is not int or not 0 <= restart_count <= 8):
        raise DiffusionContractError("invalid issue-local repair frame or budget")
    assess = assessment_adapter(request, source_text)
    before = assess(coordinate_text)
    if before.admission_failures:
        raise DiffusionContractError("issue-local repair cannot authorize admission failures")
    scope = plan_gap_repairs(request, before)
    residues = tuple(dict.fromkeys(r for i in scope.gap_indices for r in request.gaps[i].generated_residues))
    identities = {(r.chain, r.residue_number, r.insertion_code) for r in residues}
    atoms = tuple(a for a in request.generated_atoms if (a.chain, a.residue_number, a.insertion_code) in identities)
    text, diagnostics = coordinate_text, ()
    if atoms:
        def score(xyz):
            assessment = assess(rewrite_generated_coordinates(coordinate_text, xyz))
            return (len(assessment.admission_failures), sum(i.actionable for i in assessment.issues))

        result = refine_generated_region(coordinate_text, generated_residues=residues, generated_atoms=atoms,
            max_iterations=max_iterations, restart_count=restart_count, random_seed=refinement_seed,
            perturbation_angstrom=0.25, platform_name=platform_name, geometry_preconditioning=True,
            candidate_score=score)
        text = rewrite_generated_coordinates(coordinate_text, result.coordinates_angstrom)
        diagnostics = result.search_diagnostics
    return PortfolioLocalRepair(text, scope, atoms, digest(coordinate_text), digest(text), assess(text),
        refinement_seed, max_iterations, restart_count, diagnostics)


def audit_portfolio_materialization(
    request: DiffusionRequest,
    source_text: str,
    parents: tuple[GapParent, ...],
    choices: tuple[GapChoice, ...],
    candidate_text: str,
    *,
    repaired_text: str | None = None,
) -> RepairAssessment:
    """Independent identity/lineage/mobility/geometry audit, no runner authority.

    Optional issue-local output is authorized only by recomputed candidate
    findings. Neither a claimed clean mask nor a claimed selected score is read.
    """
    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    if digest(source_text) != request.normalized_pdb.sha256:
        raise DiffusionContractError("portfolio audit source digest mismatch")
    if tuple(c.gap_index for c in choices) != tuple(range(len(request.gaps))):
        raise DiffusionContractError("portfolio lineage must select each complete gap exactly once")
    if not parents or any(type(c.parent_index) is not int or not 0 <= c.parent_index < len(parents) for c in choices):
        raise DiffusionContractError("portfolio lineage references an absent parent")
    assess = assessment_adapter(request, source_text)
    axes = []
    for parent in parents:
        if (parent.request_sha256 != digest(request.to_json()) or parent.source_sha256 != digest(source_text)
                or parent.seed not in request.seeds or digest(parent.pdb_text) != parent.coordinate_sha256):
            raise DiffusionContractError("portfolio audit parent binding differs")
        if assess(parent.pdb_text).admission_failures:
            raise DiffusionContractError("portfolio audit parent changed source records")
        axes.append(_coordinates(parent.pdb_text)[0])
    actual, rows = _coordinates(candidate_text)
    _, base_rows = _coordinates(parents[0].pdb_text)
    expected_coordinates = {}
    for choice in choices:
        residues = {(r.chain, r.residue_number, r.insertion_code) for r in request.gaps[choice.gap_index].generated_residues}
        for atom in request.generated_atoms:
            if (atom.chain, atom.residue_number, atom.insertion_code) not in residues:
                continue
            expected = axes[choice.parent_index].get(atom)
            if atom not in actual or expected is None or tuple(actual[atom]) != tuple(expected):
                raise DiffusionContractError("portfolio output differs from its selected parent coordinates")
            if rows[atom][:30] != base_rows[atom][:30] or rows[atom][54:] != base_rows[atom][54:]:
                raise DiffusionContractError("portfolio output changed generated atom metadata")
            expected_coordinates[atom] = expected
    if candidate_text != rewrite_generated_coordinates(parents[0].pdb_text, expected_coordinates):
        raise DiffusionContractError("portfolio output changed non-coordinate records")
    before = assess(candidate_text)
    text = candidate_text
    if repaired_text is not None:
        scope = plan_gap_repairs(request, before)
        mobile = {r for i in scope.gap_indices for r in request.gaps[i].generated_residues}
        _, repaired_rows = _coordinates(repaired_text)
        if set(repaired_rows) != set(rows):
            raise DiffusionContractError("issue-local output changed the atom inventory")
        for atom, line in rows.items():
            new_line = repaired_rows[atom]
            identity = (atom.chain, atom.residue_number, atom.insertion_code)
            allowed = any(identity == (r.chain, r.residue_number, r.insertion_code) for r in mobile)
            if not allowed and line != new_line:
                raise DiffusionContractError("issue-local output moved a clean gap or observation")
            if allowed and (line[:30] != new_line[:30] or line[54:] != new_line[54:]):
                raise DiffusionContractError("issue-local output changed atom metadata")
        repaired_coordinates, _ = _coordinates(repaired_text)
        allowed_coordinates = {a: xyz for a, xyz in repaired_coordinates.items()
            if (a.chain, a.residue_number, a.insertion_code) in
               {(r.chain, r.residue_number, r.insertion_code) for r in mobile}}
        if repaired_text != rewrite_generated_coordinates(candidate_text, allowed_coordinates):
            raise DiffusionContractError("issue-local output changed non-coordinate records")
        text = repaired_text
    structure = _structure(text, "audited-portfolio.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    return assess(text)
