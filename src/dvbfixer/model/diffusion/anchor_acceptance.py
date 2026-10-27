"""Independent, research-only audit of staged observed-anchor repair.

No existing runner contract or publication gate is weakened. Movement limits are
explicit caller policy, never runner scores, B-derived radii or implicit defaults.
The auditor determines eligible anchors itself from stage-1 geometry findings.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import numpy as np
from openmm.app import PDBFile

from dvbfixer.diagnose.chemistry import (
    check_backbone_bond_angles,
    check_bond_lengths,
    check_peptide_omegas,
)
from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovement,
    AnchorMovementLimits,
    atom_residue,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    DiffusionRequest,
    Metric,
    ResidueIdentity,
    ValidationSummary,
)
from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12
from dvbfixer.model.diffusion.refinement_stages import (
    NumericalFailureEvidence,
    RepairAssessment,
    RepairIssue,
    problematic_anchors,
)
from dvbfixer.model.diffusion.validate import (
    ValidationThresholds,
    _parse_pdb_text,
    _Structure,
    _validate_candidate,
)

ANCHOR_ACCEPTANCE_REVISION = "independent-staged-anchor-audit-v2"
_DRIFT_GATES = frozenset({"fixed-heavy-atom-rmsd", "fixed-heavy-atom-max-displacement"})
_GEOMETRY_GATES = frozenset({"junction-peptide-connectivity", "generated-backbone-break", "d-ca-chirality",
                            "severe-steric-overlap", "generated-or-junction-bond-length",
                            "generated-or-junction-bond-angle", "generated-or-junction-amide-planarity",
                            "generated-or-junction-ramachandran", "generated-or-junction-sidechain-chi12"})


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _structure(text: str, name: str) -> _Structure:
    atoms, residues, links = _parse_pdb_text(text)
    return _Structure(Path(name), PDBFile(StringIO(text)), atoms, residues, frozenset(links), text)


def _observed_lines(text: str) -> dict[AtomIdentity, str]:
    rows = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in rows:
            raise DiffusionContractError("duplicate atom identity in staged anchor evidence")
        rows[atom] = line
    return rows


def _observed_auxiliary(text: str, observed: set[AtomIdentity]) -> tuple[str, ...]:
    """Preserve experimental atom-associated metadata even for mobile anchors."""
    lines = []
    for line in text.splitlines():
        if line[:6].strip() not in {"ANISOU", "SIGATM", "SIGUIJ"}:
            continue
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if atom in observed:
            lines.append(line)
    return tuple(sorted(lines))


def _neighborhood(request: DiffusionRequest) -> set[ResidueIdentity]:
    """Same generated+all-direct-anchor+outer-neighbor neighborhood at every stage."""
    active = {atom_residue(atom) for atom in request.generated_atoms}
    active.update(r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None)
    result = set(active)
    for placement in request.sequence_placements:
        positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for gap in request.gaps:
            if gap.chain == placement.chain:
                positions.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
        indices = {i for i, r in positions.items() if r in active}
        expanded = {j for i in indices for j in (i - 1, i, i + 1)}
        result.update(r for i, r in positions.items() if i in expanded)
    return result


def _affected(finding: Finding, known: set[ResidueIdentity]) -> tuple[ResidueIdentity, ...]:
    affected = set()
    for atom in finding.extra.get("atoms", []):
        affected.add(ResidueIdentity(atom["chain"], str(atom["resid"]), atom.get("icode", "")))
    # Older diagnostics and clash partners use display IDs. Match exact known
    # identities, never strip insertion codes or casefold chain IDs.
    partner = finding.extra.get("clash_partner", {})
    for chain, display in ((finding.chain, finding.resid), (partner.get("chain"), partner.get("resid"))):
        affected.update(r for r in known if r.chain == chain and r.residue_number + r.insertion_code == display)
    return tuple(sorted(affected))


def _findings(structure: _Structure, request: DiffusionRequest) -> tuple[RepairIssue, ...]:
    local = _neighborhood(request)
    known = set(structure.residues) | local
    issues = []
    findings = []
    for check in (check_bond_lengths, check_backbone_bond_angles, check_peptide_omegas,
                  check_ramachandran, check_sidechain_chi12):
        findings.extend(check(structure.pdb.topology, structure.pdb.positions))
    findings.extend(clashes_python(structure.pdb.topology, structure.pdb.positions,
                                    clash_warn_a=0.9, clash_error_a=0.9))
    for finding in findings:
        if finding.severity is not Severity.ERROR and finding.category != "non_planar_amide":
            continue
        affected = _affected(finding, known)
        if set(affected) & local:
            issues.append(RepairIssue(finding.category, affected, finding.message))
    # Explicit request junctions also need distance checks: OpenMM can omit the
    # very C-N bond whose absence caused the defect, so topology-only checks miss it.
    thresholds = ValidationThresholds()
    anchors = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    generated = {atom_residue(atom) for atom in request.generated_atoms}
    for placement in request.sequence_placements:
        positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
        for gap in request.gaps:
            if gap.chain == placement.chain:
                positions.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
        for index, left in positions.items():
            right = positions.get(index + 1)
            if right is None or not {left, right} & (anchors | generated):
                continue
            c = structure.atoms.get(AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C"))
            n = structure.atoms.get(AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N"))
            if c is not None and n is not None:
                distance = float(np.linalg.norm(c.coordinates_angstrom - n.coordinates_angstrom))
                if {left, right} & anchors and not thresholds.junction_cn_angstrom_min <= distance <= thresholds.junction_cn_angstrom_max:
                    issues.append(RepairIssue("junction-peptide-connectivity", (left, right), f"junction C-N {distance:.4f} Å"))
                elif not {left, right} & anchors and distance > thresholds.generated_backbone_break_angstrom_max:
                    issues.append(RepairIssue("generated-backbone-break", (left, right), f"generated C-N {distance:.4f} Å"))
    return tuple(issues)


def _assessment(request: DiffusionRequest, source: _Structure, candidate: _Structure) -> tuple[RepairAssessment, ValidationSummary]:
    strict, _ = _validate_candidate(request, source, candidate, ValidationThresholds())
    issues = list(_findings(candidate, request))
    # Preserve an independently failing gate even if an older detector provides
    # no detailed affected-atom evidence. Such a gate cannot guess an anchor.
    issues.extend(RepairIssue(gate, (), gate) for gate in strict.hard_gate_failures if gate in _GEOMETRY_GATES)
    admission = tuple(gate for gate in strict.hard_gate_failures if gate not in _GEOMETRY_GATES)
    return RepairAssessment(tuple(issues), admission), strict


@dataclass(frozen=True)
class AnchorRepairAudit:
    summary: ValidationSummary
    historical_strict_summary: ValidationSummary
    authorized_anchors: tuple[ResidueIdentity, ...]
    movements: tuple[AnchorMovement, ...]
    movement_limits: AnchorMovementLimits
    stage1_failure: NumericalFailureEvidence | None
    source_sha256: str
    raw_sha256: str
    stage1_sha256: str | None
    candidate_sha256: str
    request_sha256: str
    policy_revision: str = ANCHOR_ACCEPTANCE_REVISION
    publication_qualified: bool = False


def audit_staged_anchor_repair(
    request: DiffusionRequest,
    *,
    source_text: str,
    raw_text: str,
    stage1_text: str | None,
    candidate_text: str,
    stage1_input_sha256: str,
    stage2_input_sha256: str,
    limits: AnchorMovementLimits,
    stage1_failure: NumericalFailureEvidence | None = None,
) -> AnchorRepairAudit:
    """Recompute authorization and full geometry, not runner-provided scores/masks.

    Caller-specified limits are experimental policy for this audit, not shipping
    defaults. A passed summary remains NOT qualified for public publication.
    All coordinates/metadata outside the authorized mobility mask stay exact.
    """
    source_digest, raw_digest = _digest(source_text), _digest(raw_text)
    if source_digest != request.normalized_pdb.sha256:
        raise DiffusionContractError("staged anchor source digest mismatch")
    if stage1_text is not None and stage1_failure is not None:
        raise DiffusionContractError("stage 1 cannot supply both coordinates and numerical failure evidence")
    source = _structure(source_text, "source.pdb")
    raw = _structure(raw_text, "raw.pdb")
    stage1 = _structure(stage1_text, "stage1.pdb") if stage1_text is not None else raw
    candidate = _structure(candidate_text, "candidate.pdb")
    raw_assessment, _ = _assessment(request, source, raw)
    first_assessment, _ = _assessment(request, source, stage1)
    final_assessment, strict = _assessment(request, source, candidate)
    failures = []
    if stage1_text is None:
        if stage1_failure is None:
            failures.append("anchor-repair-stage1-failure-evidence-missing")
        elif stage1_failure.input_sha256 != raw_digest:
            failures.append("anchor-repair-stage1-failure-input-mismatch")
    if stage1_input_sha256 != raw_digest or stage2_input_sha256 != raw_digest:
        failures.append("anchor-repair-not-restarted-from-raw")
    if raw_assessment.admission_failures:
        failures.append("anchor-repair-raw-admission")
    if first_assessment.admission_failures:
        failures.append("anchor-repair-stage1-admission")
    if raw_assessment.geometry_passed:
        failures.append("anchor-repair-unnecessary-raw-clean")
    if first_assessment.geometry_passed:
        failures.append("anchor-repair-unnecessary-stage1-clean")
    anchors = problematic_anchors(request, first_assessment)
    if not anchors:
        failures.append("anchor-repair-no-problematic-anchors")
    # Existing source-only problems do not authorize observed motion. Keep their
    # identities separate from junction defects introduced by reconstruction.
    baseline = {(issue.category, issue.residues) for issue in _findings(source, request)}
    eligible = tuple(issue for issue in first_assessment.issues
                     if (issue.category, issue.residues) not in baseline)
    justified = problematic_anchors(request, RepairAssessment(eligible))
    if set(anchors) - set(justified):
        failures.append("anchor-repair-preexisting-source-defect")
    anchors = justified
    restraints = build_anchor_restraints(source_text, anchors, observed_atoms=request.fixed_atoms,
                                         bfactor_source="disabled")
    allowed = {r.atom for r in restraints}
    source_lines, raw_lines = _observed_lines(source_text), _observed_lines(raw_text)
    first_text = stage1_text if stage1_text is not None else raw_text
    first_lines, final_lines = _observed_lines(first_text), _observed_lines(candidate_text)
    observed = source.atoms.keys() - set(request.generated_atoms)
    original_auxiliary = _observed_auxiliary(source_text, observed)
    if original_auxiliary != _observed_auxiliary(raw_text, observed) or original_auxiliary != _observed_auxiliary(first_text, observed):
        failures.append("anchor-repair-input-observed-auxiliary-changed")
    if original_auxiliary != _observed_auxiliary(candidate_text, observed):
        failures.append("anchor-repair-observed-auxiliary-changed")
    for atom in observed:
        original = source_lines.get(atom)
        raw_line, first_line, final_line = raw_lines.get(atom), first_lines.get(atom), final_lines.get(atom)
        if original is None or raw_line != original or first_line != original:
            failures.append("anchor-repair-input-observed-record-changed")
        if final_line is None or original is None:
            failures.append("anchor-repair-observed-record-missing")
        elif atom not in allowed and final_line != original:
            failures.append("anchor-repair-unselected-observed-record-changed")
        elif atom in allowed and (final_line[:30] != original[:30] or final_line[54:] != original[54:]):
            failures.append("anchor-repair-observed-metadata-changed")
    try:
        movements = check_anchor_movements(restraints,
            {atom: tuple(record.coordinates_angstrom) for atom, record in candidate.atoms.items()}, limits)
        if any(not movement.within_limits for movement in movements):
            failures.append("anchor-repair-per-anchor-movement-limit")
    except ValueError:
        movements = ()
        failures.append("anchor-repair-movement-measurement-invalid")
    failures.extend(gate for gate in strict.hard_gate_failures if gate not in _DRIFT_GATES)
    if not final_assessment.geometry_passed:
        failures.append("anchor-repair-expanded-geometry")
    metrics = list(strict.metrics)
    warnings = [*strict.warnings, "Research anchor audit only; not qualified for production publication."]
    for movement in movements:
        r = movement.residue
        label = f"{r.chain}:{r.residue_number}{r.insertion_code}"
        metrics.extend((Metric(f"anchor-{label}-backbone-rmsd", movement.backbone_rmsd_angstrom, "angstrom"),
                        Metric(f"anchor-{label}-max-heavy-displacement", movement.heavy_displacement_angstrom_max, "angstrom")))
        if movement.heavy_displacement_angstrom_max > 0:
            warnings.append(f"Observed anchor {label} coordinates modified: backbone RMSD "
                            f"{movement.backbone_rmsd_angstrom:.4f} Å; maximum heavy displacement "
                            f"{movement.heavy_displacement_angstrom_max:.4f} Å.")
    summary = ValidationSummary(not failures, tuple(dict.fromkeys(failures)), tuple(metrics), tuple(warnings))
    return AnchorRepairAudit(summary, strict, anchors, movements, limits, stage1_failure, source_digest, raw_digest,
                             _digest(stage1_text) if stage1_text is not None else None, _digest(candidate_text),
                             hashlib.sha256(request.to_json().encode("utf-8")).hexdigest())
