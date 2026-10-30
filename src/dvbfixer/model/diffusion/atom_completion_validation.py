"""Independent identity and geometry audits for private completion candidates.

Passing these checks is not runner/checkpoint qualification or publication.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from io import StringIO

from openmm.app import PDBFile

from dvbfixer.diagnose.chemistry import (
    check_backbone_bond_angles,
    check_bond_lengths,
    check_ca_chirality,
    check_peptide_omegas,
)
from dvbfixer.diagnose.report import Finding, Severity
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.diagnose.structural import check_coincident_heavy_atoms
from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.atom_completion_contract import (
    AtomCompletionRequest,
    audit_completion_request_source,
)
from dvbfixer.model.diffusion.atom_completion_materialize import _source_layout
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    Metric,
    ResidueIdentity,
    ValidationSummary,
)
from dvbfixer.model.diffusion.pdb_materialize import _coordinate_identity
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE
from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12
from dvbfixer.model.diffusion.scope import _parse_coordinate_records

COMPLETION_VALIDATION_REVISION = "atom-completion-validation-v1"


@dataclass(frozen=True, slots=True)
class CompletionCandidateAudit:
    source_sha256: str
    candidate_sha256: str
    generated_coordinates: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True, slots=True)
class CompletionValidation:
    summary: ValidationSummary
    source_findings: tuple[Finding, ...]
    candidate_findings: tuple[Finding, ...]
    source_sha256: str
    candidate_sha256: str
    revision: str = COMPLETION_VALIDATION_REVISION


def audit_completion_candidate(
    source: bytes, candidate: bytes, request: AtomCompletionRequest,
) -> CompletionCandidateAudit:
    """Audit raw bytes independently of the materializer and runner's claims.

    Original atoms (including hydrogen observations), ANISOU and explicit link
    records must be verbatim and ordered. Only the exact requested heavy atoms
    may be added. Generated names, elements and coordinates are read from PDB.
    """
    audit_completion_request_source(request, source)
    source_lines, _, _ = _source_layout(source.decode("utf-8"))
    lines, _, _ = _source_layout(candidate.decode("utf-8"))
    original, _, _ = _parse_coordinate_records(source_lines)
    atoms, names, altlocs = _parse_coordinate_records(lines)
    generated = set(request.generated_atoms)
    if altlocs:
        raise DiffusionContractError("completion candidate introduces alternate locations")
    expected = {atom.identity for atom in original} | generated
    if {atom.identity for atom in atoms} != expected:
        raise DiffusionContractError("completion candidate has missing or unrequested atom identities")
    protected_records = {"ATOM", "HETATM", "ANISOU"}
    source_observations = [line for line in source_lines if line[:6].strip() in protected_records]
    retained: list[str] = []
    for line in lines:
        record = line[:6].strip()
        if record == "ANISOU":
            retained.append(line)
        elif record in {"ATOM", "HETATM"}:
            identity = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
            if identity not in generated:
                retained.append(line)
    if retained != source_observations:
        raise DiffusionContractError("completion candidate changed original observation records")
    links = {"CONECT", "LINK", "SSBOND"}
    if ([line for line in lines if line[:6].strip() in links]
            != [line for line in source_lines if line[:6].strip() in links]):
        raise DiffusionContractError("completion candidate changed explicit source links")
    expected_names = {partial.residue: partial.residue_name for partial in request.masks.partial_residues}
    targets = {target.chain: target.sequence for target in request.target_sequences}
    for gap in request.gaps:
        expected_names.update((residue, _CANONICAL_ONE_TO_THREE[targets[gap.chain][index]])
                              for index, residue in zip(range(gap.target_interval.start, gap.target_interval.stop),
                                                        gap.generated_residues, strict=True))
    by_identity = {atom.identity: atom for atom in atoms}
    for identity in request.generated_atoms:
        atom = by_identity[identity]
        residue = ResidueIdentity(identity.chain, identity.residue_number, identity.insertion_code)
        if (atom.record_name != "ATOM" or names[residue] != expected_names[residue]
                or atom.element != identity.atom_name[0]):
            raise DiffusionContractError("completion generated atom chemistry differs from canonical request")
    coordinates = _coordinates(lines)
    return CompletionCandidateAudit(hashlib.sha256(source).hexdigest(), hashlib.sha256(candidate).hexdigest(),
                                    tuple(coordinates[atom] for atom in request.generated_atoms))


def validate_completion_candidate(
    source: bytes, candidate: bytes, request: AtomCompletionRequest,
) -> CompletionValidation:
    """Retain all source/candidate findings; reject local severe defects and any D-Cα.

    Identity/record failures raise rather than producing a geometry-only pass.
    Existing defects are never waived in the changed residue/junction cohort.
    Outside findings remain visible, not silently erased from the report.
    """
    audit = audit_completion_candidate(source, candidate, request)
    original = PDBFile(StringIO(source.decode("utf-8")))
    completed = PDBFile(StringIO(candidate.decode("utf-8")))
    missing_findings = tuple(Finding(
        Severity.ERROR, "missing_atoms", partial.residue.chain,
        partial.residue.residue_number + partial.residue.insertion_code, partial.residue_name,
        "source is missing canonical heavy atoms: " + ", ".join(atom.atom_name for atom in partial.missing_atoms),
    ) for partial in request.masks.partial_residues)
    gap_findings = tuple(Finding(
        Severity.ERROR, "missing_residues", gap.chain, gap.generated_residues[0].residue_number,
        "", f"source is missing {len(gap.generated_residues)} whole residues at target positions "
        f"{gap.target_interval.start + 1}-{gap.target_interval.stop}",
    ) for gap in request.gaps)
    source_findings = missing_findings + gap_findings + _geometry_findings(original)
    affected = _affected_residue_locations(request)
    peptide_findings = _peptide_findings(candidate, request, affected)
    findings = _geometry_findings(completed) + peptide_findings
    local_errors = [finding for finding in findings
                    if (finding.severity is Severity.ERROR or finding.category == "non_planar_amide")
                    and _touches(finding, affected)]
    failures = list(dict.fromkeys(f"completion-{finding.category}" for finding in local_errors))
    if any(finding.category == "coincident_atoms" and "truncated" in finding.message.lower() for finding in findings):
        failures.append("completion-diagnostics-incomplete")
    try:
        assert_all_l(completed.topology, completed.positions)
        d_count = 0
    except ChiralityError as exc:
        d_count = len(exc.residues)
        failures.append("completion-d-ca-chirality")
    metrics = (
        Metric("completion-generated-atoms", float(len(request.generated_atoms)), "count"),
        Metric("completion-source-missing-atoms", float(len(request.masks.missing_atoms)), "count"),
        Metric("completion-fixed-record-displacement", 0.0, "angstrom"),
        Metric("completion-source-findings", float(len(source_findings)), "count"),
        Metric("completion-candidate-findings", float(len(findings)), "count"),
        Metric("completion-local-errors", float(len(local_errors)), "count"),
        Metric("completion-detectable-d-ca", float(d_count), "count"),
    )
    warnings = tuple(f"source {finding.severity.value}: {finding.chain}/{finding.resid} "
                     f"{finding.category}: {finding.message}" for finding in source_findings)
    return CompletionValidation(ValidationSummary(not failures, tuple(dict.fromkeys(failures)), metrics, warnings),
                                source_findings, findings, audit.source_sha256, audit.candidate_sha256)


def _geometry_findings(pdb: PDBFile) -> tuple[Finding, ...]:
    # Shared diagnose defaults, not relaxed thresholds selected for a candidate.
    return tuple([
        *check_bond_lengths(pdb.topology, pdb.positions),
        *check_backbone_bond_angles(pdb.topology, pdb.positions),
        *check_peptide_omegas(pdb.topology, pdb.positions),
        *check_ca_chirality(pdb.topology, pdb.positions),
        *check_ramachandran(pdb.topology, pdb.positions),
        *check_sidechain_chi12(pdb.topology, pdb.positions),
        *check_coincident_heavy_atoms(pdb.topology, pdb.positions),
        *clashes_python(pdb.topology, pdb.positions, clash_warn_a=0.6, clash_error_a=0.9),
    ])


def _affected_residue_locations(request: AtomCompletionRequest) -> set[tuple[str, str]]:
    residues = {partial.residue for partial in request.masks.partial_residues}
    for gap in request.gaps:
        residues.update(gap.generated_residues)
        residues.update(anchor for anchor in (gap.left_anchor, gap.right_anchor) if anchor is not None)
    # Adjacent target positions cover geometry reported on either side of a
    # generated backbone bond, even if diagnose labels its observed neighbour.
    positions = _target_positions(request)
    selected = {position for position, residue in positions.items() if residue in residues}
    for chain, index in selected:
        for neighbour in (index - 1, index + 1):
            residue = positions.get((chain, neighbour))
            if residue is not None:
                residues.add(residue)
    return {(residue.chain, residue.residue_number + residue.insertion_code) for residue in residues}


def _target_positions(request: AtomCompletionRequest) -> dict[tuple[str, int], ResidueIdentity]:
    positions: dict[tuple[str, int], ResidueIdentity] = {}
    for placement in request.sequence_placements:
        positions.update(((placement.chain, index), residue) for index, residue in zip(
            placement.observed_target_indices, placement.observed_residues, strict=True))
    for gap in request.gaps:
        positions.update(((gap.chain, index), residue) for index, residue in zip(
            range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
    return positions


def _touches(finding: Finding, locations: set[tuple[str, str]]) -> bool:
    if (finding.chain, finding.resid) in locations:
        return True
    for atom in finding.extra.get("atoms", ()):
        if (atom.get("chain"), atom.get("resid", "") + atom.get("icode", "")) in locations:
            return True
    partner = finding.extra.get("clash_partner", {})
    return (partner.get("chain"), partner.get("resid")) in locations


def _coordinates(lines: list[str]) -> dict[AtomIdentity, tuple[float, float, float]]:
    return {identity: (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            for line in lines if (identity := _coordinate_identity(line)) is not None}


def _peptide_findings(candidate: bytes, request: AtomCompletionRequest,
                      affected: set[tuple[str, str]]) -> tuple[Finding, ...]:
    coordinates = _coordinates(candidate.decode("utf-8").splitlines())
    positions = _target_positions(request)
    findings = []
    for target in request.target_sequences:
        for index in range(len(target.sequence) - 1):
            left, right = positions[target.chain, index], positions[target.chain, index + 1]
            if not any((r.chain, r.residue_number + r.insertion_code) in affected for r in (left, right)):
                continue
            c = coordinates[AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")]
            n = coordinates[AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")]
            distance = math.dist(c, n)
            if not 1.20 <= distance <= 1.45:
                findings.append(Finding(
                    Severity.ERROR, "peptide_connectivity", left.chain,
                    left.residue_number + left.insertion_code, "",
                    f"target-neighbour C-N distance {distance:.3f} Å is outside 1.20-1.45 Å; "
                    f"right residue {right.chain}/{right.residue_number}{right.insertion_code}",
                    extra={"observed_angstrom": distance, "minimum_angstrom": 1.20, "maximum_angstrom": 1.45},
                ))
    return tuple(findings)
