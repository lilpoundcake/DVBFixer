"""Independent completion audits and real generated-only OpenMM controls.

Reference coordinates below are a control, never diffusion-generated evidence.
"""

import hashlib
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.atom_completion import plan_atom_completion
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
from dvbfixer.model.diffusion.atom_completion_validation import (
    _affected_residue_locations,
    _touches,
    audit_completion_candidate,
    validate_completion_candidate,
)
from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    AtomIdentity,
    DiffusionContractError,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
)
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE
from dvbfixer.model.diffusion.scope import _parse_coordinate_records


@pytest.fixture
def completion_control(fixtures_root):
    original = b"".join(line for line in (fixtures_root / "protein_context/1UBQ.pdb").read_bytes().splitlines(keepends=True)
                        if line.startswith(b"ATOM  ")) + b"END\n"
    missing = {"CG", "CD", "CE", "NZ"}
    source = b"".join(line for line in original.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and line[22:26].strip() == b"27"
                              and line[12:16].strip().decode() in missing))
    records, names, _ = _parse_coordinate_records(source.decode().splitlines())
    residues = tuple(dict.fromkeys(ResidueIdentity(atom.identity.chain, atom.identity.residue_number,
                                                   atom.identity.insertion_code) for atom in records))
    reverse = {three: one for one, three in _CANONICAL_ONE_TO_THREE.items()}
    targets = (TargetSequence("A", "".join(reverse[names[residue]] for residue in residues)),)
    placements = (SequencePlacement("A", len(residues), tuple(range(len(residues))), residues),)
    digest = hashlib.sha256(source).hexdigest()
    masks = plan_atom_completion(source, source_sha256=digest, target_sequences=targets,
                                 sequence_placements=placements, fixed_atoms=tuple(atom.identity for atom in records))
    request = AtomCompletionRequest(1, ArtifactReference("source.pdb", digest), targets, placements,
                                    (), masks, (7,), "private-control-not-sampling")
    coordinates = {AtomIdentity(line[21:22].decode(), line[22:26].strip().decode(), line[26:27].strip().decode(),
                                line[12:16].strip().decode()): tuple(float(line[i:i + 8]) for i in (30, 38, 46))
                   for line in original.splitlines() if line.startswith(b"ATOM  ")}
    assert len(request.generated_atoms) == 4
    candidate = materialize_atom_completion(source, request, request.generated_atoms,
                                            [coordinates[atom] for atom in request.generated_atoms])
    return source, candidate, request


def test_independent_audit_and_real_reference_geometry(completion_control):
    source, candidate, request = completion_control
    audit = audit_completion_candidate(source, candidate, request)
    assert audit.source_sha256 == request.normalized_pdb.sha256
    assert audit.candidate_sha256 == hashlib.sha256(candidate).hexdigest()
    assert len(audit.generated_coordinates) == 4
    report = validate_completion_candidate(source, candidate, request)
    assert report.summary.passed, report.summary.hard_gate_failures
    assert report.source_findings  # Existing source issues are not hidden.
    assert len(report.summary.warnings) == len(report.source_findings)
    assert any(metric.name == "completion-fixed-record-displacement" and metric.value == 0
               for metric in report.summary.metrics)


@pytest.mark.parametrize("mutation,match", [
    (lambda text: text.replace(b" CG ", b" XX ", 1), "atom identities"),
    (lambda text: text.replace(b"ATOM      1", b"ATOM  99999", 1), "observation records"),
    (lambda text: text.replace(b"  27      ", b"  27A     ", 1), "atom identities"),
    (lambda text: text.replace(b"END\n", b"CONECT    1    2\nEND\n"), "explicit source links"),
])
def test_forged_candidates_refused_independently(completion_control, mutation, match):
    source, candidate, request = completion_control
    changed = mutation(candidate)
    assert changed != candidate
    with pytest.raises((DiffusionContractError, ValueError), match=match):
        audit_completion_candidate(source, changed, request)


def test_fixed_coordinate_change_refused_without_tolerance(completion_control):
    source, candidate, request = completion_control
    lines = candidate.splitlines(keepends=True)
    lines[0] = lines[0][:30] + f"{float(lines[0][30:38]) + .001:8.3f}".encode() + lines[0][38:]
    with pytest.raises(DiffusionContractError, match="observation records"):
        validate_completion_candidate(source, b"".join(lines), request)


def test_wrong_generated_chemistry_and_extra_atoms_refused(completion_control):
    source, candidate, request = completion_control
    lines = candidate.splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line[12:16].strip() == b"NZ" and line[22:26].strip() == b"27")
    wrong_element = lines.copy()
    wrong_element[index] = lines[index][:76] + b" C" + lines[index][78:]
    with pytest.raises(DiffusionContractError, match="chemistry"):
        audit_completion_candidate(source, b"".join(wrong_element), request)
    wrong_record = lines.copy()
    wrong_record[index] = b"HETATM" + lines[index][6:]
    with pytest.raises(DiffusionContractError, match="chemistry"):
        audit_completion_candidate(source, b"".join(wrong_record), request)
    extra = lines[index][:6] + b"99999" + b"  XX " + lines[index][16:]
    with pytest.raises(DiffusionContractError, match="atom identities"):
        audit_completion_candidate(source, candidate.replace(lines[index], lines[index] + extra), request)


def test_duplicate_missing_and_nonfinite_atoms_fail_closed(completion_control):
    source, candidate, request = completion_control
    generated = next(line for line in candidate.splitlines(keepends=True)
                     if line[12:16].strip() == b"NZ" and line[22:26].strip() == b"27")
    for forged in (candidate.replace(generated, b""), candidate.replace(b"END\n", generated + b"END\n"),
                   candidate.replace(generated, generated[:30] + b"     nan" + generated[38:])):
        with pytest.raises((DiffusionContractError, ValueError, RuntimeError)):
            audit_completion_candidate(source, forged, request)


def test_bad_generated_geometry_is_reported_not_repaired_by_validator(completion_control):
    source, candidate, request = completion_control
    original = candidate
    audit = audit_completion_candidate(source, candidate, request)
    coordinates = list(audit.generated_coordinates)
    coordinates[-1] = (80, 80, 80)
    bad = materialize_atom_completion(source, request, request.generated_atoms, coordinates)
    report = validate_completion_candidate(source, bad, request)
    assert not report.summary.passed
    assert "completion-bond_length" in report.summary.hard_gate_failures
    assert report.source_findings
    assert candidate == original


def test_exact_case_and_insertion_aware_finding_scope(completion_control):
    from dvbfixer.diagnose.report import Finding, Severity

    source, _, request = completion_control
    assert ("A", "27") in _affected_residue_locations(request)
    finding = Finding(Severity.ERROR, "clash", "D", "82", "ALA", "partner clash",
                      extra={"atoms": [{"chain": "d", "resid": "82", "icode": "A"}]})
    assert _touches(finding, {("d", "82A")})
    assert not _touches(finding, {("D", "82A"), ("d", "82")})
    with pytest.raises(ValueError, match="digest mismatch"):
        audit_completion_candidate(source + b"\n", source, request)


def test_refinement_masks_and_seed_are_not_caller_mutable(completion_control, monkeypatch):
    import dvbfixer.model.diffusion.atom_completion_refinement as adapter

    source, candidate, request = completion_control
    with pytest.raises(DiffusionContractError, match="seed"):
        adapter.refine_completion_candidate(source, candidate, request, seed=101)

    def fake_kernel(text, **kwargs):
        from types import SimpleNamespace

        assert kwargs["generated_atoms"] == request.generated_atoms
        assert kwargs["generated_residues"] == (ResidueIdentity("A", "27"),)
        assert "observed_anchor_restraints" not in kwargs
        return SimpleNamespace(coordinates_angstrom={request.masks.fixed_atoms[0]: (0, 0, 0)})

    monkeypatch.setattr(adapter, "refine_generated_region", fake_kernel)
    with pytest.raises(adapter.BoundaryRefinementError, match="exact generated mask"):
        adapter.refine_completion_candidate(source, candidate, request, seed=7)


def test_actual_generated_only_reference_refinement_control(completion_control):
    from dvbfixer.model.diffusion.atom_completion_refinement import refine_completion_candidate

    source, candidate, request = completion_control
    result = refine_completion_candidate(source, candidate, request, seed=7, max_iterations=50, restart_count=0)
    assert result.platform == "Reference"
    assert result.raw_candidate == candidate
    audit_completion_candidate(source, result.refined_candidate, request)
    assert result.refined_validation.summary.passed, result.refined_validation.summary.hard_gate_failures
    assert result.final_energy_kj_mol <= result.initial_energy_kj_mol
    assert result.raw_validation.source_findings == result.refined_validation.source_findings


def test_unrelated_source_geometry_damage_stays_visible(completion_control):
    source, candidate, request = completion_control
    original_line = source.splitlines(keepends=True)[0]
    changed_line = original_line[:30] + b"  80.000  80.000  80.000" + original_line[54:]
    source = source.replace(original_line, changed_line)
    candidate = candidate.replace(original_line, changed_line)
    digest = hashlib.sha256(source).hexdigest()
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                      masks=replace(request.masks, source_sha256=digest))
    report = validate_completion_candidate(source, candidate, request)
    assert any(finding.category == "bond_length" and finding.resid == "1" for finding in report.source_findings)
    assert any("source ERROR: A/1 bond_length" in warning for warning in report.summary.warnings)
    assert any(finding.category == "bond_length" and finding.resid == "1" for finding in report.candidate_findings)


def test_mixed_gap_and_partial_reference_control(completion_control, tmp_path):
    from dvbfixer.model.diffusion.atom_completion import plan_request_atom_completion
    from dvbfixer.model.diffusion.atom_completion_refinement import refine_completion_candidate
    from dvbfixer.model.diffusion.atom_completion_validation import _coordinates
    from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

    source, reference, partial_request = completion_control
    source = b"".join(line for line in source.splitlines(keepends=True)
                      if not (line.startswith(b"ATOM  ") and 10 <= int(line[22:26]) <= 13))
    path = tmp_path / "mixed-source.pdb"
    path.write_bytes(source)
    gap_request = build_cli_diffusion_request(path, {"A": partial_request.target_sequences[0].sequence},
                                              seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    masks = plan_request_atom_completion(gap_request, source)
    request = AtomCompletionRequest(1, gap_request.normalized_pdb, gap_request.target_sequences,
                                    gap_request.sequence_placements, gap_request.gaps, masks, (7,),
                                    "private-control-not-sampling")
    coordinates = _coordinates(reference.decode().splitlines())
    candidate = materialize_atom_completion(source, request, request.generated_atoms,
                                            [coordinates[atom] for atom in request.generated_atoms])
    raw = validate_completion_candidate(source, candidate, request)
    assert raw.summary.passed, raw.summary.hard_gate_failures
    result = refine_completion_candidate(source, candidate, request, seed=7, max_iterations=50, restart_count=0)
    assert result.refined_validation.summary.passed, result.refined_validation.summary.hard_gate_failures
    assert len(request.gaps) == 1 and len(request.masks.partial_residues) == 1
    audit_completion_candidate(source, result.refined_candidate, request)


def test_raw_d_candidate_refused_without_observed_sidechain_reflection(completion_control):
    from io import StringIO

    import numpy as np
    from openmm.app import PDBFile
    from openmm.unit import nanometer

    from dvbfixer.ffutils.geometry import _ca_triple
    from dvbfixer.model.diffusion.atom_completion_refinement import (
        BoundaryRefinementError,
        refine_completion_candidate,
    )

    source, candidate, request = completion_control
    # Reflect all observed sidechain atoms of MET1 to manufacture an explicitly
    # labelled defective source, not an authorized completion-generated change.
    pdb = PDBFile(StringIO(candidate.decode()))
    residue = next(pdb.topology.residues())
    atoms = {atom.name: atom for atom in residue.atoms()}
    xyz = np.asarray(pdb.positions.value_in_unit(nanometer))
    ca = xyz[atoms["CA"].index]
    normal = np.cross(xyz[atoms["N"].index] - ca, xyz[atoms["C"].index] - ca)
    normal /= np.linalg.norm(normal)
    original_lines = candidate.splitlines(keepends=True)
    replacements = {}
    for line in original_lines:
        if line.startswith(b"ATOM  ") and line[22:26].strip() == b"1" and line[12:16].strip() not in {b"N", b"CA", b"C", b"O"}:
            point = np.array([float(line[i:i + 8]) / 10 for i in (30, 38, 46)])
            reflected = (point - 2 * np.dot(point - ca, normal) * normal) * 10
            replacements[line] = line[:30] + "".join(f"{value:8.3f}" for value in reflected).encode() + line[54:]
    for old, new in replacements.items():
        source, candidate = source.replace(old, new), candidate.replace(old, new)
    digest = hashlib.sha256(source).hexdigest()
    request = replace(request, normalized_pdb=replace(request.normalized_pdb, sha256=digest),
                      masks=replace(request.masks, source_sha256=digest))
    # Ensure the control is genuinely D before relying on the refusal result.
    damaged = PDBFile(StringIO(candidate.decode()))
    first = {atom.name: atom for atom in next(damaged.topology.residues()).atoms()}
    assert _ca_triple(damaged.positions, first["N"], first["CA"], first["C"], first["CB"])[0] < 0
    report = validate_completion_candidate(source, candidate, request)
    assert "completion-d-ca-chirality" in report.summary.hard_gate_failures
    with pytest.raises(BoundaryRefinementError, match="raw D-Cα") as failure:
        refine_completion_candidate(source, candidate, request, seed=7)
    assert failure.value.raw_candidate == candidate
    assert not failure.value.raw_validation.summary.passed


def test_failed_kernel_retains_raw_attempt_and_diagnostics(completion_control, monkeypatch):
    import dvbfixer.model.diffusion.atom_completion_refinement as adapter

    source, candidate, request = completion_control

    def fail_kernel(*args, **kwargs):
        raise adapter.BoundaryRefinementError("declared failure control")

    monkeypatch.setattr(adapter, "refine_generated_region", fail_kernel)
    with pytest.raises(adapter.CompletionRefinementError, match="declared failure control") as failure:
        adapter.refine_completion_candidate(source, candidate, request, seed=7)
    assert failure.value.raw_candidate == candidate
    assert failure.value.raw_validation.source_findings
    assert failure.value.__cause__ is not None


def test_refined_geometry_rejection_is_retained_not_silently_accepted(completion_control, monkeypatch):
    from types import SimpleNamespace

    import dvbfixer.model.diffusion.atom_completion_refinement as adapter

    source, candidate, request = completion_control
    original = audit_completion_candidate(source, candidate, request)
    coordinates = dict(zip(request.generated_atoms, original.generated_coordinates, strict=True))
    coordinates[request.generated_atoms[-1]] = (80, 80, 80)
    monkeypatch.setattr(adapter, "refine_generated_region", lambda *args, **kwargs: SimpleNamespace(
        coordinates_angstrom=coordinates, initial_energy_kj_mol=1.0, final_energy_kj_mol=0.0, platform="mock-negative-control"))
    result = adapter.refine_completion_candidate(source, candidate, request, seed=7)
    assert result.raw_candidate == candidate
    assert result.raw_validation.summary.passed
    assert not result.refined_validation.summary.passed
    assert "completion-bond_length" in result.refined_validation.summary.hard_gate_failures
    assert result.refined_candidate != candidate
    audit_completion_candidate(source, result.refined_candidate, request)


def test_truncated_diagnostics_cannot_grant_a_private_pass(completion_control, monkeypatch):
    import dvbfixer.model.diffusion.atom_completion_validation as validation
    from dvbfixer.diagnose.report import Finding, Severity

    source, candidate, request = completion_control
    monkeypatch.setattr(validation, "_geometry_findings", lambda pdb: (
        Finding(Severity.INFO, "coincident_atoms", "", "", "", "Heavy-atom coincidence report truncated after 1000 pairs"),
    ))
    report = validation.validate_completion_candidate(source, candidate, request)
    assert not report.summary.passed
    assert "completion-diagnostics-incomplete" in report.summary.hard_gate_failures
