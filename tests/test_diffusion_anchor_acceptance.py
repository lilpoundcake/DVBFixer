"""Independent staged-repair auditing must not trust runner masks or scores."""

import hashlib
from dataclasses import replace

import pytest

pytest.importorskip("openmm")
pytest.importorskip("pdbfixer")

from dvbfixer.model.diffusion.anchor_acceptance import (
    _findings,
    _observed_auxiliary,
    _structure,
    audit_staged_anchor_repair,
)
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    ResidueIdentity,
    SequencePlacement,
    TargetSequence,
)
from dvbfixer.model.diffusion.refinement_stages import NumericalFailureEvidence
from tests.test_diffusion_validate import _request, _source_and_candidate


def shift(text, number, amount, *, atom_name=None):
    lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith("ATOM  ") and line[21] == "C" and line[22:26].strip() == number and (
            atom_name is None or line[12:16].strip() == atom_name
        ):
            line = line[:30] + f"{float(line[30:38]) + amount:8.3f}" + line[38:]
        lines.append(line)
    return "".join(lines)


def evidence():
    source, good, fixed, generated = _source_and_candidate()
    request = _request(source.encode(), fixed, generated)
    raw = shift(good, "69", 1.5, atom_name="C")
    final = shift(good, "70", 0.04)
    sha = hashlib.sha256(raw.encode()).hexdigest()
    return dict(request=request, source_text=source, raw_text=raw, stage1_text=raw,
                candidate_text=final, stage1_input_sha256=sha, stage2_input_sha256=sha,
                limits=AnchorMovementLimits(0.1, 0.2))


def test_authorization_is_recomputed_and_historical_strict_remains_failed():
    audit = audit_staged_anchor_repair(**evidence())
    assert audit.authorized_anchors == (ResidueIdentity("C", "70"),)
    assert audit.summary.passed, audit.summary.hard_gate_failures
    assert not audit.historical_strict_summary.passed
    assert "fixed-heavy-atom-max-displacement" in audit.historical_strict_summary.hard_gate_failures
    assert not audit.publication_qualified
    assert audit.movements[0].backbone_rmsd_angstrom == pytest.approx(0.04)
    assert any("C:70" in warning and "modified" in warning for warning in audit.summary.warnings)


def test_explicit_per_anchor_bounds_are_enforced():
    audit = audit_staged_anchor_repair(**(evidence() | {"limits": AnchorMovementLimits(0.02, 0.02)}))
    assert "anchor-repair-per-anchor-movement-limit" in audit.summary.hard_gate_failures


def test_default_audit_records_approved_limits_and_ca_movement():
    values = evidence()
    values.pop("limits")
    audit = audit_staged_anchor_repair(**values)
    assert audit.summary.passed
    assert audit.movement_limits == AnchorMovementLimits()
    assert audit.movements[0].ca_displacement_angstrom == pytest.approx(0.04)
    assert any("Cα displacement" in warning for warning in audit.summary.warnings)
    assert audit.policy_revision == "independent-staged-anchor-audit-v3"


@pytest.mark.parametrize("field", ["stage1_input_sha256", "stage2_input_sha256"])
def test_stage_inputs_must_be_bound_to_raw(field):
    audit = audit_staged_anchor_repair(**(evidence() | {field: "0" * 64}))
    assert "anchor-repair-not-restarted-from-raw" in audit.summary.hard_gate_failures


def test_unselected_anchor_motion_rejected_even_below_historical_tolerance():
    values = evidence()
    values["candidate_text"] = shift(values["candidate_text"], "64", 0.001, atom_name="CA")
    audit = audit_staged_anchor_repair(**values)
    assert "anchor-repair-unselected-observed-record-changed" in audit.summary.hard_gate_failures


def test_stage_one_may_not_move_observed_atoms_even_below_old_tolerance():
    values = evidence()
    values["stage1_text"] = shift(values["stage1_text"], "70", 0.001, atom_name="CA")
    audit = audit_staged_anchor_repair(**values)
    assert "anchor-repair-input-observed-record-changed" in audit.summary.hard_gate_failures


def test_observed_bfactor_metadata_must_not_be_rewritten():
    values = evidence()
    values["candidate_text"] = "".join(l[:60] + " 99.00" + l[66:] if l.startswith("ATOM  ") and l[21] == "C"
        and l[22:26].strip() == "70" else l for l in values["candidate_text"].splitlines(keepends=True))
    audit = audit_staged_anchor_repair(**values)
    assert "anchor-repair-observed-metadata-changed" in audit.summary.hard_gate_failures


def test_clean_stage_one_cannot_authorize_optional_anchor_motion():
    values = evidence()
    _source, good, _fixed, _generated = _source_and_candidate()
    audit = audit_staged_anchor_repair(**(values | {"stage1_text": good}))
    assert "anchor-repair-unnecessary-stage1-clean" in audit.summary.hard_gate_failures
    assert not audit.authorized_anchors


def test_clean_raw_cannot_be_escalated_by_fabricated_stage_one_errors():
    values = evidence()
    _source, good, _fixed, _generated = _source_and_candidate()
    sha = hashlib.sha256(good.encode()).hexdigest()
    audit = audit_staged_anchor_repair(**(values | {"raw_text": good,
        "stage1_input_sha256": sha, "stage2_input_sha256": sha}))
    assert "anchor-repair-unnecessary-raw-clean" in audit.summary.hard_gate_failures


def test_experimental_auxiliary_records_keep_exact_case_and_insertion_identity():
    source, _good, _fixed, _generated = _source_and_candidate()
    atom = AtomIdentity("C", "70", "", "CA")
    line = next(l for l in source.splitlines() if l.startswith("ATOM  ") and l[21] == "C"
                and l[22:26].strip() == "70" and l[12:16].strip() == "CA")
    auxiliary = "ANISOU" + line[6:28] + "   1000   1000   1000      0      0      0       C  "
    assert _observed_auxiliary(auxiliary, {atom}) == (auxiliary,)
    assert _observed_auxiliary(auxiliary, {AtomIdentity("c", "70", "", "CA")}) == ()
    assert _observed_auxiliary(auxiliary, {AtomIdentity("C", "70", "A", "CA")}) == ()


def test_source_digest_is_independently_verified():
    with pytest.raises(DiffusionContractError, match="source digest mismatch"):
        audit_staged_anchor_repair(**(evidence() | {"source_text": "REMARK changed\n" + evidence()["source_text"]}))


def test_typed_numerical_failure_uses_independent_raw_geometry_eligibility():
    values = evidence()
    failure = NumericalFailureEvidence(values["stage1_input_sha256"], "BoundaryNumericalError", "No finite all-L minimum")
    audit = audit_staged_anchor_repair(**(values | {"stage1_text": None, "stage1_failure": failure}))
    assert audit.summary.passed, audit.summary.hard_gate_failures
    assert ResidueIdentity("C", "70") in audit.authorized_anchors
    assert audit.stage1_sha256 is None and audit.stage1_failure == failure
    assert not audit.publication_qualified


def test_missing_or_misbound_failure_evidence_does_not_authorize_recovery():
    values = evidence() | {"stage1_text": None}
    missing = audit_staged_anchor_repair(**values)
    assert "anchor-repair-stage1-failure-evidence-missing" in missing.summary.hard_gate_failures
    wrong = NumericalFailureEvidence("0" * 64, "NumericalRefinementError", "No finite minimum")
    audit = audit_staged_anchor_repair(**(values | {"stage1_failure": wrong}))
    assert "anchor-repair-stage1-failure-input-mismatch" in audit.summary.hard_gate_failures


def test_failure_and_stage_one_coordinates_are_mutually_exclusive():
    values = evidence()
    failure = NumericalFailureEvidence(values["stage1_input_sha256"], "BoundaryNumericalError", "No finite minimum")
    with pytest.raises(DiffusionContractError, match="both coordinates and numerical failure"):
        audit_staged_anchor_repair(**(values | {"stage1_failure": failure}))


@pytest.mark.parametrize("kwargs", [{"error_type": "ValueError"}, {"message": "x" * 2049},
    {"input_sha256": "x" * 64}, {"revision": 2}])
def test_failure_evidence_is_bounded_and_rejects_unclassified_errors(kwargs):
    defaults = {"input_sha256": "0" * 64, "error_type": "BoundaryNumericalError", "message": "No finite minimum"}
    with pytest.raises(ValueError):
        NumericalFailureEvidence(**(defaults | kwargs))


def test_outer_anchor_junction_is_explicitly_checked_even_beyond_generated_region():
    values = evidence()
    request = replace(values["request"], target_sequences=(TargetSequence("C", "FSGSKSGAA"),),
        sequence_placements=(SequencePlacement("C", 9, (0, 6, 7, 8),
            tuple(ResidueIdentity("C", n) for n in ("64", "70", "71", "72"))),))
    _source, good, _fixed, _generated = _source_and_candidate()
    candidate = _structure(shift(good, "70", 5, atom_name="C"), "outer-defect.pdb")
    findings = _findings(candidate, request)
    assert any(issue.category == "junction-peptide-connectivity" and issue.residues ==
               (ResidueIdentity("C", "70"), ResidueIdentity("C", "71")) for issue in findings)
