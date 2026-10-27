"""Two-stage controller invariants without importing any scientific backend."""

from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import pytest

from dvbfixer.model.diffusion.anchor_policy import NumericalRefinementError
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.refinement_stages import (
    NumericalFailureEvidence,
    RepairAssessment,
    RepairIssue,
    problematic_anchors,
    run_refinement_stages,
)

SOURCE = Path(__file__).parent / "fixtures/protein_context/9cqd_split.pdb"
LEFT, RIGHT = ResidueIdentity("W", "25"), ResidueIdentity("W", "29")
NEW = ResidueIdentity("W", "28")


def request():
    text = SOURCE.read_text()
    atoms = tuple(AtomIdentity(l[21], l[22:26].strip(), l[26].strip(), l[12:16].strip())
                  for l in text.splitlines() if l.startswith("ATOM  "))
    return SimpleNamespace(gaps=(SimpleNamespace(left_anchor=LEFT, right_anchor=RIGHT, generated_residues=(NEW,)),), fixed_atoms=atoms)


def bad(residues=(NEW, RIGHT), category="non_planar_amide", actionable=True):
    return RepairAssessment((RepairIssue(category, residues, "junction defect", actionable),))


def test_clean_raw_does_not_call_minimizer_or_warn():
    def forbidden(*_args):
        pytest.fail("clean input must not be minimized")
    result = run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda _: RepairAssessment(), minimize=forbidden, notify=forbidden)
    assert result.coordinate_text == "raw" and result.selected_stage == "raw" and len(result.attempts) == 1


def test_stage_one_success_never_releases_anchors():
    calls = []
    def minimize(text, restraints):
        calls.append((text, restraints))
        return "stage1"
    result = run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda t: bad() if t == "raw" else RepairAssessment(), minimize=minimize)
    assert calls == [("raw", ())]
    assert result.coordinate_text == "stage1" and not result.requires_anchor_acceptance_contract


def test_default_stage_two_restarts_from_raw_and_selects_only_problem_flank():
    calls, warnings = [], []
    def minimize(text, restraints):
        calls.append((text, restraints))
        return "stage2" if restraints else "stage1-mutated"
    result = run_refinement_stages("raw", SOURCE.read_text(), request(),
        assess=lambda t: RepairAssessment() if t == "stage2" else bad(), minimize=minimize, notify=warnings.append)
    assert [text for text, _ in calls] == ["raw", "raw"]
    assert result.selected_anchors == (RIGHT,)
    assert result.requires_anchor_acceptance_contract
    assert all(r.atom.residue_number == "29" for r in calls[1][1])
    assert len(warnings) == 1 and "W:29" in warnings[0]
    assert all(a.input_sha256 == result.attempts[0].input_sha256 for a in result.attempts)


def test_disable_prevents_stage_two_and_failed_output_is_not_selected():
    calls = []
    result = run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda _: bad(),
        minimize=lambda text, restraints: calls.append(text) or "failed", relax_anchors=False)
    assert calls == ["raw"] and result.coordinate_text is None


@pytest.mark.parametrize("assessment", [bad((NEW,), "bond_angle"), bad((RIGHT,), "clash"),
    bad((NEW, ResidueIdentity("X", "999")), "clash"), bad(actionable=False), bad(category="chirality")])
def test_internal_unrelated_warning_and_chirality_findings_do_not_release(assessment):
    assert problematic_anchors(request(), assessment) == ()


def test_generated_anchor_clash_selects_only_direct_named_anchor():
    assert problematic_anchors(request(), bad((NEW, LEFT), "clash")) == (LEFT,)


def test_terminal_shared_and_case_specific_anchor_identity():
    left = ResidueIdentity("w", "29", "A")
    generated = ResidueIdentity("w", "30")
    gap = SimpleNamespace(left_anchor=left, right_anchor=None, generated_residues=(generated,))
    view = SimpleNamespace(gaps=(gap, gap))
    assessment = bad((generated, left), "bond_angle")
    assert problematic_anchors(view, assessment) == (left,)
    assert problematic_anchors(view, bad((generated, RIGHT), "non_planar_amide")) == ()


def test_actual_ramachandran_diagnostic_category_can_release_named_anchor():
    assert problematic_anchors(request(), bad((RIGHT,), "ramachandran_outlier")) == (RIGHT,)


def test_non_geometry_admission_failure_blocks_minimization():
    def forbidden(*_args):
        pytest.fail("admission failure must not be minimized")
    result = run_refinement_stages("raw", SOURCE.read_text(), request(),
        assess=lambda _: RepairAssessment(admission_failures=("identity-mismatch",)), minimize=forbidden)
    assert result.selected_stage is None and len(result.attempts) == 1


def test_numerical_failure_uses_raw_diagnostics_but_never_raw_fallback_acceptance():
    calls = []
    def fail(text, restraints):
        calls.append(text)
        raise NumericalRefinementError("numerical failure")
    result = run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda _: bad(), minimize=fail)
    assert calls == ["raw", "raw"] and result.coordinate_text is None
    assert result.attempts[1].error.startswith("NumericalRefinementError:")
    assert result.attempts[2].error.startswith("NumericalRefinementError:")
    assert result.attempts[1].failure_kind == "numerical"
    failure = result.attempts[1].numerical_failure
    assert failure is not None and failure.input_sha256 == result.attempts[0].input_sha256
    assert asdict(failure)["error_type"] == "NumericalRefinementError"


@pytest.mark.parametrize("error", [RuntimeError("numerical-looking message"), ValueError("invalid template")])
def test_unclassified_failure_never_releases_anchors(error):
    def fail(*_args):
        raise error
    result = run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda _: bad(), minimize=fail)
    assert len(result.attempts) == 2 and result.selected_stage is None
    assert result.attempts[1].failure_kind == "unclassified"


def test_numerical_failure_codec_is_pure_bounded_and_versioned():
    data = {"input_sha256": "a" * 64, "error_type": "BoundaryNumericalError", "message": "No finite minimum"}
    record = NumericalFailureEvidence(**data)
    assert NumericalFailureEvidence(**asdict(record)) == record
    for mutation in ({"error_type": "ValueError"}, {"revision": True}, {"message": 42},
                     {"message": "x" * 2049}, {"input_sha256": "A" * 64}):
        with pytest.raises(ValueError):
            NumericalFailureEvidence(**(data | mutation))


def test_programming_errors_and_cancellation_are_not_repair_triggers():
    def fail(*_args):
        raise KeyError("bug")
    with pytest.raises(KeyError):
        run_refinement_stages("raw", SOURCE.read_text(), request(), assess=lambda _: bad(), minimize=fail)
