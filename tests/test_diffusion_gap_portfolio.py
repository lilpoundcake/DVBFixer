"""A locally clean gap is never whole-complex acceptance or a new seed."""
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.contract import (
    DiffusionContractError,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
)
from dvbfixer.model.diffusion.gap_portfolio import (
    GapChoice,
    GapParent,
    audit_portfolio_materialization,
    digest,
    plan_gap_repairs,
    repair_portfolio_gaps,
    search_gap_portfolio,
)
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, RepairIssue
from tests.test_diffusion_anchor_acceptance import shift
from tests.test_diffusion_validate import _atom_identity, _request, _source_and_candidate


def two_gaps():
    _, full, _, _ = _source_and_candidate()
    good = "".join(line for line in full.splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[21] == "C" and 64 <= int(line[22:26]) <= 70) + "END\n"
    source = "".join(line for line in good.splitlines(keepends=True)
        if not line.startswith("ATOM  ") or int(line[22:26]) not in (65, 69))
    generated = tuple(_atom_identity(line) for line in good.splitlines()
        if line.startswith("ATOM  ") and int(line[22:26]) in (65, 69))
    fixed = tuple(_atom_identity(line) for line in source.splitlines() if line.startswith("ATOM  "))
    request = _request(source.encode(), fixed, generated)
    def residue(n):
        return ResidueIdentity("C", str(n))
    request = replace(request,
        sequence_placements=(SequencePlacement("C", 7, (0, 2, 3, 4, 6), tuple(residue(n) for n in (64, 66, 67, 68, 70))),),
        gaps=tuple(GapRegion("C", TargetInterval(i, i + 1), residue(n - 1), residue(n + 1),
            (residue(n),), tuple(residue(k) for k in (n - 1, n, n + 1))) for n, i in ((65, 1), (69, 5))))
    return request, source, good


def parent(request, source, text, seed):
    return GapParent(f"seed-{seed}", seed, digest(request.to_json()), digest(source), digest(text),
        "a" * 64, "b" * 64, text)


def test_complementary_rejected_parents_can_yield_a_valid_multiparent_complex():
    request, source, good = two_gaps()
    parents = (parent(request, source, shift(good, "65", 1.5, atom_name="C"), 7),
               parent(request, source, shift(good, "69", 1.5, atom_name="C"), 11))
    retained = []
    result = search_gap_portfolio(request, source, parents, max_evaluations=8,
        retain_trial=lambda trial, text: retained.append((trial, text)))
    assert result.selected.score == (0, 0)
    assert [choice.parent_index for choice in result.selected.choices] == [1, 0]
    assert result.coordinate_text == good
    assert result.kind == "experimental-multiparent-gap-portfolio"
    assert not result.publication_qualified
    assert not result.repair_scope.gap_indices
    assert len(retained) == len(result.trials) <= 8
    assert all(digest(text) == trial.coordinate_sha256 for trial, text in retained)
    assert not hasattr(result, "seed")


@pytest.mark.parametrize("mutation", ["frame", "metadata", "observation", "duplicate", "stage"])
def test_incompatible_donors_fail_before_search(mutation):
    request, source, good = two_gaps()
    first = parent(request, source, good, 7)
    second = parent(request, source, good, 11)
    with pytest.raises(DiffusionContractError):
        if mutation == "frame":
            second = replace(second, source_sha256="0" * 64)
        elif mutation == "metadata":
            changed = good.replace("  1.00", "  0.50", 1)
            second = replace(second, pdb_text=changed, coordinate_sha256=digest(changed))
        elif mutation == "observation":
            changed = shift(good, "64", 0.1)
            second = replace(second, pdb_text=changed, coordinate_sha256=digest(changed))
        elif mutation == "duplicate":
            second = first
        else:
            second = replace(second, stage="generated-and-anchors")
        search_gap_portfolio(request, source, (first, second))


def test_budget_exhaustion_does_not_grant_local_acceptance():
    request, source, good = two_gaps()
    parents = (parent(request, source, shift(good, "65", 1.5, atom_name="C"), 7),
               parent(request, source, shift(good, "69", 1.5, atom_name="C"), 11))
    result = search_gap_portfolio(request, source, parents, max_evaluations=2)
    assert result.budget_exhausted and result.selected.score != (0, 0)
    assert result.repair_scope.gap_indices


def test_repair_scope_preserves_clean_gaps_but_unknown_findings_fail_closed():
    request, _, _ = two_gaps()
    issue = RepairIssue("bond_length", (ResidueIdentity("C", "65"),), "bad bond")
    aggregate = RepairIssue("generated-or-junction-bond-length", (), "gate")
    scope = plan_gap_repairs(request, RepairAssessment((issue, aggregate)))
    assert scope.gap_indices == (0,) and scope.fixed_gap_indices == (1,)
    scope = plan_gap_repairs(request, RepairAssessment((RepairIssue("unknown", (), "unattributed"),)))
    assert scope.gap_indices == (0, 1)


def test_source_digest_is_checked_even_when_all_parents_share_a_wrong_frame():
    request, source, good = two_gaps()
    wrong = source + "REMARK changed source\n"
    with pytest.raises(DiffusionContractError, match="source differs"):
        search_gap_portfolio(request, wrong, (parent(request, wrong, good, 7),))


def test_independent_audit_rejects_wrong_lineage_and_changed_headers():
    request, source, good = two_gaps()
    bad = shift(good, "65", 1.5, atom_name="C")
    parents = (parent(request, source, bad, 7), parent(request, source, good, 11))
    choices = (GapChoice(0, 1), GapChoice(1, 0))
    assert audit_portfolio_materialization(request, source, parents, choices, good).geometry_passed
    with pytest.raises(DiffusionContractError, match="selected parent"):
        audit_portfolio_materialization(request, source, parents,
            (GapChoice(0, 0), GapChoice(1, 0)), good)
    with pytest.raises(DiffusionContractError, match="non-coordinate"):
        audit_portfolio_materialization(request, source, parents, choices, "REMARK invented\n" + good)


def test_issue_local_refinement_cannot_move_the_other_generated_gap(monkeypatch):
    from types import SimpleNamespace

    from dvbfixer.model.diffusion.gap_portfolio import _coordinates

    request, source, good = two_gaps()
    bad = shift(good, "65", 1.5, atom_name="C")
    observed_calls = []

    def minimize(text, **kwargs):
        assert text == bad
        assert kwargs["generated_residues"] == (ResidueIdentity("C", "65"),)
        assert all(a.residue_number == "65" for a in kwargs["generated_atoms"])
        xyz, _ = _coordinates(good)
        selected = {a: xyz[a] for a in kwargs["generated_atoms"]}
        assert kwargs["candidate_score"](selected) == (0, 0)
        observed_calls.append(kwargs)
        return SimpleNamespace(coordinates_angstrom=selected, search_diagnostics=())

    monkeypatch.setattr("dvbfixer.model.diffusion.boundary_refinement.refine_generated_region", minimize)
    result = repair_portfolio_gaps(request, source, bad, refinement_seed=19)
    assert len(observed_calls) == 1 and result.coordinate_text == good
    assert result.scope.fixed_gap_indices == (1,)
    parents = (parent(request, source, bad, 7),)
    choices = (GapChoice(0, 0), GapChoice(1, 0))
    assert audit_portfolio_materialization(request, source, parents, choices, bad,
        repaired_text=good).geometry_passed
    with pytest.raises(DiffusionContractError, match="clean gap"):
        audit_portfolio_materialization(request, source, parents, choices, bad,
            repaired_text=shift(good, "69", 0.04))


def test_scope_keeps_chain_case_and_insertion_codes_distinct():
    from types import SimpleNamespace

    residues = (ResidueIdentity("D", "82", "A"), ResidueIdentity("d", "82", "A"), ResidueIdentity("D", "82"))
    request = SimpleNamespace(gaps=tuple(SimpleNamespace(generated_residues=(r,), left_anchor=None,
        right_anchor=None) for r in residues))
    scope = plan_gap_repairs(request, RepairAssessment((RepairIssue("clash", (residues[0],), "clash"),)))
    assert scope.gap_indices == (0,) and scope.fixed_gap_indices == (1, 2)
