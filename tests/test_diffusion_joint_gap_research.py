"""Joint transport is a separate research contract, not two batched predictions."""
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.contract import DiffusionContractError, GapKind
from dvbfixer.model.diffusion.joint_gap_research import (
    COMPLETE_REVISION,
    MOTIF_COMPONENT_REVISION,
    REVISION,
    build_joint_gap_view,
    complete_context_admission,
    motif_component_admission,
    plan_complete_context,
    plan_motif_component,
)
from dvbfixer.model.diffusion.protein_context import prepare_cc94_request, validate_cc94_request
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request
from tests.test_diffusion_input_preparation import _request


def two_gaps(tmp_path):
    path, source, old = _request(tmp_path)
    source = b"".join(line for line in source.splitlines(keepends=True)
        if not (line.startswith(b"ATOM  ") and 20 <= int(line[22:26]) <= 22))
    path.write_bytes(source)
    parent = build_cli_diffusion_request(path, {t.chain: t.sequence for t in old.target_sequences},
        seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    parent = prepare_cc94_request(parent, source, "auto")
    assert len(parent.gaps) == 2
    return parent, source


def test_one_shared_axis_preserves_both_gaps_and_all_fixed_identities(tmp_path):
    parent, source = two_gaps(tmp_path)
    original = parent.to_json()
    view = build_joint_gap_view(parent, source, (0, 1), seed=7)
    assert parent.to_json() == original
    assert len(view.request.gaps) == len(view.engine_view.request.gaps) == 2
    assert view.request.generated_atoms == parent.generated_atoms
    assert view.request.fixed_atoms == parent.fixed_atoms
    assert view.request.gap_contexts is None
    assert view.request.protein_context.policy.name == REVISION
    assert len(view.positions) <= 512
    assert set(view.chains) == {0}
    mapping = dict(view.engine_view.atom_to_source)
    assert len(mapping) == len(set(mapping.values()))
    assert {mapping[a] for a in view.engine_view.request.generated_atoms} == set(parent.generated_atoms)
    assert not set(parent.generated_atoms) & {c.atom for c in view.request.protein_context.fixed_conditioning}
    expected = {c.atom for p in parent.gap_contexts for c in p.fixed_conditioning}
    assert {c.atom for c in view.request.protein_context.fixed_conditioning} == expected
    # Production admission remains deliberately unchanged.
    with pytest.raises(ValueError):
        validate_cc94_request(view.request)


@pytest.mark.parametrize("indices", [(0,), (1, 0), (0, 0), (0, 2), (False, 1)])
def test_invalid_group_is_not_silently_split_or_reordered(tmp_path, indices):
    parent, source = two_gaps(tmp_path)
    with pytest.raises(DiffusionContractError):
        build_joint_gap_view(parent, source, indices, seed=7)


def test_changed_source_and_implicit_seed_fail_before_sampling(tmp_path):
    parent, source = two_gaps(tmp_path)
    with pytest.raises(DiffusionContractError, match="digest"):
        build_joint_gap_view(parent, b"REMARK changed\n" + source, (0, 1), seed=7)
    with pytest.raises(DiffusionContractError):
        build_joint_gap_view(parent, source, (0, 1), seed=True)


def test_terminal_gap_is_not_relabelled_internal(tmp_path):
    parent, source = two_gaps(tmp_path)
    gap = replace(parent.gaps[0], gap_kind=GapKind.N_TERMINAL, left_anchor=None)
    # Removing precomputed contexts prevents treating the malformed pilot as
    # an authorized frozen independent request in the first place.
    parent = replace(parent, gaps=(gap, parent.gaps[1]), gap_contexts=None)
    with pytest.raises(ValueError):
        build_joint_gap_view(parent, source, (0, 1), seed=7)


def distant_partner(tmp_path):
    parent, source = two_gaps(tmp_path)
    observed = [line for line in source.decode().splitlines(keepends=True) if line.startswith("ATOM  ")]
    partner = []
    for line in observed:
        if 1 <= int(line[22:26]) <= 8:
            partner.append(line[:6] + f"{int(line[6:11]) + 2000:5d}" + line[11:21] + "a" + line[22:30]
                           + f"{float(line[30:38]) + 40:8.3f}" + line[38:])
    source = ("".join(observed + partner) + "END\n").encode()
    path = tmp_path / "with-partner.pdb"
    path.write_bytes(source)
    sequences = {t.chain: t.sequence for t in parent.target_sequences}
    sequences["a"] = sequences["A"][:8]
    parent = build_cli_diffusion_request(path, sequences, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    parent = prepare_cc94_request(parent, source, "auto")
    assert not any(c.atom.chain == "a" for p in parent.gap_contexts for c in p.fixed_conditioning)
    return parent, source


def test_complete_context_keeps_distant_body_and_case_sensitive_partner(tmp_path):
    parent, source = distant_partner(tmp_path)
    before = parent.to_json()
    view = plan_complete_context(parent, source, (0, 1), seed=7)
    plan = view.request.protein_context
    assert plan.policy.name == COMPLETE_REVISION and plan.preserved_only == ()
    assert {c.atom for c in plan.fixed_conditioning} == set(parent.fixed_atoms)
    assert {s for _, s in view.engine_view.atom_to_source} == set(parent.fixed_atoms) | set(parent.generated_atoms)
    assert {c.atom.chain for c in plan.fixed_conditioning} == {"A", "a"}
    assert set(view.chains) == {0, 1}
    assert len(view.positions) == sum(len(t.sequence) for t in parent.target_sequences)
    assert complete_context_admission(view)["supported"]
    assert parent.to_json() == before
    with pytest.raises(ValueError):
        validate_cc94_request(view.request)


def test_complete_plan_refuses_execution_over_budget_instead_of_dropping_partner(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.joint_gap_research as module

    parent, source = distant_partner(tmp_path)
    view = plan_complete_context(parent, source, (0, 1), seed=7)
    monkeypatch.setattr(module, "COMPLETE_EXECUTION_TOKEN_LIMIT", len(view.positions) - 1)
    admission = complete_context_admission(view)
    assert not admission["supported"]
    assert "complete-context-exceeds-explicit-research-token-budget" in admission["reasons"]
    assert admission["conditioning_fixed_atoms"] == len(parent.fixed_atoms)
    assert admission["output_only_observed_atoms"] == 0


def test_complete_admission_rejects_hidden_observed_atom_even_with_rehashed_plan(tmp_path):
    parent, source = distant_partner(tmp_path)
    view = plan_complete_context(parent, source, (0, 1), seed=7)
    plan = view.request.protein_context
    hidden = plan.fixed_conditioning[-1].atom
    plan = replace(plan, fixed_conditioning=plan.fixed_conditioning[:-1], preserved_only=(hidden,))
    forged = replace(view, request=replace(view.request, protein_context=plan))
    assert "observed-context-incomplete" in complete_context_admission(forged)["reasons"]


def test_complete_partial_group_never_conditions_on_other_gap_predictions(tmp_path):
    parent, source = distant_partner(tmp_path)
    view = plan_complete_context(parent, source, (0,), seed=7)
    expected = {a for a in parent.generated_atoms if a in parent.gap_contexts[0].generated
                or a in parent.gap_contexts[0].movable_junction}
    assert set(view.request.generated_atoms) == expected
    assert {c.atom for c in view.request.protein_context.fixed_conditioning} == set(parent.fixed_atoms)
    assert not {c.atom for c in view.request.protein_context.fixed_conditioning} & set(parent.generated_atoms)


def test_old_local_union_cannot_pass_complete_context_admission(tmp_path):
    parent, source = distant_partner(tmp_path)
    legacy = build_joint_gap_view(parent, source, (0, 1), seed=7)
    assert "complete-context-revision-required" in complete_context_admission(legacy)["reasons"]


def test_complete_context_rejects_changed_source(tmp_path):
    parent, source = distant_partner(tmp_path)
    with pytest.raises(ValueError, match="digest"):
        plan_complete_context(parent, source + b"\n", (0, 1), seed=7)


def test_complete_context_rejects_untracked_source_atom(tmp_path):
    import hashlib

    parent, source = distant_partner(tmp_path)
    line = next(line for line in source.decode().splitlines(keepends=True) if line.startswith("ATOM  "))
    extra = line[:6] + "99999" + line[11:12] + " XX " + line[16:]
    changed = source + extra.encode()
    reference = replace(parent.normalized_pdb, sha256=hashlib.sha256(changed).hexdigest())
    parent = replace(parent, normalized_pdb=reference,
        gap_contexts=tuple(replace(p, source_sha256=reference.sha256) for p in parent.gap_contexts))
    with pytest.raises(ValueError, match="untracked observed atoms"):
        plan_complete_context(parent, changed, (0, 1), seed=7)


def test_complete_engine_mapping_rejects_missing_observation(tmp_path):
    parent, source = distant_partner(tmp_path)
    view = plan_complete_context(parent, source, (0, 1), seed=7)
    mapping = view.engine_view.atom_to_source
    hidden = next(s for _, s in mapping if s in parent.fixed_atoms)
    with pytest.raises(ValueError, match="cover every represented atom"):
        replace(view.engine_view, atom_to_source=tuple((e, s) for e, s in mapping if s != hidden))


def test_neighbor_motif_preserves_local_conditioning_and_joint_generated_mask(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.gap_neighbors as neighbors

    parent, source = distant_partner(tmp_path)
    monkeypatch.setattr(neighbors, "source_gap_components", lambda *_: (((0, 1),), ()))
    view = plan_motif_component(parent, source, (0, 1), seed=7)
    assert view.request.protein_context.policy.name == MOTIF_COMPONENT_REVISION
    assert set(view.request.generated_atoms) == set(parent.generated_atoms)
    expected = {c.atom for plan in parent.gap_contexts for c in plan.fixed_conditioning}
    assert {c.atom for c in view.request.protein_context.fixed_conditioning} == expected
    assert not any(c.atom.chain == "a" for c in view.request.protein_context.fixed_conditioning)
    assert motif_component_admission(view)["supported"]
    with pytest.raises(ValueError, match="complete neighboring-gap component"):
        plan_motif_component(parent, source, (0,), seed=7)


def test_separate_components_never_reuse_other_generated_predictions(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.gap_neighbors as neighbors

    parent, source = two_gaps(tmp_path)
    monkeypatch.setattr(neighbors, "source_gap_components", lambda *_: (((0,), (1,)), ()))
    views = [plan_motif_component(parent, source, (i,), seed=7) for i in range(2)]
    assert set(views[0].request.generated_atoms).isdisjoint(views[1].request.generated_atoms)
    for view in views:
        assert set(view.request.generated_atoms).isdisjoint(c.atom for c in view.request.protein_context.fixed_conditioning)
        assert {s for _, s in view.engine_view.atom_to_source} <= set(parent.fixed_atoms) | set(view.request.generated_atoms)
    with pytest.raises(ValueError, match="complete neighboring-gap component"):
        plan_motif_component(parent, source, (0, 1), seed=7)


def test_neighbor_component_budget_refusal_never_splits(tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.gap_neighbors as neighbors
    import dvbfixer.model.diffusion.joint_gap_research as transport

    parent, source = two_gaps(tmp_path)
    monkeypatch.setattr(neighbors, "source_gap_components", lambda *_: (((0, 1),), ()))
    view = plan_motif_component(parent, source, (0, 1), seed=7)
    monkeypatch.setattr(transport, "MOTIF_COMPONENT_TOKEN_LIMIT", len(view.positions) - 1)
    assert not motif_component_admission(view)["supported"]
    assert len(view.request.gaps) == 2
