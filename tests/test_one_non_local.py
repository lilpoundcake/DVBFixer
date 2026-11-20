"""Private local shell and movement-policy controls; no engines/checkpoints."""
from types import SimpleNamespace

import pytest
from scripts.local_gap_rosetta_worker import (
    checked_budget,
    checked_junction_restraint,
    resume_coordinates,
)
from scripts.one_non_local_research import local_mask, movements, scientific_status

from dvbfixer.model.diffusion.contract import ResidueIdentity


def case():
    rs = [ResidueIdentity("D", str(i), "") for i in range(1, 7)]
    request = SimpleNamespace(target_sequences=[SimpleNamespace(chain="D", sequence="AAAAAA")],
        sequence_placements=[SimpleNamespace(chain="D", observed_target_indices=[0, 1, 4, 5],
                                            observed_residues=[rs[i] for i in (0, 1, 4, 5)])],
        gaps=[SimpleNamespace(chain="D", left_anchor=rs[1], right_anchor=rs[4], generated_residues=rs[2:4],
                              target_interval=SimpleNamespace(start=2, stop=4))])
    source = {}
    for r in (rs[0], rs[1], rs[4], rs[5], ResidueIdentity("d", "2", "A"), ResidueIdentity("d", "90", "")):
        p = [float(r.residue_number), 0., 0.] if r.chain == "D" else [2., 1. if r.insertion_code else 100., 0.]
        for name in ("N", "CA", "C", "O", "CB"):
            source[r.chain, r.residue_number, r.insertion_code, name] = {"xyz": p.copy()}
    return request, source


def test_shell_is_source_only_case_sensitive_and_includes_contacting_chain():
    request, source = case()
    selected = set(local_mask(request, source))
    assert ("D", "2", "", "CB") in selected
    assert ("d", "2", "A", "CB") in selected
    assert ("d", "90", "", "CB") not in selected
    assert not any(k[1] in {"3", "4"} and k[0] == "D" for k in selected)


def test_bounds_are_per_residue_not_complex_average_and_no_fit():
    _, source = case()
    selected = {a for a in source if a[:3] == ("D", "2", "")}
    candidate = {a: {"xyz": r["xyz"].copy()} for a, r in source.items()}
    candidate["D", "2", "", "CA"]["xyz"][0] += .51
    rows, outside = movements(source, candidate, selected)
    assert outside == 0 and rows[0]["CA_A"] == pytest.approx(.51)
    assert rows[0]["within_limits"] is False
    candidate["d", "90", "", "CB"]["xyz"][2] += .001
    _, outside = movements(source, candidate, selected)
    assert outside == pytest.approx(.001)


def test_sidechain_movement_does_not_reclassify_observations_as_generated():
    _, source = case()
    selected = {a for a in source if a[:3] == ("D", "2", "")}
    candidate = {a: {"xyz": r["xyz"].copy()} for a, r in source.items()}
    candidate["D", "2", "", "CB"]["xyz"][0] += 1.26
    rows, _ = movements(source, candidate, selected)
    assert rows[0]["backbone_rms_A"] == 0 and rows[0]["within_limits"] is False


@pytest.mark.parametrize("geometry,coordinates,failures", [(False, True, ()), (True, False, ()),
                                                         (True, True, ("identity-changed",))])
def test_completed_process_cannot_override_scientific_rejection(geometry, coordinates, failures):
    assert scientific_status(geometry, coordinates, failures) == (False, "rejected")


def test_acceptance_requires_all_scientific_gates_not_just_no_engine_failure():
    assert scientific_status(True, True, ()) == (True, "accepted")


def test_ordinary_multichain_template_range_includes_last_case_sensitive_chain():
    from scripts.run_retained_loop_modeller_controls import multichain_alignment

    request, _ = case()
    request.target_sequences.append(SimpleNamespace(chain="d", sequence="AAAAAA"))
    request.sequence_placements.append(SimpleNamespace(chain="d", observed_target_indices=list(range(6)),
        observed_residues=[ResidueIdentity("d", str(i), "") for i in range(1, 7)]))
    text = multichain_alignment(request)
    assert "structureX:source:FIRST:D:LAST:d::::" in text
    assert "FIRST:@" not in text


def test_extra_budget_is_one_explicit_local_resume_not_a_default_change():
    assert checked_budget({"policy": "fixed", "seed": 101, "max_iterations": 500}) == (False, 500)
    assert checked_budget({"policy": "local", "seed": 101, "max_iterations": 1000, "resume_from": "saved"}) == (True, 1000)
    for payload in ({"policy": "fixed", "seed": 101, "max_iterations": 1000, "resume_from": "saved"},
                    {"policy": "local", "seed": 101, "max_iterations": 1000}):
        with pytest.raises(ValueError, match="undeclared"):
            checked_budget(payload)


def test_resume_requires_exact_native_axis_and_finite_complete_state():
    metadata = [{"identity": ["d", "82", "A", "CA"], "mobile": False}]
    points = [[1., 2., 3.]]
    assert resume_coordinates(metadata, {"atoms": metadata}, {"max_iterations": 500, "coordinates": points}) == points
    with pytest.raises(ValueError, match="identical native"):
        resume_coordinates(metadata, {"atoms": [{**metadata[0], "mobile": True}]}, {"max_iterations": 500, "coordinates": points})
    with pytest.raises(ValueError, match="invalid native"):
        resume_coordinates(metadata, {"atoms": metadata}, {"max_iterations": 500, "coordinates": [[float("nan"), 2, 3]]})


def test_junction_control_is_only_the_declared_case_sensitive_peptide_pair():
    spec = {"first": ["B", "90", "", "C"], "second": ["B", "91", "", "N"],
            "target_A": 1.33, "sd_A": .02, "weight": 1.}
    payload = {"policy": "local", "seed": 101, "max_iterations": 1000, "resume_from": "polished",
               "junction_length_restraint": spec}
    assert checked_budget(payload) == (True, 1000)
    assert checked_junction_restraint(payload) == spec
    for changed in ({**spec, "sd_A": .1}, {**spec, "first": ["b", "90", "", "C"]}):
        with pytest.raises(ValueError, match="undeclared junction"):
            checked_junction_restraint({**payload, "junction_length_restraint": changed})
    with pytest.raises(ValueError, match="undeclared refinement"):
        checked_budget({**payload, "resume_from": None, "max_iterations": 500})


def test_targeted_resume_declares_its_parent_iteration_stage_explicitly():
    metadata = [{"identity": ["B", "90", "", "C"], "mobile": True}]
    state = {"max_iterations": 1000, "coordinates": [[1., 2., 3.]]}
    with pytest.raises(ValueError, match="declared parent"):
        resume_coordinates(metadata, {"atoms": metadata}, state)
    assert resume_coordinates(metadata, {"atoms": metadata}, state, expected_parent_cap=1000) == [[1., 2., 3.]]


def test_budget_control_is_explicit_bounded_and_cannot_add_a_junction_term():
    payload = {"policy": "local", "seed": 101, "max_iterations": 3000, "resume_from": "polished",
               "budget_control": {"passes": 3, "iterations_per_pass": 3000}}
    assert checked_budget(payload) == (True, 3000)
    for change in ({"budget_control": {"passes": 4, "iterations_per_pass": 3000}},
                   {"resume_from": None}, {"policy": "fixed"}):
        with pytest.raises(ValueError, match="undeclared budget"):
            checked_budget({**payload, **change})
    with pytest.raises(ValueError, match="undeclared refinement"):
        checked_budget({k: v for k, v in payload.items() if k != "budget_control"})
    with pytest.raises(ValueError, match="undeclared budget"):
        checked_budget({**payload, "junction_length_restraint": {
            "first": ["B", "90", "", "C"], "second": ["B", "91", "", "N"],
            "target_A": 1.33, "sd_A": .02, "weight": 1.}})
