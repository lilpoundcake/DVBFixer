"""Matched physical audit, not own-generator policy scoring of the comparator."""
from __future__ import annotations

import hashlib
import math
from dataclasses import replace

import pytest
from scripts.branched_template_baseline import preflight, template_proposal
from scripts.compare_branched_completion import measure, score
from scripts.summarize_branched_comparison import summarize

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.geometry import extend
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate


def control():
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    source = ("".join("ATOM  "+_format_generated_atom(serial=i, name=n, element=n[0], position=xyz,
                 residue_name="THR", chain_id="A", sequence_number="1", insertion_code="")[6:]
                 for i, (n, xyz) in enumerate(points.items(), 1))+"TER\nEND\n").encode()
    plan = plan_branched_source(source, ComponentInstanceRef(1, "A", "1"))
    return source, plan.request


def test_template_resource_preflight_is_actual_usability_not_import_only():
    evidence = preflight()
    assert evidence["platform"] == "Reference"
    assert len(evidence["resource_sha256"]) == 31
    assert {"soft.xml", "templates/VAL.pdb", "templates/THR.pdb"} <= evidence["resource_sha256"].keys()


def test_nonstaggered_proposal_is_not_penalized_by_the_own_sampling_grid():
    source, request = control()
    proposal = complete_sidechain(request).attempts[1].candidate  # +60 control
    points = {o.atom.name: o.xyz for o in request.observations} | {a.atom.name: a.xyz for a in proposal.additions}
    changed = []
    for a in proposal.additions:
        xyz = a.xyz
        if a.atom.name == "OG1":
            xyz = extend(points["N"], points["CA"], points["CB"], 1.43, 109.47122063449069, 61.)
        if a.atom.name == "CG2":
            xyz = extend(points["N"], points["CA"], points["CB"], 1.53, 109.47122063449069, -59.)
        changed.append(replace(a, xyz=xyz))
    candidate = replace(proposal, additions=tuple(changed), torsion_degrees=61.)
    assert not audit_candidate(request, candidate).accepted
    audit, output = measure(source, request, tuple(changed), anchored=False)
    assert audit.accepted and audit.observed_preserved
    assert audit.rdkit_target_cip == {"CA": "S", "CB": "R"}
    assert b"ATOM  " in output


def test_shared_seeded_baseline_retains_states_and_keeps_observed_atoms_fixed():
    source, request = control()
    states = []
    proposal = template_proposal(source, request, retain_snapshot=states.append)
    assert states and tuple(states) == proposal.snapshots
    assert [s.seed for s in states] == list(range(1, len(states)+1))
    assert len(states) <= 5 and all(s.observed_max_movement_angstrom < 1e-6 for s in states)
    assert {a.atom.name for a in proposal.additions} == {"CB", "OG1", "CG2"}
    assert all(a.name == "OXT" for a in proposal.temporary_extra_atoms)
    audit, output = measure(source, request, proposal.additions, anchored=False)
    assert isinstance(audit.accepted, bool) and audit.observed_preserved
    assert output.count(b"ATOM  ") == 7  # Unrequested OXT never enters the diagnostic proposal.
    assert hashlib.sha256(source).hexdigest() == request.source_sha256


def test_exact_name_score_is_separate_and_never_symmetry_minimized():
    _, request = control()
    additions = complete_sidechain(request).selected.additions
    truth = {a.atom.name: a.xyz for a in additions}
    assert score(additions, truth)["exact_name_rmsd_angstrom"] == 0
    shifted = dict(truth)
    shifted["CB"] = tuple(x+1 for x in shifted["CB"])
    result = score(additions, shifted)
    assert result["exact_name_rmsd_angstrom"] == pytest.approx(1.)
    assert result["maximum_error_angstrom"] == pytest.approx(math.sqrt(3))
    assert not result["used_for_selection"] and not result["symmetry_minimized"]


def test_baseline_checks_source_hash_before_any_rebuild():
    source, request = control()
    with pytest.raises(ValueError, match="SHA mismatch"):
        template_proposal(source, replace(request, source_sha256="0"*64))


def test_summary_keeps_failures_and_reports_both_arms_on_the_same_passing_subset():
    def row(own_error, template_error, passed):
        audit = {"accepted": True, "raw": {"findings": []}}
        return {"own": {"common_audit": audit, "scores": {"exact_name_rmsd_angstrom": own_error}},
                "template": {"common_audit": {**audit, "accepted": passed},
                             "scores": {"exact_name_rmsd_angstrom": template_error},
                             "helper_attempts": 1, "observed_temporary_max_movement_angstrom": 0.}}

    record = {"results": {"paired-1": row(.2, .1, True), "paired-2": row(.3, 9., False),
                           "anchored-1": row(.4, .2, True)}}
    result = summarize(record)
    pair = result["modes"]["paired"]
    assert pair["arms"]["template"]["raw_proposal_median_rmsd_angstrom_including_rejections"] == pytest.approx(4.55)
    assert pair["arms"]["own"]["both_pass_same_subset_median_rmsd_angstrom"] == .2
    assert pair["arms"]["template"]["both_pass_same_subset_median_rmsd_angstrom"] == .1
    assert pair["both_pass_count"] == 1
    assert pair["own_lower_raw_error_count_including_rejections"] == 1
    assert pair["own_lower_error_count_in_both_pass_subset"] == 0


def test_summary_does_not_count_all_nearplanar_gate_failures_as_physical_inversions():
    def row(cip):
        return {"own": {"common_audit": {"accepted": True}, "scores": {"exact_name_rmsd_angstrom": .1}},
                "template": {"common_audit": {"accepted": False,
                    "raw": {"findings": ["non_L_or_nearplanar_ca", "thr_cb_not_R_or_nearplanar"]},
                    "rdkit_target_cip": cip}, "scores": {"exact_name_rmsd_angstrom": .5},
                    "helper_attempts": 1, "observed_temporary_max_movement_angstrom": 0.}}

    summary = summarize({"results": {"paired-nearplanar": row({"CA": "S"}),
                                       "anchored-inverted": row({"CA": "R", "CB": "S"})}})
    assert summary["template_raw_non_L_or_nearplanar_ca_gate_count"] == 2
    assert summary["template_rounded_independent_CA_R_count"] == 1
    assert summary["template_raw_non_R_or_nearplanar_thr_cb_gate_count"] == 2
    assert summary["template_rounded_independent_CB_S_count"] == 1
