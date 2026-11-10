"""Explicit anchored singleton kernel; old capability refusals remain intact."""
from __future__ import annotations

import math
from dataclasses import replace

import pytest

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.anchored_branch import (
    audit_anchored_branch,
    complete_anchored_branch,
)
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition, CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import (
    add,
    cross,
    dot,
    extend,
    l_cb,
    scaled,
    sub,
    two_angle_tip,
    unit,
)
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate

TARGET = ComponentInstanceRef(1, "a", "82", "A")
BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def anchored_request(name, missing, chi=37.2):
    branch = BRANCHES[name]
    points = dict(BACKBONE)
    cb = points["CB"] = l_cb(points["N"], points["CA"], points["C"])
    points[branch.first] = extend(points["N"], points["CA"], cb, branch.first_length, TETRAHEDRAL_ANGLE, chi)
    points[branch.second] = extend(points["N"], points["CA"], cb, branch.second_length,
                                   TETRAHEDRAL_ANGLE, chi+branch.second_chi_offset)
    truth = points.pop(missing)
    request = CompletionRequest("a"*64, TARGET, name,
        tuple(Observation(ExactAtomRef(TARGET, n), n[0], p) for n, p in points.items()), (missing,))
    return request, truth


@pytest.mark.parametrize("name,missing", (("VAL", "CG1"), ("VAL", "CG2"), ("THR", "OG1"), ("THR", "CG2")))
@pytest.mark.parametrize("chi", (-73.2, 37.2, 140.1, 178.9))
def test_one_tip_is_recovered_from_fixed_partner_without_chi_grid(name, missing, chi):
    req, truth = anchored_request(name, missing, chi)
    before = req.digest
    result = complete_anchored_branch(req)
    assert len(result.attempts) == 1 and req.digest == before
    assert result.revision == "canonical-val-thr-observed-branch-anchor-v1"
    assert result.model_digest is None and "two-angle" in result.sampling_policy
    candidate = result.attempts[0].candidate
    assert len(candidate.additions) == 1 and candidate.additions[0].atom == ExactAtomRef(TARGET, missing)
    assert math.dist(candidate.additions[0].xyz, truth) < 1e-12
    audit = audit_anchored_branch(req, candidate)
    assert audit == result.attempts[0].audit
    if audit.accepted:
        assert result.status == "completed" and result.selected == candidate
    else:
        # A cone solution can be stereochemically correct but clash with the
        # retained backbone. Never loosen physical gates to make a math test pass.
        assert result.status == "rejected" and result.selected is None
        assert all(f.startswith("generated_clash:") for f in audit.findings)
    assert complete_sidechain(req).status == "unsupported"
    assert not audit_candidate(req, candidate).accepted


@pytest.mark.parametrize("name,missing", (("VAL", "CG1"), ("VAL", "CG2"), ("THR", "OG1"), ("THR", "CG2")))
def test_opposite_cone_solution_is_rejected_even_with_matching_measured_chi(name, missing):
    req, _ = anchored_request(name, missing)
    candidate = complete_anchored_branch(req).attempts[0].candidate
    fixed = {o.atom.name: o.xyz for o in req.observations}
    branch = BRANCHES[name]
    partner = branch.second if missing == branch.first else branch.first
    cb = fixed["CB"]
    normal = unit(cross(sub(fixed["CA"], cb), sub(fixed[partner], cb)))
    point = candidate.additions[0].xyz
    wrong = add(point, scaled(normal, -2*dot(sub(point, cb), normal)))
    points = fixed | {missing: wrong}
    chi = signed_dihedral_degrees(*(points[n] for n in ("N", "CA", "CB", branch.first)))
    forged = replace(candidate, additions=(Addition(ExactAtomRef(TARGET, missing), missing[0], wrong),), torsion_degrees=chi)
    audit = audit_anchored_branch(req, forged)
    expected = "thr_cb_not_R_or_nearplanar" if name == "THR" else "val_branch_name_orientation"
    assert not audit.accepted and expected in audit.findings
    assert "torsion_metadata_or_grid_mismatch" not in audit.findings


@pytest.mark.parametrize("angle", (0., 175.))
def test_bad_observed_anchor_geometry_refuses_without_moving_anchors(angle):
    req, _ = anchored_request("THR", "CG2")
    fixed = {o.atom.name: o.xyz for o in req.observations}
    axis = unit(sub(fixed["CA"], fixed["CB"]))
    perpendicular = unit(cross(axis, (1., 0., 0.)))
    vector = add(scaled(axis, 1.43*math.cos(math.radians(angle))),
                 scaled(perpendicular, 1.43*math.sin(math.radians(angle))))
    changed = replace(req, observations=tuple(replace(o, xyz=add(fixed["CB"], vector))
                                              if o.atom.name == "OG1" else o for o in req.observations))
    before = changed.digest
    result = complete_anchored_branch(changed)
    assert result.status == "rejected" and "observed_anchor_geometry:" in result.reason
    assert not result.attempts and changed.digest == before


def test_failed_clash_attempt_is_retained_not_optimized_or_resampled():
    req, _ = anchored_request("THR", "OG1")
    tip = complete_anchored_branch(req).selected.additions[0].xyz
    blocker = Observation(ExactAtomRef(ComponentInstanceRef(1, "W", "1"), "O"), "O", tip)
    changed = replace(req, observations=req.observations+(blocker,))
    result = complete_anchored_branch(changed)
    assert result.status == "rejected" and len(result.attempts) == 1 and result.selected is None
    assert any(f.startswith("generated_clash:") for f in result.attempts[0].audit.findings)


@pytest.mark.parametrize("tamper", (
    lambda c: replace(c, additions=()),
    lambda c: replace(c, request_digest="0"*64),
    lambda c: replace(c, torsion_degrees=math.nan),
    lambda c: replace(c, torsion_degrees=999.),
    lambda c: replace(c, torsion_degrees=True),
    lambda c: replace(c, torsion_degrees=c.torsion_degrees+1),
))
def test_anchored_audit_does_not_trust_provenance_or_measured_chi(tamper):
    req, _ = anchored_request("THR", "CG2")
    candidate = complete_anchored_branch(req).selected
    assert not audit_anchored_branch(req, tamper(candidate)).accepted


def test_paired_source_mode_is_not_implicitly_expanded_to_singletons():
    req, _ = anchored_request("THR", "OG1")
    lines = []
    for i, o in enumerate(req.observations, 1):
        line = _format_generated_atom(serial=i, name=o.atom.name, element=o.element, position=o.xyz,
                                      residue_name="THR", chain_id=TARGET.chain_id,
                                      sequence_number=TARGET.sequence_number, insertion_code=TARGET.insertion_code)
        lines.append("ATOM  "+line[6:])
    with pytest.raises(ValueError, match="partial_beta_branch"):
        plan_branched_source(("".join(lines)+"TER\nEND\n").encode(), TARGET)


@pytest.mark.parametrize("mask", ((), ("CB",), ("OG1", "CG2")))
def test_anchored_mode_requires_exactly_one_distal_atom(mask):
    req, _ = anchored_request("THR", "OG1")
    result = complete_anchored_branch(replace(req, missing_atoms=mask))
    assert result.status == "unsupported" and result.reason == "anchored_branch_requires_one_distal_atom"


@pytest.mark.parametrize("orientation", (0, True, 2))
def test_cone_primitive_requires_an_explicit_valid_hemisphere(orientation):
    with pytest.raises(ValueError):
        two_angle_tip((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), length=1.5,
                      first_angle=109., second_angle=109., orientation=orientation)
