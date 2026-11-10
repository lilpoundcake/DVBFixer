"""Private VAL/THR paired branches; no empirical distribution or public capability."""
from __future__ import annotations

import hashlib
from dataclasses import replace

import pytest

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.branched import BRANCHES, named_beta_volume
from dvbfixer.protein_completion.cohort import plan_cohort
from dvbfixer.protein_completion.contracts import Candidate, CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import add, cross, dot, l_cb, scaled, sub, unit
from dvbfixer.protein_completion.rotamers import LibraryAuthority, parse_single_chi
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import plan_source
from dvbfixer.protein_completion.validation import audit_candidate

TARGET = ComponentInstanceRef(1, "a", "82", "A")
BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def request(name, with_cb=False):
    points = dict(BACKBONE)
    if with_cb:
        points["CB"] = l_cb(points["N"], points["CA"], points["C"])
    branch = BRANCHES[name]
    missing = (branch.first, branch.second) if with_cb else ("CB", branch.first, branch.second)
    return CompletionRequest("a"*64, TARGET, name,
                             tuple(Observation(ExactAtomRef(TARGET, n), n[0], p) for n, p in points.items()), missing)


@pytest.mark.parametrize("name", ("VAL", "THR"))
@pytest.mark.parametrize("with_cb", (False, True))
def test_paired_branches_keep_observations_and_named_stereo(name, with_cb):
    req = request(name, with_cb)
    before = req.digest
    result = complete_sidechain(req)
    assert result.status == "completed" and len(result.attempts) == 3 and req.digest == before
    assert result.revision == "canonical-branched-sidechain-staggered-v1"
    assert result.model_digest is None and "no empirical" in result.sampling_policy
    branch = BRANCHES[name]
    for attempt in result.attempts:
        assert attempt.audit.accepted
        candidate = attempt.candidate
        assert {a.atom.name for a in candidate.additions} == set(req.missing_atoms)
        points = {o.atom.name: o.xyz for o in req.observations} | {a.atom.name: a.xyz for a in candidate.additions}
        volume = named_beta_volume(points["CA"], points["CB"], points[branch.first], points[branch.second])
        assert branch.expected_named_volume_sign*volume > .5
        chi = signed_dihedral_degrees(*(points[n] for n in ("N", "CA", "CB", branch.first)))
        assert abs((chi-candidate.torsion_degrees+180) % 360-180) < 1e-9
        assert audit_candidate(req, candidate).accepted


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_reflected_second_branch_fails_without_changing_ca_or_primary_chi(name):
    req = request(name)
    candidate = complete_sidechain(req).attempts[0].candidate
    points = {o.atom.name: o.xyz for o in req.observations} | {a.atom.name: a.xyz for a in candidate.additions}
    branch = BRANCHES[name]
    cb = points["CB"]
    normal = unit(cross(sub(points["CA"], cb), sub(points[branch.first], cb)))
    second = points[branch.second]
    reflected = add(second, scaled(normal, -2*dot(sub(second, cb), normal)))
    additions = tuple(replace(a, xyz=reflected) if a.atom.name == branch.second else a for a in candidate.additions)
    audit = audit_candidate(req, replace(candidate, additions=additions))
    finding = "thr_cb_not_R_or_nearplanar" if name == "THR" else "val_branch_name_orientation"
    assert not audit.accepted and finding in audit.findings
    assert "non_L_or_nearplanar_ca" not in audit.findings
    assert "torsion_metadata_or_grid_mismatch" not in audit.findings
    assert not any(f.startswith(("bond_length:", "bond_angle:")) for f in audit.findings)
    # Fixed observations are not silently reflected or normalized on a no-op.
    fixed = tuple(Observation(a.atom, a.element, a.xyz) for a in additions)
    complete = replace(req, observations=req.observations+fixed, missing_atoms=())
    before = complete.digest
    result = complete_sidechain(complete)
    assert result.status == "rejected" and finding in result.reason and complete.digest == before


@pytest.mark.parametrize("name,observed", (("VAL", "CG1"), ("VAL", "CG2"), ("THR", "OG1"), ("THR", "CG2")))
def test_single_missing_branch_refuses_instead_of_overwriting_fixed_partner(name, observed):
    req = request(name, with_cb=True)
    addition = next(a for a in complete_sidechain(req).selected.additions if a.atom.name == observed)
    changed = replace(req, observations=req.observations+(Observation(addition.atom, addition.element, addition.xyz),),
                      missing_atoms=tuple(n for n in req.missing_atoms if n != observed))
    result = complete_sidechain(changed)
    assert result.status == "unsupported" and result.reason == "partial_beta_branch_requires_observed_anchor_solver"


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_missing_cb_with_fixed_distal_branch_still_refuses(name):
    req = request(name)
    candidate = complete_sidechain(req).selected
    kept = tuple(Observation(a.atom, a.element, a.xyz) for a in candidate.additions if a.atom.name != "CB")
    result = complete_sidechain(replace(req, observations=req.observations+kept, missing_atoms=("CB",)))
    assert result.status == "unsupported" and "requires_closure" in result.reason


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_complete_observations_noop_is_independently_audited(name):
    req = request(name)
    candidate = complete_sidechain(req).selected
    fixed = tuple(Observation(a.atom, a.element, a.xyz) for a in candidate.additions)
    complete = replace(req, observations=req.observations+fixed, missing_atoms=())
    assert complete_sidechain(complete).status == "unchanged"
    assert audit_candidate(complete, Candidate(complete.digest, (), None)).accepted


def test_both_branches_enter_clash_context_checks():
    req = request("THR")
    first = complete_sidechain(req).attempts[0].candidate
    cg2 = next(a.xyz for a in first.additions if a.atom.name == "CG2")
    blocker = Observation(ExactAtomRef(ComponentInstanceRef(1, "W", "1"), "O"), "O", cg2)
    changed = replace(req, observations=req.observations+(blocker,))
    candidate = replace(first, request_digest=changed.digest)
    audit = audit_candidate(changed, candidate)
    assert not audit.accepted and any(f.startswith("generated_clash:CG2:") for f in audit.findings)


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_new_capability_does_not_accept_a_missing_mask_that_omits_branches(name):
    result = complete_sidechain(replace(request(name), missing_atoms=()))
    assert result.status == "unsupported" and result.reason == "missing_mask_does_not_match_observations"


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_ser_authority_is_not_silently_used_for_branched_residues(name):
    data = b"".join(f"SER -60 -40 0 {i} 0 0 0 {p} {mean} 0 0 0 10 0 0 0\n".encode()
                    for i, p, mean in ((1, .5, 60), (2, .3, 180), (3, .2, -60)))
    authority = LibraryAuthority("synthetic-test", "SER", "authored numeric control", "synthetic-control", "test", "0"*64)
    library = parse_single_chi(data, expected_sha256=hashlib.sha256(data).hexdigest(), authority=authority)
    result = complete_sidechain(request(name), library=library)
    assert result.status == "unsupported" and result.reason == "branched_empirical_authority_not_supported"
    assert not result.attempts


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_typed_kernel_does_not_implicitly_expand_source_cohort_or_publication(name):
    lines = []
    for i, (n, p) in enumerate(BACKBONE.items(), 1):
        line = _format_generated_atom(serial=i, name=n, element=n[0], position=p, residue_name=name,
                                      chain_id=TARGET.chain_id, sequence_number=TARGET.sequence_number,
                                      insertion_code=TARGET.insertion_code)
        lines.append("ATOM  "+line[6:])
    src = ("".join(lines)+"TER\nEND\n").encode()
    assert complete_sidechain(request(name)).status == "completed"
    for planner, selection in ((plan_source, TARGET), (plan_cohort, (TARGET,))):
        with pytest.raises(ValueError, match="target capability remains"):
            planner(src, selection)
