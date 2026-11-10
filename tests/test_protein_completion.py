"""Private first-slice acceptance; no scientific toolkit or PDBFixer required."""
from __future__ import annotations

import hashlib
import math
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.diagnose.geometry import angle_degrees, signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import (
    Addition,
    CompletionRequest,
    Observation,
)
from dvbfixer.protein_completion.geometry import extend, l_cb
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate

COMPONENT = ComponentInstanceRef(1, "a", "82", "A", occurrence=2)
BACKBONE = {"N": (-1.46, 0.0, 0.0), "CA": (0.0, 0.0, 0.0),
            "C": (0.55, 1.43, 0.0), "O": (1.76, 1.66, 0.0)}
ROOT = Path(__file__).resolve().parents[1]


def observation(name, xyz, component=COMPONENT, element=None):
    return Observation(ExactAtomRef(component, name), element or name[0], tuple(xyz))


def request(name="ALA", points=None, missing=None, context=(), links=()):
    points = BACKBONE if points is None else points
    expected = {"ALA": {"CB"}, "SER": {"CB", "OG"}, "CYS": {"CB", "SG"}}
    missing = tuple(sorted(expected.get(name, set())-points.keys())) if missing is None else missing
    return CompletionRequest(hashlib.sha256(b"source-control").hexdigest(), COMPONENT, name,
                             tuple(observation(n, p) for n, p in points.items())+context, missing, links)


@pytest.mark.parametrize("name", ("ALA", "SER", "CYS"))
def test_whole_supported_sidechain_has_l_geometry_and_exact_ownership(name):
    req = request(name)
    before = req.digest
    result = complete_sidechain(req)
    assert result.status == "completed"
    assert req.digest == before
    assert {a.atom.name for a in result.selected.additions} == set(req.missing_atoms)
    assert all(a.atom.component == COMPONENT for a in result.selected.additions)
    assert all(a.audit.ca_volume_angstrom3 > 0.5 for a in result.attempts)
    assert audit_candidate(req, result.selected).accepted
    assert len(result.attempts) == (1 if name == "ALA" else 3)
    assert "no empirical rotamer probabilities" in result.sampling_policy


@pytest.mark.parametrize("name,terminal", (("SER", "OG"), ("CYS", "SG")))
def test_only_terminal_side_atom_is_generated_and_cb_is_unchanged(name, terminal):
    cb = l_cb(BACKBONE["N"], BACKBONE["CA"], BACKBONE["C"])
    points = BACKBONE | {"CB": cb}
    req = request(name, points)
    result = complete_sidechain(req)
    assert result.status == "completed"
    assert [a.atom.name for a in result.selected.additions] == [terminal]
    assert next(o for o in req.observations if o.atom.name == "CB").xyz == cb


@pytest.mark.parametrize("chi", (-60.0, 60.0, 180.0))
def test_internal_coordinates_have_declared_bond_angle_and_signed_chi(chi):
    n, ca, c = (BACKBONE[k] for k in ("N", "CA", "C"))
    cb = l_cb(n, ca, c)
    og = extend(n, ca, cb, 1.417, 110.8, chi)
    assert math.dist(cb, og) == pytest.approx(1.417)
    assert angle_degrees(ca, cb, og) == pytest.approx(110.8)
    observed = signed_dihedral_degrees(n, ca, cb, og)
    assert abs((observed-chi+180) % 360-180) < 1e-9


def test_noop_does_not_certify_invalid_existing_geometry():
    cb = l_cb(BACKBONE["N"], BACKBONE["CA"], BACKBONE["C"])
    assert complete_sidechain(request(points=BACKBONE | {"CB": cb})).status == "unchanged"
    bad = request(points=BACKBONE | {"CB": (cb[0], cb[1], -cb[2])})
    assert complete_sidechain(bad).status == "rejected"


@pytest.mark.parametrize("name", ("GLY", "ILE", "LYS", "PRO", "NAG", "DA", "NLN", "CYX"))
def test_other_residues_are_not_silently_routed_to_protein_builder(name):
    result = complete_sidechain(request(name, missing=()))
    assert result.status == "unsupported"
    assert result.reason == "unsupported_residue"
    assert result.selected is None and not result.attempts


@pytest.mark.parametrize("backbone", ("N", "CA", "C", "O"))
def test_backbone_defects_are_refused_before_generation(backbone):
    result = complete_sidechain(request(points={k: v for k, v in BACKBONE.items() if k != backbone},
                                       missing=("CB", backbone)))
    assert result.status == "unsupported" and result.reason == "missing_backbone_not_supported"


def test_retained_distal_atom_does_not_get_repositioned_to_hide_missing_cb():
    result = complete_sidechain(request("SER", BACKBONE | {"OG": (0.0, 0.0, 3.0)}))
    assert result.status == "unsupported"
    assert result.reason == "missing_cb_with_observed_distal_atoms_requires_closure"


def test_covalent_sulfur_link_is_not_treated_as_free_cysteine():
    partner = ComponentInstanceRef(1, "B", "100")
    outside = observation("SG", (0, 0, 5), partner)
    req = request("CYS", context=(outside,),
                  links=((ExactAtomRef(COMPONENT, "SG"), outside.atom),))
    result = complete_sidechain(req)
    assert result.status == "unsupported" and "link_aware" in result.reason


def test_explicit_peptide_link_is_supported_but_spurious_distant_bond_is_not():
    previous = observation("C", (-2.1, -1.15, 0), ComponentInstanceRef(1, "a", "81"))
    peptide = (previous.atom, ExactAtomRef(COMPONENT, "N"))
    req = request(context=(previous,), links=(peptide,))
    assert complete_sidechain(req).status == "completed"
    distant = replace(previous, xyz=(20, 0, 0))
    result = complete_sidechain(request(context=(distant,), links=(peptide,)))
    assert result.status == "unsupported" and result.reason == "invalid_peptide_link_length"


@pytest.mark.parametrize("element", ("ZN", "FE", "MG"))
@pytest.mark.parametrize("z", (3, -6))
def test_nearby_metal_is_not_mistaken_for_generic_nonbonded_environment(element, z):
    metal = observation(element, (0, 0, z), ComponentInstanceRef(1, "Z", "1"), element)
    result = complete_sidechain(request("CYS", context=(metal,)))
    assert result.status == "unsupported" and "coordination" in result.reason


def test_context_blocker_rejects_all_candidates_and_preserves_attempts():
    cb = complete_sidechain(request()).selected.additions[0].xyz
    blocker = observation("C1", cb, ComponentInstanceRef(1, "A", "82", "A", occurrence=2), "C")
    req = request(context=(blocker,))
    result = complete_sidechain(req)
    assert result.status == "rejected" and result.selected is None
    assert len(result.attempts) == 1
    assert any("generated_clash" in f for f in result.attempts[0].audit.findings)
    assert req.observations[-1] == blocker  # A versus a has not collapsed


def test_context_changes_ranking_without_reference_coordinates():
    original = complete_sidechain(request("SER"))
    blocked = original.selected
    og = next(a.xyz for a in blocked.additions if a.atom.name == "OG")
    cb = next(a.xyz for a in blocked.additions if a.atom.name == "CB")
    # Put the blocker beyond OG, not on top of the common Cβ anchor: the latter
    # correctly makes *all* orientations infeasible, rather than testing ranking.
    outward = tuple(o+1.3*(o-c)/math.dist(og, cb) for o, c in zip(og, cb))
    blocker = observation("O1", outward, ComponentInstanceRef(1, "L", "1"), "O")
    result = complete_sidechain(request("SER", context=(blocker,)))
    assert result.status == "completed"
    assert result.selected.torsion_degrees != blocked.torsion_degrees
    assert any(not attempt.audit.accepted for attempt in result.attempts)


@pytest.mark.parametrize("change,reason", (
    ({"missing_atoms": ("CB", "OXT")}, "missing_mask"),
    ({"missing_atoms": ()}, "missing_mask"),
    ({"residue": replace(COMPONENT, alternate_location="B")}, "alternate_location"),
))
def test_incompatible_masks_and_selections_refuse(change, reason):
    result = complete_sidechain(replace(request(), **change))
    assert result.status == "unsupported" and reason in result.reason


@pytest.mark.parametrize("mutation", ("digest", "extra", "missing", "duplicate", "element", "nan", "mirror", "planar"))
def test_independent_validator_rejects_forged_candidates(mutation):
    req = request()
    candidate = complete_sidechain(req).selected
    a = candidate.additions[0]
    if mutation == "digest":
        candidate = replace(candidate, request_digest="0"*64)
    elif mutation == "extra":
        candidate = replace(candidate, additions=(a, Addition(ExactAtomRef(COMPONENT, "CA"), "C", (0, 0, 0))))
    elif mutation == "missing":
        candidate = replace(candidate, additions=())
    elif mutation == "duplicate":
        candidate = replace(candidate, additions=(a, a))
    elif mutation == "element":
        candidate = replace(candidate, additions=(replace(a, element="O"),))
    elif mutation == "nan":
        candidate = replace(candidate, additions=(replace(a, xyz=(math.nan, 0, 0)),))
    else:
        xyz = (a.xyz[0], a.xyz[1], -a.xyz[2] if mutation == "mirror" else 0.0)
        candidate = replace(candidate, additions=(replace(a, xyz=xyz),))
    assert not audit_candidate(req, candidate).accepted


def test_candidate_is_bound_to_context_and_source_bytes():
    req = request()
    candidate = complete_sidechain(req).selected
    assert not audit_candidate(replace(req, source_sha256="0"*64), candidate).accepted
    changed = replace(req, observations=req.observations+(observation("C1", (20, 20, 20),
                                                                   ComponentInstanceRef(1, "B", "2")),))
    assert not audit_candidate(changed, candidate).accepted


def test_auditor_recomputes_chi_instead_of_trusting_candidate_annotation():
    req = request("SER")
    result = complete_sidechain(req)
    assert not audit_candidate(req, replace(result.selected, torsion_degrees=15)).accepted
    assert not audit_candidate(req, replace(result.selected, torsion_degrees=[])).accepted


def test_explicit_unsupported_context_and_model_do_not_disappear():
    req = request()
    for context in (
        observation("XE", (10, 10, 10), ComponentInstanceRef(1, "B", "1"), "XE"),
        observation("C1", (10, 10, 10), ComponentInstanceRef(2, "B", "1"), "C"),
        observation("C1", (10, 10, 10), ComponentInstanceRef(1, "B", "1", alternate_location="A"), "C"),
    ):
        assert complete_sidechain(replace(req, observations=req.observations+(context,))).status == "unsupported"


@pytest.mark.parametrize("kind", ("nan", "duplicate", "mutable", "occupancy", "link"))
def test_contract_rejects_malformed_inputs(kind):
    with pytest.raises(ValueError):
        if kind == "nan":
            observation("CB", (0, math.nan, 0))
        elif kind == "duplicate":
            req = request()
            replace(req, observations=req.observations+(req.observations[0],))
        elif kind == "mutable":
            replace(request(), missing_atoms=["CB"])
        elif kind == "occupancy":
            replace(request().observations[0], occupancy=0)
        else:
            a = ExactAtomRef(COMPONENT, "CA")
            replace(request(), links=((a, a),))


def test_partial_occupancy_and_measured_hydrogen_are_not_silently_discarded():
    req = request()
    altered = replace(req.observations[0], occupancy=0.5)
    assert complete_sidechain(replace(req, observations=(altered,)+req.observations[1:])).status == "unsupported"
    hydrogen = observation("H", (-2, 0, 0))
    assert complete_sidechain(replace(req, observations=req.observations+(hydrogen,))).status == "unsupported"


def test_degenerate_anchors_reject_without_publishing_additions():
    points = BACKBONE | {"N": (0, 0, 0)}
    result = complete_sidechain(request(points=points))
    assert result.status == "rejected" and result.selected is None


def test_selected_source_d_ca_is_not_reflected_during_terminal_completion():
    cb = l_cb(BACKBONE["N"], BACKBONE["CA"], BACKBONE["C"])
    points = BACKBONE | {"CB": (cb[0], cb[1], -cb[2])}
    result = complete_sidechain(request("SER", points))
    assert result.status == "rejected" and result.selected is None
    assert all("non_L_or_nearplanar_ca" in a.audit.findings for a in result.attempts)


def test_rigid_motion_and_observation_order_do_not_change_proposals():
    req = request("SER")
    baseline = complete_sidechain(req)
    def motion(p):
        return -p[1]+5, p[0]-7, p[2]+3
    moved = replace(req, observations=tuple(replace(o, xyz=motion(o.xyz)) for o in reversed(req.observations)))
    actual = complete_sidechain(moved)
    assert actual.status == baseline.status
    for a, b in zip(actual.attempts, baseline.attempts):
        for x, y in zip(a.candidate.additions, b.candidate.additions):
            assert x.xyz == pytest.approx(motion(y.xyz), abs=1e-12)
        assert a.audit.accepted == b.audit.accepted


@pytest.mark.parametrize("residue", ("ALA", "SER", "CYS"))
def test_one_hundred_proper_rotations_preserve_all_candidate_geometry(residue):
    req = request(residue)
    reference = complete_sidechain(req)
    # Rodrigues rotations about an oblique axis, including translations. This
    # exercises the actual private service, not just the earlier research helper.
    axis = (2/3, -1/3, 2/3)
    def motion(p, theta):
        x, y, z = p
        a, b, c = axis
        perpendicular = (b*z-c*y, c*x-a*z, a*y-b*x)
        projection = sum(u*v for u, v in zip(axis, p))
        return tuple(u*math.cos(theta)+v*math.sin(theta)+w*projection*(1-math.cos(theta))+t
                     for u, v, w, t in zip(p, perpendicular, axis, (17, -8, 21)))
    for i in range(100):
        theta = i*0.731
        rotated = replace(req, observations=tuple(replace(o, xyz=motion(o.xyz, theta)) for o in req.observations))
        result = complete_sidechain(rotated)
        assert result.status == reference.status
        for actual, expected in zip(result.attempts, reference.attempts):
            assert actual.audit.accepted == expected.audit.accepted
            assert actual.audit.ca_volume_angstrom3 == pytest.approx(expected.audit.ca_volume_angstrom3)
            for generated, original in zip(actual.candidate.additions, expected.candidate.additions):
                assert generated.xyz == pytest.approx(motion(original.xyz, theta), abs=1e-10)


def test_request_identity_model_case_icode_and_occurrence_survive():
    req = request()
    result = complete_sidechain(req)
    assert result.selected.additions[0].atom == ExactAtomRef(COMPONENT, "CB")
    for field, value in (("model", 2), ("chain_id", "A"), ("insertion_code", "a"), ("occurrence", 1)):
        different = replace(COMPONENT, **{field: value})
        other = replace(req, residue=different,
                        observations=tuple(replace(o, atom=ExactAtomRef(different, o.atom.name)) for o in req.observations))
        assert other.digest != req.digest


def test_private_service_imports_without_scientific_stack():
    code = """
import sys
sys.path.insert(0, 'src')
from dvbfixer.protein_completion.service import complete_sidechain
assert not {'numpy', 'openmm', 'pdbfixer', 'Bio', 'rdkit'} & sys.modules.keys()
"""
    subprocess.run([sys.executable, "-S", "-c", code], cwd=ROOT, check=True)


@pytest.mark.parametrize("resid,name,missing", ((28, "ALA", ("CB",)), (20, "SER", ("OG",)),
                                              (20, "SER", ("CB", "OG"))))
def test_withheld_real_fixture_atoms_are_not_proposal_inputs(resid, name, missing):
    source = (ROOT/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    target_component = ComponentInstanceRef(1, "A", str(resid))
    points = {}
    for line in source.decode().splitlines():
        if line.startswith("ATOM  ") and line[21] == "A" and int(line[22:26]) == resid:
            points[line[12:16].strip()] = tuple(float(line[i:i+8]) for i in (30, 38, 46))
    observed = tuple(observation(n, p, target_component) for n, p in points.items() if n not in missing)
    req = CompletionRequest(hashlib.sha256(source).hexdigest(), target_component, name, observed, missing)
    result = complete_sidechain(req)
    assert result.status == "completed"
    # Deposited coordinates are consulted only after candidate generation.
    error = math.sqrt(sum(math.dist(a.xyz, points[a.atom.name])**2
                          for a in result.selected.additions)/len(missing))
    assert math.isfinite(error)  # no reference-fit gate or oracle candidate selection
    assert all(a.atom.name in missing for a in result.selected.additions)
