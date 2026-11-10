"""Explicit paired-branch PDB admission, global THR-CB and rounded readback."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.branched_source import (
    audit_branched_materialized_source,
    materialize_branched_source,
    plan_branched_source,
)
from dvbfixer.protein_completion.cohort import plan_cohort
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import add, cross, dot, scaled, sub, unit
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import plan_source

TARGET = ComponentInstanceRef(1, "a", "82", "A")
BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def atom(serial, name, xyz, residue="THR", chain="a", number="82", code="A"):
    line = _format_generated_atom(serial=serial, name=name, element=name[0], position=xyz,
                                  residue_name=residue, chain_id=chain, sequence_number=number, insertion_code=code)
    return ("ATOM  "+line[6:]).encode()


def source(name="THR", with_cb=False):
    points = dict(BACKBONE)
    if with_cb:
        from dvbfixer.protein_completion.geometry import l_cb
        points["CB"] = l_cb(points["N"], points["CA"], points["C"])
    return b"HEADER    EXPLICIT BRANCHED SOURCE\n"+b"".join(
        atom(i, n, p, residue=name) for i, (n, p) in enumerate(points.items(), 1))+b"TER\nEND\n"


def full_thr_context(*, mirrored=False, planar=False, missing=None, chain="B", code="", number="1"):
    target = ComponentInstanceRef(1, chain, number, code)
    req = CompletionRequest("a"*64, target, "THR",
        tuple(Observation(ExactAtomRef(target, n), n[0], p) for n, p in BACKBONE.items()), ("CB", "OG1", "CG2"))
    candidate = complete_sidechain(req).selected
    points = BACKBONE | {a.atom.name: a.xyz for a in candidate.additions}
    if mirrored or planar:
        cb = points["CB"]
        normal = unit(cross(sub(points["CA"], cb), sub(points["OG1"], cb)))
        factor = -2 if mirrored else -1
        points["CG2"] = add(points["CG2"], scaled(normal, factor*dot(sub(points["CG2"], cb), normal)))
    return b"".join(atom(i+20, n, tuple(x+30 for x in p), chain=chain, number=number, code=code)
                    for i, (n, p) in enumerate(points.items(), 1) if n != missing)+b"TER\n"


@pytest.mark.parametrize("name", ("VAL", "THR"))
@pytest.mark.parametrize("with_cb", (False, True))
def test_explicit_source_has_exact_ownership_and_roundtrip(name, with_cb, tmp_path):
    src = source(name, with_cb)
    plan = plan_branched_source(src, TARGET)
    result = complete_sidechain(plan.request)
    out = materialize_branched_source(src, plan, result.selected)
    assert out.audit == audit_branched_materialized_source(src, plan, result.selected, out.pdb)
    assert {a.atom.name for a in out.audit.generated} == set(plan.request.missing_atoms)
    assert out.audit.plan_digest == plan.digest
    if name == "THR":
        assert plan.source_thr_beta[0].state == "target_pending"
        assert all(c.state == "observed_R" and c.volume_angstrom3 > .5 for c in out.audit.thr_beta)
    record = out.dat_record()
    path = tmp_path/"result.dat"
    record.save(path)
    assert record.load(path).added_atoms == record.added_atoms
    assert all(a["chain"] == "a" and a["icode"] == "A" and a["resname"] == name for a in record.added_atoms)
    for planner, target in ((plan_source, TARGET), (plan_cohort, (TARGET,))):
        with pytest.raises(ValueError, match="target capability remains"):
            planner(src, target)


@pytest.mark.parametrize("defect", ("mirrored", "planar"))
def test_global_thr_cb_gate_rejects_wrong_context_despite_correct_alpha(defect):
    src = source("VAL").replace(b"END\n", full_thr_context(**{defect: True})+b"END\n")
    with pytest.raises(ValueError, match="non-R/nearplanar global THR CB"):
        plan_branched_source(src, TARGET)


@pytest.mark.parametrize("missing", ("OG1", "CG2"))
def test_unselected_incomplete_thr_is_not_treated_as_a_checked_stereocentre(missing):
    src = source("VAL").replace(b"END\n", full_thr_context(missing=missing)+b"END\n")
    with pytest.raises(ValueError, match="undefined global THR CB"):
        plan_branched_source(src, TARGET)


@pytest.mark.parametrize("chain,code", (("A", "A"), ("a", "B")))
def test_chain_case_and_icode_do_not_extend_pending_target_exemption(chain, code):
    context = full_thr_context(missing="CG2", chain=chain, number="82", code=code)
    src = source().replace(b"END\n", context+b"END\n")
    with pytest.raises(ValueError, match="undefined global THR CB"):
        plan_branched_source(src, TARGET)


@pytest.mark.parametrize("tamper", (
    lambda p: replace(p, request=replace(p.request, source_sha256="0"*64)),
    lambda p: replace(p, request=replace(p.request, observations=p.request.observations[:-1])),
    lambda p: replace(p, source_thr_beta=()),
))
def test_full_request_context_and_beta_audits_are_rederived(tamper):
    src = source()
    plan = plan_branched_source(src, TARGET)
    candidate = complete_sidechain(plan.request).selected
    with pytest.raises(ValueError, match="binding"):
        materialize_branched_source(src, tamper(plan), candidate)
    with pytest.raises(ValueError, match="immutable branched-source"):
        replace(plan, policy="unqualified-new-policy")


@pytest.mark.parametrize("name", ("VAL", "THR"))
def test_rounding_clash_refuses_before_returning_an_accepted_structure(name):
    water = atom(10, "O", (.589, -.6, 1.405), residue="HOH", chain="W", number="1", code="")
    src = source(name).replace(b"END\n", water+b"END\n")
    plan = plan_branched_source(src, TARGET)
    candidate = complete_sidechain(plan.request).selected
    assert candidate is not None
    with pytest.raises(ValueError, match="rounded branched PDB rejected"):
        materialize_branched_source(src, plan, candidate)


def test_independent_readback_rejects_changed_source_records_and_new_atom_shifts():
    src = source()
    plan = plan_branched_source(src, TARGET)
    chosen = complete_sidechain(plan.request).selected
    out = materialize_branched_source(src, plan, chosen).pdb
    for bad in (out.replace(b"  -1.460", b"  -1.461", 1), out.replace(b"   0.536", b"   0.537", 1),
                out.replace(b"END\n", b"REMARK injected\nEND\n")):
        assert bad != out
        with pytest.raises(ValueError):
            audit_branched_materialized_source(src, plan, chosen, bad)


def test_explicit_branch_graph_support_and_noncanonical_edge_refusal():
    src = source(with_cb=True).replace(b"END\n", b"CONECT    2    5\nEND\n")
    plan = plan_branched_source(src, TARGET)
    out = materialize_branched_source(src, plan, complete_sidechain(plan.request).selected)
    assert b"CONECT    2    5\n" in out.pdb
    wrong = src.replace(b"CONECT    2    5", b"CONECT    4    5")
    with pytest.raises(ValueError, match="noncanonical"):
        plan_branched_source(wrong, TARGET)


def test_anisou_and_sparse_ter_serials_survive_branched_readback():
    src = source()
    oxygen = next(l for l in src.splitlines(keepends=True) if l.startswith(b"ATOM  ") and l[12:16].strip() == b"O")
    anisou = b"ANISOU"+oxygen[6:28]+b"   1000   1000   1000      0      0      0"+oxygen[70:]
    src = src.replace(oxygen, oxygen+anisou).replace(b"TER\n", b"TER   99996\n")
    plan = plan_branched_source(src, TARGET)
    out = materialize_branched_source(src, plan, complete_sidechain(plan.request).selected)
    assert anisou in out.pdb and b"TER   99996\n" in out.pdb
    assert all(f"ATOM  {serial:5d}".encode() in out.pdb for serial in (99997, 99998, 99999))


def test_serial_exhaustion_does_not_return_a_partial_branched_output():
    src = source().replace(b"TER\n", b"TER   99997\n")
    plan = plan_branched_source(src, TARGET)
    chosen = complete_sidechain(plan.request).selected
    assert chosen is not None and len(chosen.additions) == 3
    with pytest.raises(ValueError, match="limits"):
        materialize_branched_source(src, plan, chosen)


@pytest.mark.parametrize("resid,name", ((5, "VAL"), (26, "VAL"), (7, "THR"), (55, "THR")))
@pytest.mark.parametrize("with_cb", (False, True))
def test_real_ubiquitin_paired_branches_with_global_thr_stereo(resid, name, with_cb):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    missing = {b"CG1", b"CG2"} if name == "VAL" else {b"OG1", b"CG2"}
    if not with_cb:
        missing.add(b"CB")
    src = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72
                   and not (int(l[22:26]) == resid and l[12:16].strip() in missing))+b"TER\nEND\n"
    plan = plan_branched_source(src, ComponentInstanceRef(1, "A", str(resid)))
    out = materialize_branched_source(src, plan, complete_sidechain(plan.request).selected)
    assert len(out.audit.thr_beta) == 7 and all(c.state == "observed_R" for c in out.audit.thr_beta)
    assert set(a.atom.name.encode() for a in out.audit.generated) == missing
