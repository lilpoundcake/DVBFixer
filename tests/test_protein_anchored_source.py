"""Fixed-partner source/policy binding and all-or-nothing singleton readback."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import (
    audit_anchored_materialized_source,
    materialize_anchored_source,
    plan_anchored_source,
)
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.cohort import plan_cohort
from dvbfixer.protein_completion.geometry import add, cross, dot, extend, l_cb, scaled, sub, unit
from dvbfixer.protein_completion.source import plan_source

TARGET = ComponentInstanceRef(1, "a", "82", "A")
BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}


def full_points(name):
    points = dict(BACKBONE)
    cb = points["CB"] = l_cb(points["N"], points["CA"], points["C"])
    branch = BRANCHES[name]
    points[branch.first] = extend(points["N"], points["CA"], cb, branch.first_length, TETRAHEDRAL_ANGLE, 60.)
    points[branch.second] = extend(points["N"], points["CA"], cb, branch.second_length,
                                   TETRAHEDRAL_ANGLE, 60.+branch.second_chi_offset)
    return points


def atom(serial, name, xyz, residue="THR", chain="a", number="82", code="A"):
    line = _format_generated_atom(serial=serial, name=name, element=name[0], position=xyz,
                                  residue_name=residue, chain_id=chain, sequence_number=number, insertion_code=code)
    return ("ATOM  "+line[6:]).encode()


def source(name="THR", missing="OG1"):
    return b"HEADER    OBSERVED PARTNER CONTROL\n"+b"".join(
        atom(i, n, p, residue=name) for i, (n, p) in enumerate(full_points(name).items(), 1)
        if n != missing)+b"TER\nEND\n"


@pytest.mark.parametrize("name,missing", (("VAL", "CG1"), ("VAL", "CG2"), ("THR", "OG1"), ("THR", "CG2")))
def test_singleton_roundtrip_and_exact_partner_ownership(name, missing, tmp_path):
    src = source(name, missing)
    plan = plan_anchored_source(src, TARGET)
    chosen = complete_anchored_branch(plan.request).selected
    assert chosen is not None
    out = materialize_anchored_source(src, plan, chosen)
    audit = audit_anchored_materialized_source(src, plan, chosen, out.pdb)
    assert out.audit == audit and audit.plan_digest == plan.digest
    assert len(audit.generated) == 1 and audit.generated[0].atom.name == missing
    if name == "THR":
        assert plan.source_thr_beta[0].state == "target_pending"
        assert all(c.state == "observed_R" for c in audit.thr_beta)
    partner = plan.fixed_partner
    original_partner = next(l for l in src.splitlines(keepends=True) if l.startswith(b"ATOM  ")
                            and l[12:16].strip().decode() == partner.atom.name)
    assert original_partner in out.pdb
    record = out.dat_record()
    record.save(tmp_path/"one.dat")
    assert record.load(tmp_path/"one.dat").to_dict() == record.to_dict()
    assert len(record.added_atoms) == 1 and record.added_atoms[0]["atom"] == missing
    assert record.added_atoms[0]["icode"] == "A" and record.added_atoms[0]["chain"] == "a"
    assert partner.atom.name not in {a["atom"] for a in record.added_atoms}
    for planner, target in ((plan_source, TARGET), (plan_cohort, (TARGET,)), (plan_branched_source, TARGET)):
        with pytest.raises(ValueError):
            planner(src, target)


def test_partner_changed_in_both_claimed_observations_and_metadata_still_refuses():
    src = source()
    plan = plan_anchored_source(src, TARGET)
    chosen = complete_anchored_branch(plan.request).selected
    partner = plan.fixed_partner
    changed = replace(partner, xyz=(partner.xyz[0]+.001, *partner.xyz[1:]))
    request = replace(plan.request, observations=tuple(changed if o.atom == partner.atom else o
                                                       for o in plan.request.observations))
    forged = replace(plan, request=request, fixed_partner=changed)
    with pytest.raises(ValueError, match="binding"):
        materialize_anchored_source(src, forged, chosen)
    with pytest.raises(ValueError, match="distal partner"):
        replace(plan, fixed_partner=next(o for o in plan.request.observations if o.atom.name == "CA"))


@pytest.mark.parametrize("tamper", (
    lambda p: replace(p, request=replace(p.request, source_sha256="0"*64)),
    lambda p: replace(p, source_thr_beta=()),
    lambda p: replace(p, policy="paired-val-thr-source-global-thr-beta-v1"),
))
def test_source_stereo_and_policy_cannot_be_substituted(tamper):
    src = source()
    plan = plan_anchored_source(src, TARGET)
    chosen = complete_anchored_branch(plan.request).selected
    with pytest.raises(ValueError):
        materialize_anchored_source(src, tamper(plan), chosen)


@pytest.mark.parametrize("chain,code,defect", (("A", "A", "missing"), ("a", "B", "missing"),
                                            ("B", "", "mirrored"), ("B", "", "planar")))
def test_pending_target_exemption_never_hides_invalid_context_thr(chain, code, defect):
    points = full_points("THR")
    if defect == "missing":
        del points["OG1"]
    else:
        cb = points["CB"]
        normal = unit(cross(sub(points["CA"], cb), sub(points["OG1"], cb)))
        factor = -2 if defect == "mirrored" else -1
        points["CG2"] = add(points["CG2"], scaled(normal, factor*dot(sub(points["CG2"], cb), normal)))
    extra = b"".join(atom(i+20, n, tuple(x+30 for x in p), chain=chain, number="82", code=code)
                     for i, (n, p) in enumerate(points.items(), 1))+b"TER\n"
    with pytest.raises(ValueError, match="global THR CB"):
        plan_anchored_source(source().replace(b"END\n", extra+b"END\n"), TARGET)


def test_generated_atom_and_partner_edits_are_independently_rejected():
    src = source()
    plan = plan_anchored_source(src, TARGET)
    chosen = complete_anchored_branch(plan.request).selected
    out = materialize_anchored_source(src, plan, chosen).pdb
    atoms = _parse_atom_lines(out.decode().splitlines(keepends=True))
    for name in ("OG1", plan.fixed_partner.atom.name):
        target = next(a for a in atoms if a.name == name)
        lines = out.decode().splitlines(keepends=True)
        line = lines[target.line_index]
        lines[target.line_index] = line[:30]+f"{target.position[0]+.001:8.3f}"+line[38:]
        with pytest.raises(ValueError):
            audit_anchored_materialized_source(src, plan, chosen, "".join(lines).encode())


def test_anchored_admission_does_not_exempt_an_observed_d_ca():
    lines = source().decode().splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith("ATOM  ") and line[12:16].strip() == "CB":
            lines[i] = line[:46]+f"{-float(line[46:54]):8.3f}"+line[54:]
    with pytest.raises(ValueError, match="non-L/nearplanar global protein CA"):
        plan_anchored_source("".join(lines).encode(), TARGET)


def test_source_admits_only_a_single_distal_atom_not_a_paired_request():
    src = source().replace(next(l for l in source().splitlines(keepends=True)
                               if l.startswith(b"ATOM  ") and l[12:16].strip() == b"CG2"), b"")
    paired = plan_branched_source(src, TARGET)
    with pytest.raises(ValueError, match="one_distal_atom"):
        plan_anchored_source(src, TARGET)
    with pytest.raises(ValueError, match="paired plans cannot opt in"):
        materialize_anchored_source(src, paired, None)


@pytest.mark.parametrize("missing,water", (("OG1", (-2.169, -1.093, -2.192)),
                                          ("CG2", (-.026, -2.343, 1.462))))
def test_single_atom_rounding_clash_rejects_without_returning_pdb(missing, water):
    src = source(missing=missing)
    water_line = atom(20, "O", water, residue="HOH", chain="W", number="1", code="")
    src = src.replace(b"END\n", water_line+b"END\n")
    plan = plan_anchored_source(src, TARGET)
    chosen = complete_anchored_branch(plan.request).selected
    assert chosen is not None
    with pytest.raises(ValueError, match="rounded anchored PDB rejected"):
        materialize_anchored_source(src, plan, chosen)


def test_partner_anisou_sparse_serials_and_conect_stay_exact():
    src = source()
    plan = plan_anchored_source(src, TARGET)
    partner = next(l for l in src.splitlines(keepends=True) if l.startswith(b"ATOM  ")
                   and l[12:16].strip().decode() == plan.fixed_partner.atom.name)
    anisou = b"ANISOU"+partner[6:28]+b"   1000   1000   1000      0      0      0"+partner[70:]
    src = src.replace(partner, partner+anisou).replace(b"TER\n", b"TER   99998\n")
    src = src.replace(b"END\n", b"CONECT    5    7\nEND\n")
    plan = plan_anchored_source(src, TARGET)
    out = materialize_anchored_source(src, plan, complete_anchored_branch(plan.request).selected)
    assert partner+anisou in out.pdb and b"CONECT    5    7\n" in out.pdb
    assert b"ATOM  99999" in out.pdb
    exhausted = src.replace(b"99998", b"99999")
    plan = plan_anchored_source(exhausted, TARGET)
    with pytest.raises(ValueError, match="limits"):
        materialize_anchored_source(exhausted, plan, complete_anchored_branch(plan.request).selected)


@pytest.mark.parametrize("number,name,missing", ((5, "VAL", "CG1"), (5, "VAL", "CG2"),
    (26, "VAL", "CG1"), (26, "VAL", "CG2"), (7, "THR", "OG1"), (7, "THR", "CG2"),
    (55, "THR", "OG1"), (55, "THR", "CG2")))
def test_real_ubiquitin_singleton_source_and_readback(number, name, missing):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    src = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72
                   and not (int(l[22:26]) == number and l[12:16].strip().decode() == missing))+b"TER\nEND\n"
    target = ComponentInstanceRef(1, "A", str(number))
    plan = plan_anchored_source(src, target)
    assert plan.request.residue_name == name
    chosen = complete_anchored_branch(plan.request).selected
    out = materialize_anchored_source(src, plan, chosen)
    assert len(out.audit.generated) == 1 and len(out.audit.thr_beta) == 7
    assert all(c.state == "observed_R" for c in out.audit.thr_beta)
