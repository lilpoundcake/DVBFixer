"""New mixed policy: old contracts, partner ownership and global stereo stay strict."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.cohort import complete_cohort, plan_cohort
from dvbfixer.protein_completion.evidence import write_cohort_evidence
from dvbfixer.protein_completion.geometry import add, cross, dot, extend, l_cb, scaled, sub, unit
from dvbfixer.protein_completion.mixed_cohort import (
    MixedCohortPlan,
    audit_materialized_mixed_cohort,
    complete_mixed_cohort,
    materialize_mixed_cohort,
    plan_mixed_cohort,
)


def points(name):
    xyz = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    cb = xyz["CB"] = l_cb(xyz["N"], xyz["CA"], xyz["C"])
    if name in BRANCHES:
        branch = BRANCHES[name]
        xyz[branch.first] = extend(xyz["N"], xyz["CA"], cb, branch.first_length, TETRAHEDRAL_ANGLE, 60.)
        xyz[branch.second] = extend(xyz["N"], xyz["CA"], cb, branch.second_length,
                                    TETRAHEDRAL_ANGLE, 60.+branch.second_chi_offset)
    elif name in {"SER", "CYS"}:
        terminal, length, angle = ("OG", 1.417, 110.8) if name == "SER" else ("SG", 1.81, 114.4)
        xyz[terminal] = extend(xyz["N"], xyz["CA"], cb, length, angle, 60.)
    return xyz


def source(specs):
    """(name, chain, icode, missing, offset, rotate180x), all residue number 1."""
    lines = ["HEADER    MIXED CONTROL\n"]
    targets = []
    for j, (name, chain, code, missing, offset, rotate) in enumerate(specs):
        target = ComponentInstanceRef(1, chain, "1", code)
        targets.append(target)
        for i, (atom, xyz) in enumerate(points(name).items(), 1):
            if atom in missing:
                continue
            xyz = tuple(((-x if rotate and k in {1, 2} else x)+offset[k]) for k, x in enumerate(xyz))
            line = _format_generated_atom(serial=j*20+i, name=atom, element=atom[0], position=xyz,
                                          residue_name=name, chain_id=chain, sequence_number="1", insertion_code=code)
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    return ("".join(lines)+"END\n").encode(), tuple(targets)


FOUR = (("ALA", "A", "", ("CB",), (0., 0., 0.), False),
        ("THR", "a", "A", ("OG1", "CG2"), (30., 0., 0.), False),
        ("VAL", "B", "", ("CG1",), (60., 0., 0.), False),
        ("SER", "C", "", ("OG",), (90., 0., 0.), False))


def test_all_three_modes_joint_rounding_and_generated_only_dat(tmp_path):
    src, targets = source(FOUR)
    plan = plan_mixed_cohort(src, targets)
    assert [t.mode for t in plan.targets] == ["standard", "paired", "anchored", "standard"]
    result = complete_mixed_cohort(src, plan)
    assert result.selected is not None and len(result.attempts) == 9
    out = materialize_mixed_cohort(src, plan, result.selected)
    assert out.audit == audit_materialized_mixed_cohort(src, plan, result.selected, out.pdb)
    assert len(out.audit.generated) == 5 and all(a.accepted for a in out.audit.targets)
    assert all(c.state == "observed_R" for c in out.audit.thr_beta)
    record = out.dat_record()
    record.save(tmp_path/"mixed.dat")
    assert record.load(tmp_path/"mixed.dat").to_dict() == record.to_dict()
    anchored = plan.targets[2]
    assert anchored.fixed_partner.atom.name == "CG2"
    assert ("B", "1", "", "CG2") not in record.added_keys()
    assert ("B", "1", "", "CB") not in record.added_keys()
    assert record.added_keys() == {(a.atom.component.chain_id, "1", a.atom.component.insertion_code, a.atom.name)
                                  for a in out.audit.generated}
    with pytest.raises(ValueError):
        plan_cohort(src, targets)


def test_request_order_is_source_canonical_not_caller_or_case_sort():
    src, targets = source(FOUR)
    a, b = plan_mixed_cohort(src, targets), plan_mixed_cohort(src, tuple(reversed(targets)))
    assert a == b and a.digest == b.digest
    ar, br = complete_mixed_cohort(src, a), complete_mixed_cohort(src, b)
    assert ar == br
    assert materialize_mixed_cohort(src, a, ar.selected).pdb == materialize_mixed_cohort(src, b, br.selected).pdb


def test_paired_and_anchored_new_atoms_can_clash_despite_individual_passes():
    specs = (("THR", "A", "", ("CB", "OG1", "CG2"), (0., 0., 0.), False),
             ("THR", "a", "A", ("CG2",), (.8, -4.5, -4.95), True))
    src, targets = source(specs)
    plan = plan_mixed_cohort(src, targets)
    assert [t.mode for t in plan.targets] == ["paired", "anchored"]
    assert sum(c.state == "target_pending" for c in plan.source_thr_beta) == 2
    result = complete_mixed_cohort(src, plan)
    rejected = next(a for a in result.attempts if not a.accepted)
    assert all(result.targets[j].attempts[i].audit.accepted for j, i in enumerate(rejected.indices))
    assert any(f.startswith("generated_clash:") for a in rejected.audits for f in a.findings)
    candidates = tuple(result.targets[j].attempts[i].candidate for j, i in enumerate(rejected.indices))
    with pytest.raises(ValueError, match="joint mixed"):
        materialize_mixed_cohort(src, plan, candidates)


def test_two_singletons_with_no_joint_solution_return_no_partial_output():
    specs = (("THR", "A", "", ("CG2",), (0., 0., 0.), False),
             ("THR", "a", "A", ("CG2",), (.8, -4.5, -4.95), True))
    src, targets = source(specs)
    plan = plan_mixed_cohort(src, targets)
    result = complete_mixed_cohort(src, plan)
    assert all(t.status == "completed" for t in result.targets)
    assert result.status == "rejected" and result.reason == "all_joint_candidates_failed"
    assert result.selected is None and len(result.attempts) == 1


def test_one_rounding_failure_refuses_entire_mixed_output():
    src, targets = source((("THR", "A", "", ("OG1",), (0., 0., 0.), False),
                           ("ALA", "B", "", ("CB",), (30., 0., 0.), False)))
    water = _format_generated_atom(serial=80, name="O", element="O", position=(-2.169, -1.093, -2.192),
                                   residue_name="HOH", chain_id="W", sequence_number="1", insertion_code="")
    src = src.replace(b"END\n", water.encode()+b"END\n")
    plan = plan_mixed_cohort(src, targets)
    chosen = complete_mixed_cohort(src, plan).selected
    assert chosen is not None
    with pytest.raises(ValueError, match="rounded mixed PDB"):
        materialize_mixed_cohort(src, plan, chosen)


@pytest.mark.parametrize("targets", ((), (ComponentInstanceRef(1, "X", "1"),),
    (ComponentInstanceRef(1, "A", "1"),)*2, tuple(ComponentInstanceRef(1, c, "1") for c in "ABCDE")))
def test_target_set_is_bounded_unique_and_exact(targets):
    src, _ = source(FOUR)
    with pytest.raises(ValueError):
        plan_mixed_cohort(src, targets)


def test_modes_partners_context_and_stereo_records_are_not_caller_authority():
    src, targets = source(FOUR)
    plan = plan_mixed_cohort(src, targets)
    with pytest.raises(ValueError):
        replace(plan.targets[2], mode="paired")
    partner = plan.targets[2].fixed_partner
    changed = replace(partner, xyz=(partner.xyz[0]+.001, *partner.xyz[1:]))
    request = replace(plan.targets[2].request, observations=tuple(changed if o.atom == partner.atom else o
                        for o in plan.targets[2].request.observations))
    altered = replace(plan.targets[2], request=request, fixed_partner=changed)
    forged = replace(plan, targets=(*plan.targets[:2], altered, plan.targets[3]))
    for bad in (forged, replace(plan, targets=tuple(reversed(plan.targets))), replace(plan, source_thr_beta=())):
        with pytest.raises(ValueError, match="binding"):
            complete_mixed_cohort(src, bad)
    with pytest.raises(ValueError):
        replace(plan, policy="canonical-sidechain-cohort-staggered-v1")


def test_unselected_pending_or_wrong_thr_is_never_exempted():
    specs = (*FOUR[:1], ("THR", "X", "", ("OG1",), (30., 0., 0.), False))
    src, targets = source(specs)
    with pytest.raises(ValueError, match="undefined global THR CB"):
        plan_mixed_cohort(src, targets[:1])
    xyz = points("THR")
    cb = xyz["CB"]
    normal = unit(cross(sub(xyz["CA"], cb), sub(xyz["OG1"], cb)))
    xyz["CG2"] = add(xyz["CG2"], scaled(normal, -2*dot(sub(xyz["CG2"], cb), normal)))
    src, targets = source(FOUR[:1])
    extra = "".join("ATOM  "+_format_generated_atom(serial=i+30, name=n, element=n[0], position=(p[0]+30, p[1], p[2]),
        residue_name="THR", chain_id="X", sequence_number="1", insertion_code="")[6:]
        for i, (n, p) in enumerate(xyz.items(), 1)).encode()
    with pytest.raises(ValueError, match="non-R/nearplanar global THR CB"):
        plan_mixed_cohort(src.replace(b"END\n", extra+b"END\n"), targets)


def test_candidate_swaps_record_edits_and_missing_new_atoms_refuse():
    src, targets = source(FOUR)
    plan = plan_mixed_cohort(src, targets)
    chosen = complete_mixed_cohort(src, plan).selected
    out = materialize_mixed_cohort(src, plan, chosen).pdb
    for bad in (chosen[:-1], tuple(reversed(chosen)), list(chosen)):
        with pytest.raises(ValueError):
            materialize_mixed_cohort(src, plan, bad)
    removed = b"".join(l for l in out.splitlines(keepends=True)
                       if not (l.startswith(b"ATOM  ") and l[21:22] == b"B" and l[12:16].strip() == b"CG1"))
    for bad in (removed, out.replace(b"  -1.460", b"  -1.461", 1), out.replace(b"END\n", b"REMARK extra\nEND\n")):
        with pytest.raises(ValueError):
            audit_materialized_mixed_cohort(src, plan, chosen, bad)


def test_budget_limits_and_noop_targets_fail_closed(monkeypatch):
    import dvbfixer.protein_completion.mixed_cohort as module

    src, targets = source(FOUR)
    plan = plan_mixed_cohort(src, targets)
    monkeypatch.setattr(module, "MAX_COMBINATIONS", 1)
    result = complete_mixed_cohort(src, plan)
    assert result.reason == "combination_budget_exceeded" and result.selected is None and not result.attempts
    monkeypatch.setattr(module, "MAX_CONTEXT_ATOMS", 1)
    with pytest.raises(ValueError, match="context exceeds"):
        plan_mixed_cohort(src, targets)
    good, targets = source((("ALA", "A", "", (), (0., 0., 0.), False),))
    with pytest.raises(ValueError, match="require_missing_atoms"):
        plan_mixed_cohort(good, targets)


def test_different_sources_cannot_be_joined():
    src, targets = source(FOUR)
    plan = plan_mixed_cohort(src, targets)
    altered = replace(plan.targets[0], request=replace(plan.targets[0].request, source_sha256="0"*64))
    with pytest.raises(ValueError, match="one source SHA"):
        MixedCohortPlan((altered, *plan.targets[1:]), plan.source_thr_beta)


def test_old_cohort_and_evidence_contracts_cannot_adopt_mixed_plans_or_results(tmp_path):
    # All-standard geometry is valid in both APIs, but policy/plan binding is not
    # interchangeable merely because output coordinates might coincide.
    src, targets = source((FOUR[0], FOUR[3]))
    mixed = plan_mixed_cohort(src, targets)
    result = complete_mixed_cohort(src, mixed)
    with pytest.raises(ValueError, match="cohort plan"):
        complete_cohort(src, mixed)
    with pytest.raises(ValueError, match="original cohort"):
        write_cohort_evidence(src, mixed, result, tmp_path/"bad-kind")
    legacy = plan_cohort(src, targets)
    with pytest.raises(ValueError, match="replay"):
        write_cohort_evidence(src, legacy, result, tmp_path/"wrong-result")
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("specs", (
    (("CYS", "A", "", ("CB", "SG"), (0., 0., 0.), False), FOUR[1], FOUR[2], FOUR[3]),
    tuple(("VAL" if j % 2 else "THR", chr(65+j), "", ("CB", "CG1", "CG2") if j % 2 else
           ("CB", "OG1", "CG2"), (40.*j, 0., 0.), False) for j in range(4)),
))
def test_free_cysteine_and_full_81_combination_branch_bound(specs):
    src, targets = source(specs)
    plan = plan_mixed_cohort(src, targets)
    result = complete_mixed_cohort(src, plan)
    assert result.selected is not None
    if all(t.mode == "paired" for t in plan.targets):
        assert len(result.attempts) == 81 and all(a.accepted for a in result.attempts)
    out = materialize_mixed_cohort(src, plan, result.selected)
    assert len(out.audit.generated) == sum(len(t.request.missing_atoms) for t in plan.targets)
    assert all(c.state == "observed_R" for c in out.audit.thr_beta)


@pytest.mark.parametrize("masks", (((28, ("CB",)), (26, ("CG1", "CG2")), (7, ("OG1",)), (57, ("OG",))),
                                   ((7, ("CB", "OG1", "CG2")), (55, ("CG2",)))))
def test_real_ubiquitin_joint_masks_have_no_pending_output_stereo(masks):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    selected = dict(masks)
    src = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72
                   and l[12:16].strip().decode() not in selected.get(int(l[22:26]), ()))+b"TER\nEND\n"
    targets = tuple(ComponentInstanceRef(1, "A", str(n)) for n in selected)
    plan = plan_mixed_cohort(src, targets)
    result = complete_mixed_cohort(src, plan)
    assert result.selected is not None
    out = materialize_mixed_cohort(src, plan, result.selected)
    assert out.audit == audit_materialized_mixed_cohort(src, plan, result.selected, out.pdb)
    assert all(c.state == "observed_R" for c in out.audit.thr_beta) and len(out.audit.thr_beta) == 7
    assert len(out.audit.generated) == sum(len(mask) for mask in selected.values())
