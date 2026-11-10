"""Private multi-target joint geometry and all-or-nothing source readback."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.cohort import (
    CohortPlan,
    audit_materialized_cohort,
    complete_cohort,
    materialize_cohort,
    plan_cohort,
)
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import materialize_source, plan_source

BACKBONE = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
TARGETS = (ComponentInstanceRef(1, "A", "1"), ComponentInstanceRef(1, "a", "1", "A"))


def source(name="SER", offset=(-3., 0., -4.)):
    lines = ["HEADER    TWO TARGET CONTROL\n"]
    for j, target in enumerate(TARGETS):
        for i, (n, p) in enumerate(BACKBONE.items(), 1):
            xyz = p if j == 0 else (p[0]+offset[0], -p[1]+offset[1], -p[2]+offset[2])
            line = _format_generated_atom(serial=j*10+i, name=n, element=n[0], position=xyz,
                                          residue_name=name, chain_id=target.chain_id,
                                          sequence_number=target.sequence_number, insertion_code=target.insertion_code)
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    return ("".join(lines)+"END\n").encode()


def test_new_atoms_of_other_targets_are_in_joint_context():
    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    assert result.status == "completed"
    assert len(result.attempts) == 6
    failed = next(a for a in result.attempts if not a.accepted)
    assert all(result.targets[j].attempts[i].audit.accepted for j, i in enumerate(failed.indices))
    assert any("generated_clash" in f for a in failed.audits for f in a.findings)
    bad = tuple(result.targets[j].attempts[i].candidate for j, i in enumerate(failed.indices))
    with pytest.raises(ValueError, match="joint generated geometry"):
        materialize_cohort(src, plan, bad)
    output = materialize_cohort(src, plan, result.selected)
    assert len(output.generated) == 4
    assert audit_materialized_cohort(src, plan, result.selected, output.pdb) == output.generated
    record = output.dat_record()
    assert {(a["chain"], a["icode"], a["atom"]) for a in record.added_atoms} == {
        (chain, code, name) for chain, code in (("A", ""), ("a", "A")) for name in ("CB", "OG")}
    assert len({a.atom for a in output.generated}) == 4


def test_target_order_is_canonicalized_before_enumeration_and_serial_allocation():
    src = source()
    a, b = plan_cohort(src, TARGETS), plan_cohort(src, tuple(reversed(TARGETS)))
    assert a == b and a.digest == b.digest
    first, second = complete_cohort(src, a), complete_cohort(src, b)
    assert first == second
    assert materialize_cohort(src, a, first.selected).pdb == materialize_cohort(src, b, second.selected).pdb


@pytest.mark.parametrize("targets", ((), (TARGETS[0], TARGETS[0]), list(TARGETS),
    tuple(ComponentInstanceRef(1, chr(65+i), "1") for i in range(5)), (ComponentInstanceRef(1, "X", "1"),)))
def test_target_admission_is_bounded_unique_and_exact(targets):
    with pytest.raises(ValueError):
        plan_cohort(source(), targets)


@pytest.mark.parametrize("tamper", (
    lambda p: CohortPlan(tuple(reversed(p.requests))),
    lambda p: CohortPlan((replace(p.requests[0], observations=p.requests[0].observations[:-1]), *p.requests[1:])),
))
def test_plan_source_context_and_order_are_rederived(tamper):
    src = source()
    plan = plan_cohort(src, TARGETS)
    with pytest.raises(ValueError, match="binding"):
        complete_cohort(src, tamper(plan))


def test_single_target_gate_does_not_inherit_cohort_exemptions():
    with pytest.raises(ValueError, match="undefined global"):
        plan_source(source(), TARGETS[0])


def test_contract_cannot_combine_different_sources():
    plan = plan_cohort(source(), TARGETS)
    with pytest.raises(ValueError, match="one source SHA"):
        CohortPlan((replace(plan.requests[0], source_sha256="0"*64), plan.requests[1]))


def test_two_individual_ala_passes_can_have_no_joint_solution():
    src = source("ALA", (0., -1.538, -5.))
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    assert all(r.status == "completed" for r in result.targets)
    assert result.status == "rejected" and result.reason == "all_joint_candidates_failed"
    assert result.selected is None and len(result.attempts) == 1


def test_context_and_combination_budgets_fail_closed(monkeypatch):
    import dvbfixer.protein_completion.cohort as module

    src = source()
    plan = plan_cohort(src, TARGETS)
    monkeypatch.setattr(module, "MAX_COMBINATIONS", 1)
    result = complete_cohort(src, plan)
    assert result.status == "rejected" and result.reason == "combination_budget_exceeded"
    assert result.selected is None and not result.attempts
    monkeypatch.setattr(module, "MAX_CONTEXT_ATOMS", 8)
    with pytest.raises(ValueError, match="context exceeds"):
        plan_cohort(src, TARGETS)


def test_four_serine_controls_enumerate_all_81_combinations():
    lines = []
    targets = tuple(ComponentInstanceRef(1, chain, "1") for chain in "ABCD")
    for j, target in enumerate(targets):
        for i, (name, xyz) in enumerate(BACKBONE.items(), 1):
            text = _format_generated_atom(serial=j*10+i, name=name, element=name[0],
                                          position=(xyz[0]+40*j, xyz[1], xyz[2]), residue_name="SER",
                                          chain_id=target.chain_id, sequence_number="1", insertion_code="")
            lines.append("ATOM  "+text[6:])
        lines.append("TER\n")
    src = ("".join(lines)+"END\n").encode()
    plan = plan_cohort(src, targets)
    result = complete_cohort(src, plan)
    assert len(result.attempts) == 81 and all(a.accepted for a in result.attempts)
    assert len(materialize_cohort(src, plan, result.selected).generated) == 8


def test_individual_target_failure_retains_attempts_without_joint_output():
    src = source().replace(b"   1.760", b"  10.760", 1)
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    assert result.status == "rejected" and result.reason == "individual_target_failed"
    assert result.selected is None and result.attempts == ()
    assert result.targets[0].attempts and all(not a.audit.accepted for a in result.targets[0].attempts)


def test_one_target_cohort_matches_existing_single_target_writer():
    src = source("ALA", (30., 0., 0.))
    # Supply CB to the other residue explicitly, then complete only the first.
    second = plan_cohort(src, TARGETS)
    ready = materialize_cohort(src, second, complete_cohort(src, second).selected).pdb
    src = b"".join(l for l in ready.splitlines(keepends=True)
                   if not (l.startswith(b"ATOM  ") and l[21:22] == b"A" and l[12:16].strip() == b"CB"))
    single = plan_source(src, TARGETS[0])
    single_candidate = complete_sidechain(single).selected
    plan = plan_cohort(src, (TARGETS[0],))
    cohort = complete_cohort(src, plan)
    assert materialize_cohort(src, plan, cohort.selected).pdb == materialize_source(src, single, single_candidate).pdb
    with pytest.raises(ValueError, match="each have missing"):
        plan_cohort(ready, TARGETS)


def test_candidate_swaps_missing_targets_and_mutation_refuse():
    src = source()
    plan = plan_cohort(src, TARGETS)
    chosen = complete_cohort(src, plan).selected
    for candidates in (tuple(reversed(chosen)), chosen[:1], list(chosen),
                       (replace(chosen[0], request_digest="0"*64), chosen[1])):
        with pytest.raises(ValueError):
            materialize_cohort(src, plan, candidates)


def test_no_partial_accepted_structure_when_rounding_breaks_one_target():
    src = source("ALA", (30., 0., 0.))
    water = _format_generated_atom(serial=25, name="O", element="O", position=(.589, -.6, 1.405),
                                   residue_name="HOH", chain_id="W", sequence_number="1", insertion_code="")
    src = src.replace(b"END\n", water.encode()+b"END\n")
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    assert result.selected is not None
    with pytest.raises(ValueError, match="rounded joint PDB"):
        materialize_cohort(src, plan, result.selected)


def test_independent_output_audit_rejects_observation_edits_and_missing_additions():
    src = source()
    plan = plan_cohort(src, TARGETS)
    chosen = complete_cohort(src, plan).selected
    out = materialize_cohort(src, plan, chosen).pdb
    removed = b"".join(l for l in out.splitlines(keepends=True)
                       if not (l.startswith(b"ATOM  ") and l[21:22] == b"A" and l[12:16].strip() == b"OG"))
    for bad in (removed, out.replace(b"  -1.460", b"  -1.461"),
                out.replace(b"END\n", b"REMARK added\nEND\n")):
        with pytest.raises(ValueError):
            audit_materialized_cohort(src, plan, chosen, bad)


@pytest.mark.parametrize("residues", ((28, 46), (20, 57, 65), (20, 28, 46, 57)))
def test_real_ubiquitin_projection_multiple_missing_sidechains(residues):
    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    src = b"".join(l for l in original.splitlines(keepends=True)
                   if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72
                   and not (int(l[22:26]) in residues and l[12:16].strip() in {b"CB", b"OG"}))+b"TER\nEND\n"
    targets = tuple(ComponentInstanceRef(1, "A", str(n)) for n in residues)
    plan = plan_cohort(src, targets)
    result = complete_cohort(src, plan)
    assert result.selected is not None
    out = materialize_cohort(src, plan, result.selected)
    assert len(out.generated) == sum(len(r.missing_atoms) for r in plan.requests)
    assert audit_materialized_cohort(src, plan, result.selected, out.pdb) == out.generated
