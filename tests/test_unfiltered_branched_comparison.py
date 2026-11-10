"""Planning does not inspect own success; scoring never exchanges fixed partners."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from scripts import unfiltered_branched_comparison as module

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import (
    Addition,
    CompletionRequest,
    CompletionResult,
    Observation,
)


def test_plan_freezes_all_48_masks_without_running_or_filtering_generators(tmp_path, monkeypatch):
    def never(*args, **kwargs):
        raise AssertionError("planning must not run admission or generation")

    for name in ("plan_anchored_source", "plan_branched_source", "complete_anchored_branch", "complete_sidechain", "template_proposal", "measure"):
        monkeypatch.setattr(module, name, never)
    module.plan(Path(__file__).resolve().parents[1], tmp_path/"plan")
    record = json.loads((tmp_path/"plan/plan.json").read_text())
    assert record["case_count"] == 48 and len(record["cases"]) == 48 and record["no_own_success_filter"]
    assert {c["fixture"] for c in record["cases"]} == {"1EMV", "1TIM", "1DQJ"}
    assert sum(c["mode"] == "anchored" for c in record["cases"]) == 24
    for fixture in module.FIXTURES:
        cases = [c for c in record["cases"] if c["fixture"] == fixture]
        assert len(cases) == 16
        assert len({tuple(c["target"].items()) for c in cases}) == 4
    for name, digest in record["artifacts"].items():
        assert hashlib.sha256((tmp_path/"plan"/name).read_bytes()).hexdigest() == digest
    # The altloc/occupancy source remains in the plan rather than being excluded.
    lines = (tmp_path/"plan/1DQJ-reference.pdb").read_text().splitlines()
    assert any(l.startswith("ATOM  ") and l[16].strip() for l in lines)


def test_frozen_input_tampering_refuses_before_execution_or_output(tmp_path, monkeypatch):
    module.plan(Path(__file__).resolve().parents[1], tmp_path/"plan")
    record = json.loads((tmp_path/"plan/plan.json").read_text())
    path = tmp_path/"plan"/record["cases"][0]["source"]
    path.write_bytes(path.read_bytes()+b"REMARK tampered\n")

    def never():
        raise AssertionError("input verification must precede preflight")

    monkeypatch.setattr(module, "preflight", never)
    with pytest.raises(ValueError, match="input changed"):
        module.run(tmp_path/"plan", tmp_path/"out", tmp_path/"summary.json")
    assert not (tmp_path/"out").exists()


def test_failed_own_arm_does_not_prevent_template_execution(tmp_path, monkeypatch):
    module.plan(Path(__file__).resolve().parents[1], tmp_path/"plan")
    calls = []

    def rejected(request):
        return CompletionResult("rejected", request.digest, "injected-own-rejection")

    def template(source, request, *, retain_snapshot):
        calls.append(request.digest)
        raise ValueError("injected mock template protocol error, not a science failure")

    monkeypatch.setattr(module, "preflight", lambda: {"mock": True})
    monkeypatch.setattr(module, "complete_sidechain", rejected)
    monkeypatch.setattr(module, "complete_anchored_branch", rejected)
    monkeypatch.setattr(module, "template_proposal", template)
    module.run(tmp_path/"plan", tmp_path/"out", tmp_path/"summary.json")
    record = json.loads((tmp_path/"summary.json").read_text())
    admitted = [r for r in record["results"].values() if r["status"] == "admitted"]
    assert admitted and len(calls) == len(admitted)
    assert all(r["own"]["status"] == "rejected" for r in admitted)
    assert all(r["template"]["status"] == "execution_or_protocol_blocker" for r in admitted)


def request_and_truth(mask):
    target = ComponentInstanceRef(1, "a", "82", "A")
    truth = {"N": (-1., 0., 0.), "CA": (0., 0., 0.), "C": (0., 1., 0.), "O": (0., 2., 0.),
             "CB": (0., 0., 1.), "CG1": (1., 0., 2.), "CG2": (-1., 0., 2.)}
    observed = tuple(Observation(ExactAtomRef(target, n), n[0], p) for n, p in truth.items() if n not in mask)
    return CompletionRequest("a"*64, target, "VAL", observed, mask), truth


def test_val_symmetry_is_reporting_only_and_requires_both_methyls_generated():
    request, truth = request_and_truth(("CG1", "CG2"))
    additions = tuple(Addition(ExactAtomRef(request.residue, n), "C", truth["CG2" if n == "CG1" else "CG1"])
                      for n in request.missing_atoms)
    result = module.reporting_scores(request, additions, truth)
    assert result["exact_name_rmsd_angstrom"] == 2.
    assert result["generated_only_val_symmetry_rmsd_angstrom"] == 0.
    assert result["val_swap_admissible"] and result["val_swapped_in_reporting"]
    assert result["reporting_only_no_reranking"] and not result["used_for_selection"]
    single, truth = request_and_truth(("CG1",))
    result = module.reporting_scores(single, additions[:1], truth)
    assert result["generated_only_val_symmetry_rmsd_angstrom"] == 2.
    assert not result["val_swap_admissible"] and not result["val_swapped_in_reporting"]


def test_native_reference_gate_failure_is_visible_and_not_called_wrong_stereo():
    from scripts.compare_branched_completion import measure

    from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines

    original = (Path(__file__).resolve().parents[1]/"tests/fixtures/protein_context/1TIM.pdb").read_bytes()
    reference = b"".join(l for l in original.splitlines(keepends=True) if l.startswith(b"ATOM  ") and l[21:22] == b"A")
    source = b"".join(l for l in reference.splitlines(keepends=True)
                      if not (l[22:26].strip() == b"39" and l[12:16].strip() == b"CG1"))+b"TER\nEND\n"
    target = ComponentInstanceRef(1, "A", "39")
    plan = module.plan_anchored_source(source, target)
    truth = {a.name: a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True))
             if _component_ref(a) == target}
    native = (Addition(ExactAtomRef(target, "CG1"), "C", truth["CG1"]),)
    audit, _ = measure(source, plan.request, native, anchored=True)
    assert not audit.accepted and audit.rdkit_target_cip == {"CA": "S"}
    assert audit.global_stereo_findings == ()
    assert audit.raw["findings"] and all(f.startswith("generated_clash:CG1:") for f in audit.raw["findings"])
    assert module.reporting_scores(plan.request, native, truth)["exact_name_rmsd_angstrom"] == 0
    own = module.complete_anchored_branch(plan.request)
    assert own.status == "rejected" and own.selected is None
    assert module.reporting_scores(plan.request, own.attempts[0].candidate.additions, truth)["exact_name_rmsd_angstrom"] < .02
