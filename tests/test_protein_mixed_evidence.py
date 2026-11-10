"""Mixed-only evidence: complete replay, fixed partners, failures and safe commit."""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.branched import TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.cohort import complete_cohort, plan_cohort
from dvbfixer.protein_completion.evidence import read_cohort_evidence, write_cohort_evidence
from dvbfixer.protein_completion.geometry import extend, l_cb
from dvbfixer.protein_completion.mixed_cohort import complete_mixed_cohort, plan_mixed_cohort
from dvbfixer.protein_completion.mixed_evidence import read_mixed_evidence, write_mixed_evidence

TARGETS = (ComponentInstanceRef(1, "A", "1"), ComponentInstanceRef(1, "a", "1", "A"))


def source(mode="accepted"):
    lines = []
    for j, target in enumerate(TARGETS):
        xyz = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
        cb = xyz["CB"] = l_cb(xyz["N"], xyz["CA"], xyz["C"])
        xyz["OG1"] = extend(xyz["N"], xyz["CA"], cb, 1.43, TETRAHEDRAL_ANGLE, 60.)
        xyz["CG2"] = extend(xyz["N"], xyz["CA"], cb, 1.53, TETRAHEDRAL_ANGLE, -60.)
        missing = {"CB", "OG1", "CG2"} if j == 0 else {"CG2"}
        if mode == "joint":
            missing = {"CG2"}
        if mode == "rounded" and j == 0:
            missing = {"OG1"}
        for i, (name, p) in enumerate(xyz.items(), 1):
            if name in missing:
                continue
            point = p
            if j:
                point = (p[0]+.8, -p[1]-4.5, -p[2]-4.95) if mode == "joint" else (p[0]+30, p[1], p[2])
            line = _format_generated_atom(serial=j*20+i, name=name, element=name[0], position=point,
                                          residue_name="THR", chain_id=target.chain_id, sequence_number="1",
                                          insertion_code=target.insertion_code)
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    if mode == "rounded":
        lines.append(_format_generated_atom(serial=80, name="O", element="O", position=(-2.169, -1.093, -2.192),
                                             residue_name="HOH", chain_id="W", sequence_number="1", insertion_code=""))
    data = ("".join(lines)+"END\n").encode()
    return data.replace(b"  1.00", b"  0.50", 1) if mode == "admission" else data


def publish(tmp_path, mode="accepted"):
    src = source(mode)
    path = write_mixed_evidence(src, TARGETS, tmp_path/"bundle")
    return src, path


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode()


def rehash(path, name, data):
    (path/name).write_bytes(data)
    index = json.loads((path/"bundle.json").read_bytes())
    index["artifacts"][name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    (path/"bundle.json").write_bytes(encoded(index))


def test_mixed_replay_preserves_modes_partner_and_exact_output(tmp_path):
    src, path = publish(tmp_path)
    evidence = read_mixed_evidence(path, expected_source_sha256=hashlib.sha256(src).hexdigest())
    assert evidence.status == "accepted" and evidence.stage == "rounded" and evidence.refusal is None
    assert [t.mode for t in evidence.plan.targets] == ["paired", "anchored"]
    partner = evidence.plan.targets[1].fixed_partner
    record = evidence.materialized.dat_record()
    assert ("a", "1", "A", partner.atom.name) not in record.added_keys()
    assert record.added_keys() == {("A", "1", "", n) for n in ("CB", "OG1", "CG2")} | {("a", "1", "A", "CG2")}
    assert (path/"source.pdb").read_bytes() == src
    assert evidence.result == complete_mixed_cohort(src, plan_mixed_cohort(src, TARGETS))


@pytest.mark.parametrize("mode,stage", (("admission", "admission"), ("joint", "raw"), ("rounded", "rounded")))
def test_all_failure_stages_retain_evidence_without_accepted_files(tmp_path, mode, stage):
    src, path = publish(tmp_path, mode)
    evidence = read_mixed_evidence(path)
    assert evidence.status == "rejected" and evidence.stage == stage and evidence.refusal
    assert evidence.materialized is None
    assert (path/"source.pdb").read_bytes() == src
    assert not (path/"candidate.pdb").exists() and not (path/"candidate.dat").exists()
    if stage == "admission":
        assert evidence.plan is None and evidence.result is None and not (path/"plan.json").exists()
    else:
        assert evidence.plan is not None and evidence.result is not None
        assert (path/"attempts.json").is_file()


def test_each_reader_rejects_the_other_kind_and_public_readers_reject_mixed(tmp_path):
    from dvbfixer.model.diffusion.atom_completion_publication import read_completion_bundle
    from dvbfixer.model.diffusion.publication import PublicationError, read_published_bundle

    _, path = publish(tmp_path)
    for reader in (read_cohort_evidence, read_completion_bundle):
        with pytest.raises(ValueError):
            reader(path)
    with pytest.raises(PublicationError):
        read_published_bundle(path)
    # A one-SER input fits the original policy; do not relabel it mixed implicitly.
    src = source().split(b"TER\n")[0].replace(b"THR", b"SER")+b"TER\nEND\n"
    plan = plan_cohort(src, TARGETS[:1])
    result = complete_cohort(src, plan)
    old = write_cohort_evidence(src, plan, result, tmp_path/"old")
    with pytest.raises(ValueError, match="index"):
        read_mixed_evidence(old)
    with pytest.raises(ValueError, match="mixed plan/result"):
        write_mixed_evidence(src, TARGETS[:1], tmp_path/"wrong-plan", plan=plan, result=result)


@pytest.mark.parametrize("tamper", (
    lambda r: replace(r, attempts=r.attempts[:-1]),
    lambda r: replace(r, selected_index=999),
    lambda r: replace(r, reason="native-qualified"),
    lambda r: replace(r, targets=(replace(r.targets[0], attempts=()), r.targets[1])),
    lambda r: replace(r, attempts=(replace(r.attempts[0], audits=()), *r.attempts[1:])),
))
def test_retained_attempt_omissions_and_selection_forgery_refuse_before_staging(tmp_path, tamper):
    src = source()
    plan = plan_mixed_cohort(src, TARGETS)
    result = complete_mixed_cohort(src, plan)
    with pytest.raises(ValueError, match="replay"):
        write_mixed_evidence(src, TARGETS, tmp_path/"forged", plan=plan, result=tamper(result))
    assert not list(tmp_path.iterdir())


def test_partner_and_observation_coherent_forgery_is_not_source_authority(tmp_path):
    src = source()
    plan = plan_mixed_cohort(src, TARGETS)
    target = plan.targets[1]
    partner = target.fixed_partner
    forged_partner = replace(partner, xyz=(partner.xyz[0]+.001, *partner.xyz[1:]))
    request = replace(target.request, observations=tuple(forged_partner if o.atom == partner.atom else o for o in target.request.observations))
    forged = replace(plan, targets=(plan.targets[0], replace(target, request=request, fixed_partner=forged_partner)))
    with pytest.raises(ValueError, match="replay"):
        write_mixed_evidence(src, TARGETS, tmp_path/"forged", plan=forged, result=complete_mixed_cohort(src, plan))
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("name", ("source.pdb", "targets.json", "plan.json", "attempts.json", "audit.json", "candidate.pdb", "candidate.dat"))
def test_every_artifact_is_digest_bound(tmp_path, name):
    _, path = publish(tmp_path)
    (path/name).write_bytes((path/name).read_bytes()+b" ")
    with pytest.raises(ValueError):
        read_mixed_evidence(path)


@pytest.mark.parametrize("name,change", (
    ("plan.json", lambda d: d["targets"][1].update(mode="paired")),
    ("plan.json", lambda d: d["targets"][1]["fixed_partner"].update(xyz=[0., 0., 0.])),
    ("plan.json", lambda d: d.update(source_thr_beta=[])),
    ("attempts.json", lambda d: d.update(attempts=[])),
    ("attempts.json", lambda d: d.update(selected_index=999)),
    ("audit.json", lambda d: d["rounded_audit"].update(thr_beta=[])),
    ("audit.json", lambda d: d.update(scientific_qualification=True)),
    ("candidate.dat", lambda d: d.update(total_added=999)),
    ("candidate.dat", lambda d: d["added_atoms"].append({"chain": "a", "resid": "1", "icode": "A", "resname": "THR", "atom": "OG1", "element": "O"})),
    ("targets.json", lambda d: d[1].update(insertion_code="B")),
))
def test_rehashed_semantic_tampering_still_fails_independent_replay(tmp_path, name, change):
    _, path = publish(tmp_path)
    value = json.loads((path/name).read_bytes())
    change(value)
    rehash(path, name, encoded(value))
    with pytest.raises(ValueError):
        read_mixed_evidence(path)


@pytest.mark.parametrize("field,value", (("kind", "experimental-protein-completion-cohort-evidence"),
    ("policy", "canonical-sidechain-cohort-staggered-v1"), ("schema_version", True),
    ("public_profile_enabled", True), ("scientific_qualification", True), ("stage", "admission")))
def test_index_cannot_promote_or_relabel_mixed_evidence(tmp_path, field, value):
    _, path = publish(tmp_path)
    index = json.loads((path/"bundle.json").read_bytes())
    index[field] = value
    (path/"bundle.json").write_bytes(encoded(index))
    with pytest.raises(ValueError):
        read_mixed_evidence(path)


def test_rehashed_changed_observation_in_pdb_refuses(tmp_path):
    _, path = publish(tmp_path)
    rehash(path, "candidate.pdb", (path/"candidate.pdb").read_bytes().replace(b"  -1.460", b"  -1.461", 1))
    with pytest.raises(ValueError, match="replay"):
        read_mixed_evidence(path)


def test_admission_reason_must_match_recomputed_refusal(tmp_path):
    _, path = publish(tmp_path, "admission")
    audit = json.loads((path/"audit.json").read_bytes())
    audit["refusal"] = "some unrelated problem"
    rehash(path, "audit.json", encoded(audit))
    with pytest.raises(ValueError, match="replay"):
        read_mixed_evidence(path)


def test_exact_inventory_readiness_source_pin_and_no_overwrite(tmp_path):
    src, path = publish(tmp_path)
    before = (path/"bundle.json").read_bytes()
    with pytest.raises(FileExistsError):
        write_mixed_evidence(src, TARGETS, path)
    assert (path/"bundle.json").read_bytes() == before
    with pytest.raises(ValueError, match="externally expected"):
        read_mixed_evidence(path, expected_source_sha256="0"*64)
    (path/"extra.log").write_text("keep")
    with pytest.raises(ValueError, match="flat declared"):
        read_mixed_evidence(path)
    (path/"bundle.json").unlink()
    with pytest.raises(OSError):
        read_mixed_evidence(path)
    with pytest.raises(FileExistsError):
        write_mixed_evidence(src, TARGETS, path)
    assert (path/"extra.log").read_text() == "keep"


def test_symlinks_hardlinks_and_traversal_refs_are_refused(tmp_path):
    _, path = publish(tmp_path)
    original = path/"source.pdb"
    outside = tmp_path/"outside.pdb"
    outside.write_bytes(original.read_bytes())
    original.unlink()
    original.symlink_to(outside)
    with pytest.raises(ValueError, match="symlink"):
        read_mixed_evidence(path)
    original.unlink()
    original.hardlink_to(outside)
    with pytest.raises(ValueError, match="single-linked"):
        read_mixed_evidence(path)
    index = json.loads((path/"bundle.json").read_bytes())
    index["artifacts"]["../outside.pdb"] = index["artifacts"].pop("source.pdb")
    (path/"bundle.json").write_bytes(encoded(index))
    with pytest.raises(ValueError, match="flat declared"):
        read_mixed_evidence(path)


def test_staging_failure_retains_diagnostics_without_readiness(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.mixed_evidence as module

    write = module._durable_write

    def fail(path, data):
        if path.name == "candidate.dat":
            raise OSError("injected I/O failure")
        write(path, data)

    monkeypatch.setattr(module, "_durable_write", fail)
    with pytest.raises(OSError, match="injected"):
        write_mixed_evidence(source(), TARGETS, tmp_path/"failed")
    assert not (tmp_path/"failed").exists()
    stages = list(tmp_path.glob(".mixed-cohort-staging-*"))
    assert len(stages) == 1 and (stages[0]/"attempts.json").is_file()
    assert not (stages[0]/"bundle.json").exists()


def test_unknown_materialization_failure_is_not_scientific_rejection(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.mixed_evidence as module

    def fail(*args, **kwargs):
        raise module.SourceCompletionError("unexpected materializer bug")

    monkeypatch.setattr(module, "materialize_mixed_cohort", fail)
    with pytest.raises(ValueError, match="unexpected"):
        write_mixed_evidence(source(), TARGETS, tmp_path/"failed")
    assert not list(tmp_path.iterdir())


def test_bool_int_forgery_and_artifact_limits_fail_before_staging(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.mixed_evidence as module

    src = source()
    plan = plan_mixed_cohort(src, TARGETS)
    result = complete_mixed_cohort(src, plan)
    attempt = result.attempts[0]
    forged = replace(result, attempts=(replace(attempt, audits=(replace(attempt.audits[0], accepted=1), *attempt.audits[1:])), *result.attempts[1:]))
    assert forged == result
    with pytest.raises(ValueError, match="replay"):
        write_mixed_evidence(src, TARGETS, tmp_path/"bad", plan=plan, result=forged)
    monkeypatch.setattr(module, "MAX_ARTIFACT_BYTES", 1)
    with pytest.raises(ValueError, match="artifact budget"):
        write_mixed_evidence(src, TARGETS, tmp_path/"oversized")
    assert not list(tmp_path.iterdir())
    monkeypatch.setattr(module, "MAX_ARTIFACT_BYTES", 64*1024*1024)
    monkeypatch.setattr(module, "MAX_TOTAL_BYTES", 1)
    with pytest.raises(ValueError, match="artifact budget"):
        write_mixed_evidence(src, TARGETS, tmp_path/"total-limit")
    assert not list(tmp_path.iterdir())


def test_competing_and_crash_left_reservations_are_not_adopted(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.mixed_evidence as module

    src = source()
    crash = tmp_path/"crash-left"
    crash.mkdir()
    (crash/"user.log").write_text("preserve")
    with pytest.raises(FileExistsError):
        write_mixed_evidence(src, TARGETS, crash)
    commit = module.commit_manifest_last

    def competing(staging, destination):
        destination.mkdir()
        (destination/"user.log").write_text("preserve")
        commit(staging, destination)

    monkeypatch.setattr(module, "commit_manifest_last", competing)
    with pytest.raises(FileExistsError):
        write_mixed_evidence(src, TARGETS, tmp_path/"race")
    assert (crash/"user.log").read_text() == (tmp_path/"race/user.log").read_text() == "preserve"
    assert not (tmp_path/"race/bundle.json").exists()
    staging = list(tmp_path.glob(".mixed-cohort-staging-*"))
    assert len(staging) == 1 and (staging[0]/"attempts.json").is_file()


@pytest.mark.parametrize("targets", ((TARGETS[0], TARGETS[0]), (ComponentInstanceRef(1, "X", "1"),)))
def test_bounded_target_admission_refusals_are_replayable(tmp_path, targets):
    path = write_mixed_evidence(source(), targets, tmp_path/"rejected")
    audit = read_mixed_evidence(path)
    assert audit.stage == "admission" and audit.status == "rejected"
    assert audit.plan is None and audit.result is None


def test_retained_plan_result_pair_and_order_are_not_optional_authority(tmp_path):
    src = source()
    plan = plan_mixed_cohort(src, TARGETS)
    result = complete_mixed_cohort(src, plan)
    for kwargs in ({"plan": plan}, {"result": result},
                   {"plan": replace(plan, targets=tuple(reversed(plan.targets))), "result": result}):
        with pytest.raises(ValueError):
            write_mixed_evidence(src, TARGETS, tmp_path/"bad", **kwargs)
    assert not list(tmp_path.iterdir())


def test_reader_rejects_a_symlink_root(tmp_path):
    _, path = publish(tmp_path)
    alias = tmp_path/"alias"
    alias.symlink_to(path, target_is_directory=True)
    with pytest.raises(ValueError, match="regular directory"):
        read_mixed_evidence(alias)


def test_individual_failure_retains_failed_proposals_without_joint_output(tmp_path):
    src = source().replace(b"   1.760", b"  10.760", 1)
    path = write_mixed_evidence(src, TARGETS, tmp_path/"individual-failure")
    evidence = read_mixed_evidence(path)
    assert evidence.status == "rejected" and evidence.stage == "raw"
    assert evidence.result.reason == "individual_target_failed"
    assert evidence.result.targets[0].attempts and evidence.result.attempts == ()
    assert all(not a.audit.accepted for a in evidence.result.targets[0].attempts)
    assert evidence.materialized is None and not (path/"candidate.pdb").exists()


@pytest.mark.parametrize("field,value", (("bytes", True), ("bytes", 64*1024*1024+1),
                                        ("sha256", "../not-a-digest")))
def test_invalid_artifact_references_are_rejected_before_payload_read(tmp_path, field, value):
    _, path = publish(tmp_path)
    index = json.loads((path/"bundle.json").read_bytes())
    index["artifacts"]["source.pdb"][field] = value
    (path/"bundle.json").write_bytes(encoded(index))
    with pytest.raises(ValueError, match="artifact reference"):
        read_mixed_evidence(path)


def test_unexpected_planning_io_failure_is_not_an_admission_refusal(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.mixed_evidence as module

    def fail(*args, **kwargs):
        raise OSError("unexpected source-planning I/O error")

    monkeypatch.setattr(module, "plan_mixed_cohort", fail)
    with pytest.raises(OSError, match="unexpected"):
        write_mixed_evidence(source(), TARGETS, tmp_path/"failed")
    assert not list(tmp_path.iterdir())


def test_rehashed_target_boolean_fields_refuse_strict_identity_transport(tmp_path):
    _, path = publish(tmp_path)
    targets = json.loads((path/"targets.json").read_bytes())
    targets[0]["model"] = True
    rehash(path, "targets.json", encoded(targets))
    with pytest.raises(ValueError, match="strict exact target"):
        read_mixed_evidence(path)
