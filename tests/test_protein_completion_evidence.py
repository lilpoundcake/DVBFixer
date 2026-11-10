"""Private cohort evidence: fail-closed replay, ownership and no-overwrite commit."""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace

import pytest

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.cohort import complete_cohort, plan_cohort
from dvbfixer.protein_completion.evidence import read_cohort_evidence, write_cohort_evidence

TARGETS = (ComponentInstanceRef(1, "A", "1"), ComponentInstanceRef(1, "a", "1", "A"))


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)+"\n").encode()


def source(name="SER", offset=(-3., 0., -4.), water=False):
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    lines = []
    for j, target in enumerate(TARGETS):
        for i, (n, p) in enumerate(points.items(), 1):
            xyz = p if j == 0 else (p[0]+offset[0], -p[1]+offset[1], -p[2]+offset[2])
            line = _format_generated_atom(serial=j*10+i, name=n, element=n[0], position=xyz,
                                          residue_name=name, chain_id=target.chain_id, sequence_number="1",
                                          insertion_code=target.insertion_code)
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    if water:
        lines.append(_format_generated_atom(serial=25, name="O", element="O", position=(.589, -.6, 1.405),
                                             residue_name="HOH", chain_id="W", sequence_number="1", insertion_code=""))
    return ("".join(lines)+"END\n").encode()


def publish(tmp_path, src=None):
    src = source() if src is None else src
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    path = write_cohort_evidence(src, plan, result, tmp_path/"evidence")
    return src, plan, result, path


def rehash(path, name, data):
    (path/name).write_bytes(data)
    index = json.loads((path/"bundle.json").read_bytes())
    index["artifacts"][name] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    (path/"bundle.json").write_bytes(encoded(index))


def test_accepted_bundle_replays_attempts_output_and_generated_sidecar(tmp_path):
    src, plan, result, path = publish(tmp_path)
    audit = read_cohort_evidence(path, expected_source_sha256=hashlib.sha256(src).hexdigest())
    assert audit.status == "accepted" and audit.stage == "rounded"
    assert audit.plan == plan and audit.result == result
    assert audit.materialized is not None and len(audit.materialized.generated) == 4
    assert (path/"source.pdb").read_bytes() == src
    assert any(not a.accepted for a in audit.result.attempts)
    assert not list(tmp_path.glob(".protein-cohort-staging-*"))


def test_public_reader_does_not_accept_experimental_cohort_kind(tmp_path):
    from dvbfixer.model.diffusion.atom_completion_publication import read_completion_bundle
    from dvbfixer.model.diffusion.publication import PublicationError, read_published_bundle

    _, _, _, path = publish(tmp_path)
    with pytest.raises(PublicationError):
        read_published_bundle(path)
    with pytest.raises(ValueError):
        read_completion_bundle(path)


@pytest.mark.parametrize("src,stage", ((source("ALA", (0., -1.538, -5.)), "raw"),
                                       (source("ALA", (30., 0., 0.), water=True), "rounded")))
def test_rejections_retain_source_and_all_attempts_without_accepted_structure(tmp_path, src, stage):
    _, plan, result, path = publish(tmp_path, src)
    audit = read_cohort_evidence(path)
    assert audit.status == "rejected" and audit.stage == stage
    assert audit.plan == plan and audit.result == result and audit.materialized is None
    assert not (path/"candidate.pdb").exists() and not (path/"candidate.dat").exists()
    assert (path/"attempts.json").is_file() and (path/"source.pdb").read_bytes() == src


@pytest.mark.parametrize("change", (
    lambda r: replace(r, selected_index=1),
    lambda r: replace(r, attempts=r.attempts[:-1]),
    lambda r: replace(r, attempts=(replace(r.attempts[0], audits=()), *r.attempts[1:])),
    lambda r: replace(r, targets=(replace(r.targets[0], attempts=()), r.targets[1])),
    lambda r: replace(r, reason="qualified-native-quality"),
))
def test_forged_or_incomplete_attempts_refuse_before_staging(tmp_path, change):
    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    with pytest.raises(ValueError, match="replay"):
        write_cohort_evidence(src, plan, change(result), tmp_path/"bad")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name", ("source.pdb", "plan.json", "targets.json", "attempts.json", "audit.json", "candidate.pdb", "candidate.dat"))
def test_every_payload_is_digest_bound(tmp_path, name):
    _, _, _, path = publish(tmp_path)
    with (path/name).open("ab") as handle:
        handle.write(b" ")
    with pytest.raises(ValueError):
        read_cohort_evidence(path)


@pytest.mark.parametrize("name,change", (
    ("attempts.json", lambda d: d.update(selected_index=1)),
    ("attempts.json", lambda d: d.update(attempts=[])),
    ("audit.json", lambda d: d.update(scientific_qualification=True)),
    ("candidate.dat", lambda d: d.update(total_added=999)),
    ("candidate.dat", lambda d: d["added_atoms"].append({"chain": "A", "resid": "1", "icode": "", "resname": "SER", "atom": "N", "element": "N"})),
    ("targets.json", lambda d: d[0].update(chain_id="a")),
    ("plan.json", lambda d: d["requests"][0].update(links=[])),
))
def test_rehashed_semantic_tampering_is_rejected_by_replay(tmp_path, name, change):
    _, _, _, path = publish(tmp_path)
    value = json.loads((path/name).read_bytes())
    before = encoded(value)
    change(value)
    after = encoded(value)
    if after == before:
        # Empty links in this synthetic source are already valid; mutate the
        # missing mask rather than pretending an unchanged payload was tampered.
        value["requests"][0]["missing_atoms"] = []
        after = encoded(value)
    assert after != before
    rehash(path, name, after)
    with pytest.raises(ValueError):
        read_cohort_evidence(path)


def test_rehashed_observation_coordinate_mutation_is_not_accepted(tmp_path):
    _, _, _, path = publish(tmp_path)
    data = (path/"candidate.pdb").read_bytes().replace(b"  -1.460", b"  -1.461", 1)
    rehash(path, "candidate.pdb", data)
    with pytest.raises(ValueError, match="replay"):
        read_cohort_evidence(path)


@pytest.mark.parametrize("field,value", (("schema_version", 2), ("schema_version", True),
    ("kind", "diffusion-bundle"), ("public_profile_enabled", True),
    ("scientific_qualification", True), ("status", "accepted-native"), ("policy", "unversioned")))
def test_index_cannot_promote_or_relabel_private_controls(tmp_path, field, value):
    _, _, _, path = publish(tmp_path)
    index = json.loads((path/"bundle.json").read_bytes())
    index[field] = value
    (path/"bundle.json").write_bytes(encoded(index))
    with pytest.raises(ValueError):
        read_cohort_evidence(path)


def test_missing_index_or_extra_files_are_not_ready_evidence(tmp_path):
    _, _, _, path = publish(tmp_path)
    (path/"unexpected.log").write_text("not declared")
    with pytest.raises(ValueError, match="exactly"):
        read_cohort_evidence(path)
    (path/"bundle.json").unlink()
    with pytest.raises(OSError):
        read_cohort_evidence(path)


def test_existing_and_crash_left_destinations_are_never_adopted(tmp_path):
    src, plan, result, path = publish(tmp_path)
    before = (path/"bundle.json").read_bytes()
    with pytest.raises(FileExistsError):
        write_cohort_evidence(src, plan, result, path)
    assert (path/"bundle.json").read_bytes() == before
    crash = tmp_path/"crash-left"
    crash.mkdir()
    (crash/"private.log").write_text("keep")
    with pytest.raises(FileExistsError):
        write_cohort_evidence(src, plan, result, crash)
    assert (crash/"private.log").read_text() == "keep"


def test_payload_symlinks_hardlinks_and_directory_symlinks_refuse(tmp_path):
    _, _, _, path = publish(tmp_path)
    original = path/"source.pdb"
    external = tmp_path/"outside.pdb"
    external.write_bytes(original.read_bytes())
    original.unlink()
    original.symlink_to(external)
    with pytest.raises(ValueError, match="symlink"):
        read_cohort_evidence(path)
    original.unlink()
    original.hardlink_to(external)
    with pytest.raises(ValueError, match="single-linked"):
        read_cohort_evidence(path)
    alias = tmp_path/"alias"
    alias.symlink_to(path, target_is_directory=True)
    with pytest.raises(ValueError, match="regular directory"):
        read_cohort_evidence(alias)


def test_external_source_pin_and_flat_inventory_are_enforced(tmp_path):
    _, _, _, path = publish(tmp_path)
    with pytest.raises(ValueError, match="externally expected"):
        read_cohort_evidence(path, expected_source_sha256="0"*64)
    index = json.loads((path/"bundle.json").read_bytes())
    index["artifacts"]["../source.pdb"] = index["artifacts"].pop("source.pdb")
    (path/"bundle.json").write_bytes(encoded(index))
    with pytest.raises(ValueError, match="flat declared"):
        read_cohort_evidence(path)


def test_numeric_boolean_forgery_is_not_accepted_as_equal_python_values(tmp_path):
    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    first = result.attempts[0]
    forged_audit = replace(first.audits[0], accepted=1)
    forged = replace(result, attempts=(replace(first, audits=(forged_audit, *first.audits[1:])), *result.attempts[1:]))
    assert forged == result  # Python equality alone is insufficient for contracts.
    with pytest.raises(ValueError, match="replay"):
        write_cohort_evidence(src, plan, forged, tmp_path/"forged")
    assert not list(tmp_path.iterdir())


def test_staging_io_failure_never_commits_readiness_and_retains_written_diagnostics(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.evidence as module

    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    original_write = module._durable_write

    def fail_sidecar(path, data):
        if path.name == "candidate.dat":
            raise OSError("injected staging I/O failure")
        original_write(path, data)

    monkeypatch.setattr(module, "_durable_write", fail_sidecar)
    with pytest.raises(OSError, match="injected"):
        write_cohort_evidence(src, plan, result, tmp_path/"out")
    assert not (tmp_path/"out").exists()
    stages = list(tmp_path.glob(".protein-cohort-staging-*"))
    assert len(stages) == 1
    assert (stages[0]/"source.pdb").read_bytes() == src
    assert (stages[0]/"attempts.json").is_file()
    assert not (stages[0]/"bundle.json").exists()
    with pytest.raises(OSError):
        read_cohort_evidence(stages[0])


def test_competing_reservation_between_preflight_and_commit_is_preserved(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.evidence as module

    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    commit = module.commit_manifest_last

    def competitor(staging, destination):
        destination.mkdir()
        (destination/"user-owned.txt").write_text("keep")
        commit(staging, destination)

    monkeypatch.setattr(module, "commit_manifest_last", competitor)
    with pytest.raises(FileExistsError):
        write_cohort_evidence(src, plan, result, tmp_path/"race")
    assert (tmp_path/"race/user-owned.txt").read_text() == "keep"
    assert not (tmp_path/"race/bundle.json").exists()
    stages = list(tmp_path.glob(".protein-cohort-staging-*"))
    assert len(stages) == 1 and (stages[0]/"source.pdb").read_bytes() == src


def test_artifact_size_and_total_budgets_fail_before_staging(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.evidence as module

    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)
    monkeypatch.setattr(module, "MAX_ARTIFACT_BYTES", 1)
    with pytest.raises(ValueError, match="artifact budget"):
        write_cohort_evidence(src, plan, result, tmp_path/"too-big")
    assert list(tmp_path.iterdir()) == []
    monkeypatch.setattr(module, "MAX_ARTIFACT_BYTES", 64*1024*1024)
    monkeypatch.setattr(module, "MAX_TOTAL_BYTES", 1)
    with pytest.raises(ValueError, match="artifact budget"):
        write_cohort_evidence(src, plan, result, tmp_path/"total-too-big")
    assert list(tmp_path.iterdir()) == []


def test_unexpected_materialization_errors_are_not_relabelled_as_scientific_rejections(tmp_path, monkeypatch):
    import dvbfixer.protein_completion.evidence as module

    src = source()
    plan = plan_cohort(src, TARGETS)
    result = complete_cohort(src, plan)

    def failure(*args, **kwargs):
        raise module.SourceCompletionError("unexpected materializer failure")

    monkeypatch.setattr(module, "materialize_cohort", failure)
    with pytest.raises(ValueError, match="unexpected materializer"):
        write_cohort_evidence(src, plan, result, tmp_path/"unexpected")
    assert list(tmp_path.iterdir()) == []
