"""Synthetic evidence tests are not checkpoint execution or qualification."""

import hashlib
import json
import os
import subprocess
import sys
import zipfile
from dataclasses import asdict, replace
from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.atom_completion_evidence import (
    _arrays,
    attest_completion_run,
    audit_completion_evidence,
)
from dvbfixer.model.diffusion.atom_completion_materialize import materialize_atom_completion
from dvbfixer.model.diffusion.atom_completion_publication import (
    CompletionBundleIndex,
    completion_dat_record,
    publish_completion_bundle,
    read_completion_bundle,
)
from dvbfixer.model.diffusion.atom_completion_refinement import COMPLETION_REFINEMENT_REVISION
from dvbfixer.model.diffusion.atom_completion_sampling import completion_engine_view
from dvbfixer.model.diffusion.atom_completion_trace import (
    CompletionSamplerTrace,
    completion_conditioning,
)
from dvbfixer.model.diffusion.atom_completion_validation import _coordinates
from dvbfixer.model.diffusion.boundary_refinement import BOUNDARY_REFINEMENT_REVISION
from dvbfixer.model.diffusion.contract import ArtifactReference, DiffusionContractError
from dvbfixer.model.diffusion.preflight import RunnerPreflightFacts, RunnerPreflightReport
from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER
from tests.test_diffusion_atom_completion_sampling import sampling_control as sampling_control
from tests.test_diffusion_atom_completion_validation import completion_control as completion_control


def _json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False))


def _reference(root, name):
    return ArtifactReference(name, hashlib.sha256((root / name).read_bytes()).hexdigest())


@pytest.fixture
def synthetic_run(sampling_control, tmp_path):
    source, reference, request = sampling_control
    root = tmp_path / "synthetic-run-not-sampling"
    root.mkdir()
    view = completion_engine_view(request, source)
    points = _coordinates(reference.decode().splitlines())
    source_points = _coordinates(source.decode().splitlines())
    raw = materialize_atom_completion(source, request, request.generated_atoms, [points[atom] for atom in request.generated_atoms])
    (root / "source.pdb").write_bytes(source)
    (root / "request.json").write_text(request.to_json())
    (root / "motif.pdb").write_text(view.motif_text)
    (root / "raw.pdb").write_bytes(raw)
    (root / "refined.pdb").write_bytes(raw)
    _json(root / "engine-view.json", asdict(view))
    native = np.zeros((72, 37, 3), dtype=np.float64)
    mask = np.zeros((72, 37))
    observed = np.zeros((72, 37, 3))
    fixed = set(view.represented_fixed_atoms)
    for atom, index in zip(view.atom_axis, view.atom_target_indices, strict=True):
        column = ATOM37_ORDER.index(atom.atom_name)
        native[index, column] = points[atom]
        if atom in fixed:
            mask[index, column] = 1
            observed[index, column] = source_points[atom]
    np.savez_compressed(root / "native.npz", atom37=native)
    np.savez_compressed(root / "reader.npz", atom_mask=mask, atom_positions=observed, residue_index=np.arange(1, 73))
    facts = RunnerPreflightFacts(platform="darwin", architecture="arm64", python_version="3.12",
        engine_revision=cc94.ENGINE_REVISION, patch_identity="cc94-native-relative-v1;source-physical-v1;apple-patch-v1",
        environment_identity=cc94.MPS_ENVIRONMENT, checkpoint_sha256=cc94.CHECKPOINT_SHA256,
        framework_version="2.6.0", accelerator_available=True, effective_device="mps", fallback_disabled=True,
        refinement_platform="OpenMM CPU")
    (root / "preflight.json").write_text(RunnerPreflightReport(cc94.MPS_PROFILE, 4, 4, facts).to_json())
    adapter = b"# synthetic unit-test evidence, NOT executed checkpoint code\n"
    progress = {"schema_version": 1, "kind": "private-cc94-atom-completion-trace", "status": "accepted-private",
        "seed": 7, "device": "mps", "ode_updates": 500, "step_scale": 1.0, "s_churn": 0.0, "chain_phase": "native",
        "checkpoint_sha256": cc94.CHECKPOINT_SHA256, "config_sha256": cc94.CONFIG_SHA256, "engine_revision": cc94.ENGINE_REVISION,
        "checkpoint_sampling_executed": True, "checkpoint_sampling_started": True,
        "published": False, "scientific_qualification": False, "public_profile_enabled": False,
        "conditioning": completion_conditioning(view.fixed_positions), "adapter_sha256": hashlib.sha256(adapter).hexdigest()}
    for role, name in (("request", "request.json"), ("source", "source.pdb"), ("motif", "motif.pdb"),
                       ("native", "native.npz"), ("raw", "raw.pdb")):
        progress[role + "_sha256"] = _reference(root, name).sha256
    _json(root / "trace.json", progress)
    _json(root / "refinement.json", {"status": "accepted-private", "revision": COMPLETION_REFINEMENT_REVISION,
        "kernel_revision": BOUNDARY_REFINEMENT_REVISION, "seed": 7, "platform": "CPU", "published": False,
        "scientific_qualification": False, "source_sha256": _reference(root, "source.pdb").sha256,
        "raw_sha256": _reference(root, "raw.pdb").sha256, "refined_sha256": _reference(root, "refined.pdb").sha256,
        "initial_energy_kj_mol": 1.0, "final_energy_kj_mol": 0.0})
    trace = attest_completion_run(root, adapter)
    return root, trace


def test_trace_roundtrip_independent_audit_and_dat_ownership(synthetic_run):
    root, trace = synthetic_run
    assert CompletionSamplerTrace.from_json(trace.to_json()) == trace
    audit = audit_completion_evidence(trace, root)
    dat = completion_dat_record(audit.request, audit.payloads["source"], audit.payloads["refined"])
    assert len(dat.added_atoms) == 4
    assert dat.added_keys() == {("A", "27", "", name) for name in ("CG", "CD", "CE", "NZ")}
    assert not dat.added_keys() & {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in audit.request.masks.fixed_atoms}
    with pytest.raises(FileExistsError):
        attest_completion_run(root, (root / trace.adapter_source.path).read_bytes())


def test_manifest_last_bundle_roundtrip_and_exclusive_reservation(synthetic_run, tmp_path):
    from dvbfixer.model.diffusion.publication import PublicationError, read_published_bundle

    root, trace = synthetic_run
    destination = tmp_path / "experimental-bundle"
    publish_completion_bundle(trace, root, destination)
    assert read_completion_bundle(destination).refined_validation.summary.passed
    record = DatRecord.load(destination / "candidate.dat")
    assert len(record.added_atoms) == 4
    before = (destination / "bundle.json").read_bytes()
    with pytest.raises(FileExistsError):
        publish_completion_bundle(trace, root, destination)
    assert (destination / "bundle.json").read_bytes() == before
    assert (root / "refined.pdb").exists()
    with pytest.raises(PublicationError, match="unsupported"):
        read_published_bundle(destination)  # Never masquerade as a public gap bundle.


@pytest.mark.parametrize("change", [{"ode_updates": 499}, {"ode_updates": True}, {"fallback_disabled": False},
                                     {"step_scale": True}, {"chain_phase": "neutral-v1"}, {"seed": True},
                                     {"checkpoint_sha256": "0" * 64}, {"refinement_mode": "observed-mobile"}])
def test_forged_trace_policy_rejected(synthetic_run, change):
    with pytest.raises(DiffusionContractError):
        replace(synthetic_run[1], **change)


@pytest.mark.parametrize("role,name,change", [
    ("progress", "trace.json", {"ode_updates": 499}),
    ("progress", "trace.json", {"scientific_qualification": True}),
    ("refinement", "refinement.json", {"platform": "CUDA"}),
    ("refinement", "refinement.json", {"seed": 19}),
])
def test_report_forgery_rejected_even_with_updated_digest(synthetic_run, role, name, change):
    root, trace = synthetic_run
    data = json.loads((root / name).read_text())
    data.update(change)
    _json(root / name, data)
    trace = replace(trace, **{role: _reference(root, name)})
    with pytest.raises(DiffusionContractError):
        audit_completion_evidence(trace, root)


def test_native_or_reader_substitution_rejected(synthetic_run):
    root, trace = synthetic_run
    native = np.load(root / "native.npz")["atom37"]
    native[26, ATOM37_ORDER.index("NZ"), 0] += 1
    np.savez_compressed(root / "native.npz", atom37=native)
    native_ref = _reference(root, "native.npz")
    progress = json.loads((root / "trace.json").read_text())
    progress["native_sha256"] = native_ref.sha256
    _json(root / "trace.json", progress)
    with pytest.raises(DiffusionContractError, match="native sample"):
        audit_completion_evidence(replace(trace, native=native_ref, progress=_reference(root, "trace.json")), root)


def test_symlink_hardlink_and_payload_corruption_fail_closed(synthetic_run, tmp_path):
    root, trace = synthetic_run
    (root / "source-link.pdb").symlink_to(root / "source.pdb")
    with pytest.raises(ValueError, match="symlink"):
        audit_completion_evidence(replace(trace, source=replace(trace.source, path="source-link.pdb")), root)
    hard = tmp_path / "external-hardlink.pdb"
    os.link(root / "source.pdb", hard)
    with pytest.raises(ValueError, match="single-linked"):
        audit_completion_evidence(trace, root)


def test_incomplete_or_crash_left_bundle_never_ready(synthetic_run, tmp_path, monkeypatch):
    import dvbfixer.model.diffusion.atom_completion_publication as publication

    root, trace = synthetic_run
    destination = tmp_path / "interrupted"

    def interrupted(staging, target):
        target.mkdir()
        assert not (target / "bundle.json").exists()
        raise OSError("declared interrupted commit")

    monkeypatch.setattr(publication, "commit_manifest_last", interrupted)
    with pytest.raises(OSError, match="interrupted"):
        publication.publish_completion_bundle(trace, root, destination)
    with pytest.raises(FileNotFoundError):
        read_completion_bundle(destination)
    with pytest.raises(FileExistsError):
        publication.publish_completion_bundle(trace, root, destination)
    assert list(tmp_path.glob(".completion-staging-*"))  # Retain failed staging.


def test_corrupt_dat_and_false_public_qualification_rejected(synthetic_run, tmp_path):
    root, trace = synthetic_run
    destination = publish_completion_bundle(trace, root, tmp_path / "bundle")
    index = CompletionBundleIndex.from_json((destination / "bundle.json").read_text())
    dat = DatRecord.load(destination / "candidate.dat")
    dat.added_atoms[0]["atom"] = "CA"
    dat.save(destination / "candidate.dat", verbose=False)
    index = replace(index, dat=_reference(destination, "candidate.dat"))
    (destination / "bundle.json").write_text(index.to_json())
    with pytest.raises(DiffusionContractError, match="ownership"):
        read_completion_bundle(destination)
    with pytest.raises(DiffusionContractError, match="qualification"):
        replace(index, scientific_qualification=True)


def test_numpy_header_bomb_and_object_arrays_refused_before_load():
    header = BytesIO()
    np.lib.format.write_array_header_1_0(header, {"shape": (10**12, 37, 3), "fortran_order": False, "descr": "<f8"})
    payload = BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr("atom37.npy", header.getvalue())
    with pytest.raises(DiffusionContractError, match="bounds"):
        _arrays(payload.getvalue(), {"atom37"})
    objects = BytesIO()
    np.savez(objects, atom37=np.array([object()], dtype=object))
    with pytest.raises(DiffusionContractError, match="bounds"):
        _arrays(objects.getvalue(), {"atom37"})


def test_trace_contract_remains_stdlib_only():
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, "-S", "-c", """
import sys
from dvbfixer.model.diffusion.atom_completion_trace import CompletionSamplerTrace
assert not any(module in sys.modules for module in ('numpy', 'openmm', 'torch', 'Bio'))
"""], env=dict(os.environ, PYTHONPATH=str(root / "src")), check=True, capture_output=True, text=True)


def test_reader_mask_forgery_rejected_even_when_rehashed(synthetic_run):
    root, trace = synthetic_run
    with np.load(root / "reader.npz") as archive:
        values = {name: archive[name] for name in archive.files}
    values["atom_mask"][26, ATOM37_ORDER.index("CG")] = 1
    np.savez_compressed(root / "reader.npz", **values)
    with pytest.raises(DiffusionContractError, match="masks"):
        audit_completion_evidence(replace(trace, reader=_reference(root, "reader.npz")), root)


def test_dat_keeps_case_insertion_and_whole_gap_ownership_exact():
    from dvbfixer.model.diffusion.atom_completion import plan_request_atom_completion
    from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
    from tests.test_diffusion_atom_completion import _request

    source, legacy = _request()
    request = AtomCompletionRequest(1, legacy.normalized_pdb, legacy.target_sequences, legacy.sequence_placements,
                                    legacy.gaps, plan_request_atom_completion(legacy, source), legacy.seeds, "data-only-control")
    candidate = materialize_atom_completion(source, request, request.generated_atoms, [(1, 2, 3)] * len(request.generated_atoms))
    dat = completion_dat_record(request, source, candidate)
    assert {("D", "82", "", "NZ"), ("D", "82", "A", "NZ"), ("d", "82", "", "NZ")} <= dat.added_keys()
    assert len(dat.added_atoms) == 16  # Three partial LYS plus one whole GLY gap.
    assert len(dat.residue_summary) == 4
