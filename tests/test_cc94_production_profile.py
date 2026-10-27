"""Checkpoint migration, physical-chain admission, and no-fallback contracts."""
import importlib
import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.contract import BackendOption
from dvbfixer.model.diffusion.preflight import DIFFUSION_PROFILES
from dvbfixer.model.diffusion.protein_context import (
    ProteinContextError,
    independent_gap_requests,
    prepare_cc94_request,
    source_context_position_features,
    validate_cc94_request,
    validate_context_capability,
)
from dvbfixer.model.diffusion.runtime import DiffusionRuntimeError, resolve_diffusion_runtime
from dvbfixer.model.diffusion.scope import assess_diffusion_scope
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def cli(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "tests"))
    return importlib.import_module("test_diffusion_cli")


@pytest.mark.parametrize("profile", sorted(cc94.PROFILES))
def test_cc94_prepared_local_request_is_stable_with_native_relative_phase(tmp_path, cli, profile):
    source = cli._gap_input(tmp_path)
    request = build_cli_diffusion_request(source, {"C": "SNRFSGSKSGNTA"}, seeds=(7,), profile=profile)
    request = prepare_cc94_request(request, source.read_bytes(), "auto")
    assert BackendOption("cc89_chain_phase", "native") in request.backend_options
    assert BackendOption("cc89_chain_phase", "neutral-v1") not in request.backend_options
    assert prepare_cc94_request(request, source.read_bytes(), "auto") == request
    validate_cc94_request(request)
    validate_context_capability(request, "protpardelle")
    assert assess_diffusion_scope(request, source.read_bytes()).supported
    conflict = replace(request, backend_options=tuple(
        BackendOption(o.name, "neutral-v1") if o.name == "cc89_chain_phase" else o for o in request.backend_options))
    with pytest.raises(ProteinContextError):
        validate_cc94_request(conflict)


def test_full_9cqd_seven_gaps_keep_independent_observed_only_context():
    source = ROOT / "tests/fixtures/protein_context/9cqd_split.pdb"
    provenance = json.loads((ROOT / "docs/research/9cqd-protenix-cc89-results.json").read_text())
    sequences = {t["chain"]: t["sequence"] for t in provenance["target_sequences"]}
    request = build_cli_diffusion_request(source, sequences, seeds=(7,), profile=cc94.CUDA_PROFILE)
    prepared = prepare_cc94_request(request, source.read_bytes(), "auto")
    assert len(prepared.gaps) == 7
    assert len(prepared.gap_contexts) == 7
    assert assess_diffusion_scope(prepared, source.read_bytes()).supported
    for single in independent_gap_requests(prepared):
        validate_cc94_request(single)
        assert {item.atom for item in single.protein_context.fixed_conditioning} <= set(prepared.fixed_atoms)
        assert not set(single.generated_atoms) & set(prepared.fixed_atoms)
        positions, chains = source_context_position_features(single)
        assert len(positions) == len(chains) <= 512


@pytest.mark.parametrize("host,machine,profile", [
    ("linux", "x86_64", cc94.CUDA_PROFILE), ("darwin", "arm64", cc94.MPS_PROFILE),
])
def test_new_requests_require_cc94_but_old_profiles_remain_explicit(monkeypatch, host, machine, profile):
    monkeypatch.setattr("dvbfixer.model.diffusion.runtime.sys.platform", host)
    monkeypatch.setattr("dvbfixer.model.diffusion.runtime.platform.machine", lambda: machine)
    monkeypatch.setattr("dvbfixer.model.diffusion.runtime._installed_runner", lambda name: "/installed/" + name)
    assert resolve_diffusion_runtime("protpardelle").profile == profile
    legacy = "protpardelle-1c-cuda" if host == "linux" else "protpardelle-1c-mps-local-v1"
    assert resolve_diffusion_runtime("protpardelle", profile=legacy).profile == legacy
    monkeypatch.setattr("dvbfixer.model.diffusion.runtime._installed_runner",
                        lambda name: "/old/" + name if "cc94" not in name else None)
    with pytest.raises(DiffusionRuntimeError, match="cc94"):
        resolve_diffusion_runtime("protpardelle")


def test_cc94_evidence_does_not_inherit_cc89_mps_or_cohort_acceptance():
    for profile in cc94.PROFILES:
        metadata = DIFFUSION_PROFILES[profile]
        assert metadata.expected_checkpoint_sha256 == cc94.CHECKPOINT_SHA256
        assert "source-disjoint-24-of-24" not in metadata.evidence_labels
        assert "apple-portability-8-of-8" not in metadata.evidence_labels
    assert "unverified" in DIFFUSION_PROFILES[cc94.MPS_PROFILE].status


def test_independent_artifact_reader_does_not_import_scientific_refinement(tmp_path):
    import hashlib
    import subprocess

    path = tmp_path / "source.pdb"
    path.write_text("END\n")
    script = """import sys
from pathlib import Path
from dvbfixer.model.diffusion.contract import ArtifactReference
from dvbfixer.model.diffusion.independent_gaps import _read
assert _read(Path(sys.argv[1]),ArtifactReference('source.pdb',sys.argv[2]))=='END\\n'
assert 'dvbfixer.model.diffusion.runner_refinement' not in sys.modules
assert 'pdbfixer' not in sys.modules
"""
    subprocess.run([sys.executable, "-c", script, str(tmp_path), hashlib.sha256(path.read_bytes()).hexdigest()], check=True)


def test_cc94_denoiser_probe_checks_without_neutralizing_chains():
    torch = pytest.importorskip("torch")
    path = ROOT / "deploy/protpardelle-1c/checkpoint_gap_smoke.py"
    spec = importlib.util.spec_from_file_location("cc94_test_adapter", path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    probe = adapter._RelativePositionProbe((1, 2, 202), (0, 0, 1), "cpu")
    features = {"residue_index": torch.tensor([[1, 2, 202]]), "chain_index": torch.tensor([[0, 0, 1]])}
    before = features["chain_index"].clone()
    probe(None, (), features)
    assert torch.equal(features["chain_index"], before)
    assert probe.calls == 1
    with pytest.raises(RuntimeError, match="physical position"):
        probe(None, (), {**features, "chain_index": torch.zeros_like(before)})


def test_new_runner_rejects_cc89_weights_before_loading_model(tmp_path, monkeypatch):
    from types import SimpleNamespace

    path = ROOT / "deploy/protpardelle-1c/cc94_runner.py"
    spec = importlib.util.spec_from_file_location("cc94_test_runner", path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    weights = tmp_path / "weights" / cc94.CHECKPOINT_NAME
    weights.parent.mkdir()
    weights.write_bytes(b"wrong weights")
    monkeypatch.setattr(runner, "_sha", lambda p: (
        DIFFUSION_PROFILES["protpardelle-1c-cuda"].expected_checkpoint_sha256
        if p == weights else cc94.CONFIG_SHA256))
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **kw: (_ for _ in ()).throw(ValueError("no source")))
    monkeypatch.setattr(runner.subprocess, "run", lambda *a, **kw: None)
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(
        __version__="2.6.0+cu124", cuda=SimpleNamespace(is_available=lambda: False)))
    report = runner.preflight_report(cc94.CUDA_PROFILE, weights, Path(sys.executable))
    assert not report.passed
    assert any(i.code.value == "checkpoint-digest-mismatch" for i in report.issues)
