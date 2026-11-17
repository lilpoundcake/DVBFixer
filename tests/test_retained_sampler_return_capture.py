import json
from types import SimpleNamespace

import numpy as np
import pytest
from scripts import retained_sampler_return_capture as capture

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_passive_capture_returns_same_tensor_and_restores_methods_on_downstream_failure(monkeypatch):
    tensor, request = object(), object()
    events = []

    class Model:
        dvbfixer_step_callback = SimpleNamespace(steps=["retained"])

        def sample_diffusion(self, **kwargs):
            events.append(("sample", kwargs))
            return tensor

    axis_result = (("a", "A"), ("ALA", "ALA"), ("C", "C"))
    legacy = SimpleNamespace(_model_atom_axis=lambda *args: axis_result)
    old_axis, old_sample = legacy._model_atom_axis, Model.sample_diffusion

    def record(output, value, **kwargs):
        assert value is tensor and kwargs["axis"] == axis_result[0]
        events.append(("capture", kwargs["steps"]))
        return value

    monkeypatch.setattr(capture, "capture_return", record)
    with pytest.raises(RuntimeError, match="downstream"), capture.passive_capture(
        legacy, Model, output="capture", request=request, request_path="request", source="source"
    ):
        assert legacy._model_atom_axis("array", request) is axis_result
        assert Model().sample_diffusion(N_sample=1) is tensor
        events.append(("downstream",))
        raise RuntimeError("downstream head failure")
    assert [row[0] for row in events] == ["sample", "capture", "downstream"]
    assert legacy._model_atom_axis is old_axis and Model.sample_diffusion is old_sample


def test_capture_refuses_partial_trace_before_writing(tmp_path):
    with pytest.raises(ValueError, match="all200"):
        capture.capture_return(tmp_path / "new", None, request=None, request_path=None,
                               source=None, axis=(), names=(), elements=(), steps=[])
    assert not (tmp_path / "new").exists()


def test_capture_commits_full_case_sensitive_axis_last_without_tensor_mutation(tmp_path, monkeypatch):
    source, request_path = tmp_path / "source.pdb", tmp_path / "request.json"
    source.write_text("synthetic-source\n")
    request_path.write_text("synthetic-request\n")
    fixed, generated = AtomIdentity("a", "82", "A", "CA"), AtomIdentity("A", "82", "A", "CA")
    request = SimpleNamespace(normalized_pdb=SimpleNamespace(sha256=capture.digest(source)),
        protein_context=SimpleNamespace(fixed_conditioning=(SimpleNamespace(atom=fixed),)),
        generated_atoms=(generated,), seeds=(11,))
    xyz = np.arange(6, dtype=np.float32).reshape(1, 2, 3)
    original = xyz.copy()
    tensor = SimpleNamespace(device=SimpleNamespace(type="mps"), dtype="torch.float32",
        detach=lambda: SimpleNamespace(cpu=lambda: SimpleNamespace(numpy=lambda: xyz)))
    monkeypatch.setattr(capture, "materialize_candidate_pdb", lambda *args, **kw: "synthetic-raw\n")
    steps = [{"step_index": i, "step_count": 200, "post_overwrite_max_angstrom": 0} for i in range(200)]
    output = tmp_path / "capture"
    args = dict(request=request, request_path=request_path, source=source,
                axis=(fixed, generated), names=("ALA", "ALA"), elements=("C", "C"), steps=steps)
    assert capture.capture_return(output, tensor, **args) is tensor
    assert np.array_equal(xyz, original)
    manifest = json.loads((output / "capture.json").read_text())
    assert not manifest["native_full_inference_success"] and not manifest["scientific_acceptance"]
    for name, expected in manifest["files"].items():
        assert capture.digest(output / name) == expected
    axis = json.loads((output / "native-axis.json").read_text())["atoms"]
    assert [atom["chain"] for atom in axis] == ["a", "A"]
    with pytest.raises(FileExistsError):
        capture.capture_return(output, tensor, **args)
