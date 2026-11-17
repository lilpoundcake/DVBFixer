import json
from types import SimpleNamespace

import numpy as np
import pytest
from scripts import run_retained_loop_context_sampling as adapter

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_private_wiring_selects_only_represented_atoms_and_restores_legacy_module(monkeypatch):
    selected, outside = AtomIdentity("a", "1", "", "CA"), AtomIdentity("A", "1", "", "CA")
    coords = {selected: np.array([1., 2., 3.]), outside: np.array([4., 5., 6.])}
    request = SimpleNamespace(fixed_atoms=(selected, outside), protein_context=SimpleNamespace(
        fixed_conditioning=(SimpleNamespace(atom=selected, coordinate_angstrom=(1., 2., 3.)),), preserved_only=(outside,)))
    calls = []

    def reader(path, atoms):
        calls.append(tuple(atoms))
        return {a: coords[a] for a in atoms}

    class Callback:
        def __init__(self, axis, coordinates, *, alignment_atoms):
            assert axis == (selected,)
            assert set(coordinates) == {selected}
            assert alignment_atoms == (selected,)

    module = SimpleNamespace(_read_fixed_coordinates=reader, MPSReinjector=Callback)
    monkeypatch.setattr(adapter, "context_evidence", lambda *args: SimpleNamespace(alignment_atoms=(selected,)))
    with pytest.raises(RuntimeError, match="test interruption"), adapter.context_wiring(module, request):
        values = module._read_fixed_coordinates("source.pdb", request.fixed_atoms)
        module.MPSReinjector((selected,), values)
        raise RuntimeError("test interruption")
    assert calls == [(selected,)]
    assert module._read_fixed_coordinates is reader and module.MPSReinjector is Callback


def test_observation_partition_cannot_be_changed_by_the_adapter():
    atom = AtomIdentity("A", "1", "", "CA")
    request = SimpleNamespace(fixed_atoms=(atom,), protein_context=SimpleNamespace(
        fixed_conditioning=(SimpleNamespace(atom=atom),), preserved_only=(atom,)))
    with pytest.raises(ValueError, match="partition"), adapter.context_wiring(SimpleNamespace(), request):
        raise AssertionError("invalid partition must not execute")


def test_launcher_pins_torch_layernorm_and_project_namespace(tmp_path, monkeypatch):
    (tmp_path / "plan.json").write_text(json.dumps({"bindings": {}, "jobs": [{"command": ["unused"]}]}))
    monkeypatch.setenv("LAYERNORM_TYPE", "cuda")

    def launch(command, *, cwd, env, check):
        assert env["LAYERNORM_TYPE"] == "torch"
        assert env["PYTORCH_ENABLE_MPS_FALLBACK"] == "0"
        assert env["PYTHONPATH"].split(adapter.os.pathsep) == [str(adapter.PROJECT), str(adapter.PROJECT / "src")]
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(adapter.subprocess, "run", launch)
    adapter.run(tmp_path)
    assert (tmp_path / "completed.json").is_file()
