"""Isolated cache hooks preserve numerical outputs, inputs and random states."""
import importlib.util
import json
from pathlib import Path

import pytest
import torch

SPEC = importlib.util.spec_from_file_location("private_mps_cache", Path(__file__).with_name("mps_cache_boundary.py"))
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def test_explicit_boundary_keeps_values_parameters_and_rng(tmp_path):
    if not torch.backends.mps.is_available():
        pytest.skip("requires MPS")
    torch.manual_seed(37)
    module = torch.nn.Sequential(torch.nn.Linear(16, 32), torch.nn.ReLU(), torch.nn.Linear(32, 16)).eval().to("mps")
    values = torch.randn(65, 16).to("mps")
    log = tmp_path / "boundaries.jsonl"
    keys = tuple(module.state_dict())
    with torch.no_grad():
        wanted = module(values)
        unchanged = values.clone()
        cpu_rng, mps_rng = torch.get_rng_state().clone(), torch.mps.get_rng_state().clone()
        handles = HELPER.install_boundary(module, log, label="test-MLP")
        actual = module(values)
    for handle in handles:
        handle.remove()
    assert torch.equal(actual, wanted) and torch.equal(values, unchanged)
    assert torch.equal(cpu_rng, torch.get_rng_state()) and torch.equal(mps_rng, torch.mps.get_rng_state())
    assert tuple(module.state_dict()) == keys
    entries = [json.loads(line) for line in log.read_text().splitlines()]
    assert [entry["phase"] for entry in entries] == ["pre", "post"]
    assert all(entry["after"]["driver"] >= entry["after"]["live"] for entry in entries)


def test_cpu_module_is_not_touched(tmp_path):
    module = torch.nn.Linear(2, 2).eval()
    log = tmp_path / "unused.jsonl"
    handles = HELPER.install_boundary(module, log, label="CPU")
    module(torch.ones(1, 2))
    for handle in handles:
        handle.remove()
    assert not log.exists()
