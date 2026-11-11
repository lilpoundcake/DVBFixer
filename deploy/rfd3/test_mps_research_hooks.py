"""Explicit isolated-runtime controls, not checkpoint/scientific qualification."""
import importlib

import numpy as np
import pytest
import torch
from biotite.structure import AtomArray
from rfd3.utils.inference import set_indices


def test_source_component_does_not_truncate_four_digit_residue_identity():
    array = AtomArray(2)
    array.chain_id[:] = "a"
    array.res_id[:] = 1
    array.coord[:] = np.array([[1, 2, 3], [4, 5, 6]])
    result = set_indices(array, "J", 1, 1, "a1037")
    assert list(result.src_component) == ["a1037", "a1037"]
    assert list(result.chain_id) == ["J", "J"]
    assert np.array_equal(result.coord, [[1, 2, 3], [4, 5, 6]])


def test_explicit_mps_seed_contract_changes_noise_and_repeats_exactly(tmp_path, monkeypatch):
    monkeypatch.setenv("DVBFIXER_RFD3_TRANSPORT_OUTPUT", str(tmp_path / "transport"))
    entry = importlib.import_module("rfd3_mps_seeded_research")
    monkeypatch.setattr(entry.transport, "traced", lambda *_args, **_kwargs: torch.randn(32, device="mps").cpu().numpy())
    monkeypatch.setenv("DVBFIXER_RFD3_MODEL_SEED", "7")
    first = entry.reseeded(None, {})
    repeat = entry.reseeded(None, {})
    monkeypatch.setenv("DVBFIXER_RFD3_MODEL_SEED", "42")
    other = entry.reseeded(None, {})
    assert np.array_equal(first, repeat)
    assert not np.array_equal(first, other)
    monkeypatch.setenv("DVBFIXER_RFD3_MODEL_SEED", "211")
    with pytest.raises(ValueError, match="development"):
        entry.reseeded(None, {})
