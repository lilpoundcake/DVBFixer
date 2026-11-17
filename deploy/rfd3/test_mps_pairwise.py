"""Isolated native runtime checks; random tensors are not protein proposals."""
import importlib.util
from pathlib import Path

import pytest
import torch
from rfd3.model.layers.blocks import SinusoidalDistEmbed

SPEC = importlib.util.spec_from_file_location("private_mps_pairwise", Path(__file__).with_name("mps_pairwise.py"))
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


@pytest.mark.parametrize("device", ["cpu", "mps"])
@pytest.mark.parametrize("batch", [False, True])
@pytest.mark.parametrize("rows", [1, 4, 32])
def test_tiled_expression_matches_complete_native_dense_result(device, batch, rows):
    if device == "mps" and not torch.backends.mps.is_available():
        pytest.skip("requires native MPS")
    torch.manual_seed(17)
    module = SinusoidalDistEmbed(c_atompair=16).eval().to(device)
    shape = (2, 11, 3) if batch else (11, 3)
    positions = torch.randn(shape).to(device)
    mask = (torch.rand((*shape[:-2], 11, 11, 1)) > .4).to(device)
    with torch.no_grad():
        wanted = module(positions, mask)
        actual = HELPER.sinusoidal_rows(module, positions, mask, rows=rows)
    torch.testing.assert_close(actual, wanted, rtol=1e-5, atol=1e-6)
    assert torch.equal(actual[~mask.expand_as(actual)], torch.zeros_like(actual[~mask.expand_as(actual)]))
    assert actual.dtype == wanted.dtype and actual.shape == wanted.shape


@pytest.mark.parametrize("rows", [0, -1, True, 1.5])
def test_invalid_tile_size_refused(rows):
    module = SinusoidalDistEmbed(c_atompair=16).eval()
    with torch.no_grad(), pytest.raises(ValueError, match="positive integer"):
        HELPER.sinusoidal_rows(module, torch.zeros(2, 3), torch.ones(2, 2, 1, dtype=torch.bool), rows=rows)


def test_training_and_gradients_are_not_silently_modified():
    module = SinusoidalDistEmbed(c_atompair=16)
    with torch.no_grad(), pytest.raises(ValueError, match="inference"):
        HELPER.sinusoidal_rows(module, torch.zeros(2, 3), torch.ones(2, 2, 1, dtype=torch.bool))
    module.eval()
    with pytest.raises(ValueError, match="inference"):
        HELPER.sinusoidal_rows(module, torch.zeros(2, 3), torch.ones(2, 2, 1, dtype=torch.bool))
