"""Small isolated native expression controls, not neural restoration evidence."""
import importlib.util
from pathlib import Path

import pytest
import torch
from rfd3.model.layers.blocks import PositionPairDistEmbedder

SPEC = importlib.util.spec_from_file_location("private_mps_extended", Path(__file__).with_name("mps_pairwise_extended.py"))
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


@pytest.mark.parametrize("device", ["cpu", "mps"])
@pytest.mark.parametrize("batch", [False, True])
@pytest.mark.parametrize("rows", [1, 4, 32])
def test_reference_tiling_matches_native_with_coincident_positions_and_mask(device, batch, rows):
    if device == "mps" and not torch.backends.mps.is_available():
        pytest.skip("requires MPS")
    torch.manual_seed(23)
    module = PositionPairDistEmbedder(16, embed_frame=False).eval().to(device)
    shape = (2, 11, 3) if batch else (11, 3)
    positions = torch.randn(shape).to(device)
    positions[..., 2, :] = positions[..., 1, :]
    mask = (torch.rand((*shape[:-2], 11, 11, 1)) > .4).to(device)
    with torch.no_grad():
        wanted = module(positions, mask)
        actual = HELPER.reference_rows(module, positions, mask, rows=rows)
        accumulator = torch.randn_like(wanted)
        expected = accumulator + wanted
        combined = HELPER.reference_rows(module, positions, mask, rows=rows, accumulate_into=accumulator)
        assert combined is accumulator
        torch.testing.assert_close(combined, expected, rtol=1e-5, atol=1e-6)
    torch.testing.assert_close(actual, wanted, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize("device", ["cpu", "mps"])
@pytest.mark.parametrize("batch", [False, True])
@pytest.mark.parametrize("rows", [1, 4, 32])
def test_pointwise_tiling_preserves_all_native_mlp_layers(device, batch, rows):
    if device == "mps" and not torch.backends.mps.is_available():
        pytest.skip("requires MPS")
    torch.manual_seed(29)
    module = torch.nn.Sequential(torch.nn.ReLU(), torch.nn.Linear(16, 16, bias=False),
                                 torch.nn.ReLU(), torch.nn.Linear(16, 16, bias=False),
                                 torch.nn.ReLU(), torch.nn.Linear(16, 16, bias=False)).eval().to(device)
    values = torch.randn((2, 11, 13, 16) if batch else (11, 13, 16)).to(device)
    with torch.no_grad():
        wanted = module(values)
        actual = HELPER.pointwise_rows(module, values, rows=rows)
    torch.testing.assert_close(actual, wanted, rtol=1e-5, atol=1e-6)


def test_unsupported_frame_embedding_and_cross_row_reduction_fail():
    module = PositionPairDistEmbedder(16, embed_frame=True).eval()
    with torch.no_grad():
        with pytest.raises(ValueError, match="non-frame"):
            HELPER.reference_rows(module, torch.zeros(2, 3), torch.ones(2, 2, 1, dtype=torch.bool))
        with pytest.raises(ValueError, match="pair axes"):
            HELPER.pointwise_rows(lambda x: x.sum(dim=-2), torch.ones(2, 2, 16))
