"""Native gated/RMSNorm Transition is pointwise across query/key pair axes."""
import importlib.util
from pathlib import Path

import pytest
import torch
from rfd3.model.layers.layer_utils import Transition

SPEC = importlib.util.spec_from_file_location("private_transition_rows", Path(__file__).with_name("mps_pairwise_extended.py"))
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


@pytest.mark.parametrize("device", ["cpu", "mps"])
@pytest.mark.parametrize("batch", [False, True])
@pytest.mark.parametrize("rows", [1, 4, 32])
def test_native_RMSNorm_gated_transition_keeps_all_last_axis_features(device, batch, rows):
    if device == "mps" and not torch.backends.mps.is_available():
        pytest.skip("requires MPS")
    torch.manual_seed(41)
    module = Transition(n=4, c=128).eval().to(device)
    values = torch.randn((2, 11, 13, 128) if batch else (11, 13, 128)).to(device)
    keys = tuple(module.state_dict())
    with torch.no_grad():
        wanted = module(values)
        actual = HELPER.pointwise_rows(module, values, rows=rows)
    torch.testing.assert_close(actual, wanted, rtol=1e-5, atol=1e-6)
    assert tuple(module.state_dict()) == keys
