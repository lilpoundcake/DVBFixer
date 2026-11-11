"""Native private MPS reinjection controls; no checkpoint sampling claims."""
import numpy as np
import pytest
import torch
from mps_research_smoke import MPSReinjector

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_mps_reinjection_preserves_observations_and_generated_rigid_geometry():
    axis = tuple(AtomIdentity("A", "1", "B", name) for name in ("N", "CA", "C", "CB"))
    target = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    fixed = dict(zip(axis[:3], target[:3]))
    callback = MPSReinjector(axis, fixed)
    rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
    mobile = torch.tensor((target @ rotation + 11)[None], device="mps", dtype=torch.float32)
    result = callback(mobile, 0, 1)
    assert result.device.type == "mps" and result.dtype == torch.float32
    assert torch.equal(result[0, :3], torch.tensor(target[:3], device="mps", dtype=torch.float32))
    assert np.allclose(result.cpu().numpy(), target[None], atol=1e-5)
    assert callback.steps[0]["post_overwrite_max_angstrom"] == 0
    assert callback.last_pre_overwrite is not None


def test_private_mps_hook_refuses_wrong_device_and_unsupported_weighting():
    axis = tuple(AtomIdentity("a", "7", "", name) for name in ("N", "CA", "C"))
    fixed = dict(zip(axis, np.eye(3)))
    with pytest.raises(ValueError, match="uniform"):
        MPSReinjector(axis, fixed, anchor_weights={a:float(i+1) for i,a in enumerate(axis)})
    callback = MPSReinjector(axis, fixed)
    with pytest.raises(ValueError, match="MPS float32"):
        callback(torch.eye(3), 0, 1)
