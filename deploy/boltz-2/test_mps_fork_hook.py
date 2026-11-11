"""Real fork sampler, fake score model: transport controls, NOT sampling evidence."""
import numpy as np
import pytest
import torch
from boltz.model.modules.diffusionv2 import AtomDiffusion
from mps_research_reinjection import MPSFixedAtomCallback
from torch import nn

from dvbfixer.model.diffusion.contract import AtomIdentity


class Probe(AtomDiffusion):
    def __init__(self, device):
        nn.Module.__init__(self)
        self.probe = nn.Parameter(torch.zeros((), device=device))
        self.sigma_min, self.sigma_max, self.sigma_data = .0004, 1., 1.
        self.rho, self.gamma_0, self.gamma_min = 7., 0., 1.
        self.noise_scale, self.step_scale = 1., 1.
        self.step_scale_random = None
        self.num_sampling_steps = 3
        self.alignment_reverse_diff = False

    @property
    def device(self):
        return self.probe.device

    def preconditioned_network_forward(self, noised_atom_coords, sigma, network_condition_kwargs):
        return torch.zeros_like(noised_atom_coords)


@pytest.mark.parametrize('device', ['cpu', 'mps'])
def test_actual_sampler_calls_hook_after_every_update_and_returns_exact_fixed(device, monkeypatch):
    monkeypatch.setenv('BOLTZ_DENOISER_CACHE', 'off')
    model = Probe(device)
    calls = []
    target = torch.tensor([[1., 2., 3.], [2., 1., 4.]], device=device)
    def callback(x, index, count):
        calls.append((index, count))
        x = x.clone()
        x[:, :2] = target
        return x
    result = model.sample(torch.ones((1, 4), device=device), step_callback=callback)['sample_atom_coords']
    assert calls == [(0, 3), (1, 3), (2, 3)]
    assert torch.equal(result[:, :2], target[None])
    assert torch.isfinite(result).all()


@pytest.mark.parametrize('fault', ['shape', 'dtype', 'finite', 'type'])
def test_hook_refuses_invalid_returns(fault, monkeypatch):
    monkeypatch.setenv('BOLTZ_DENOISER_CACHE', 'off')
    def bad(x, *_):
        return {'shape': lambda: x[:, :1], 'dtype': lambda: x.double(),
                'finite': lambda: x * float('nan'), 'type': lambda: None}[fault]()
    with pytest.raises((ValueError, TypeError)):
        Probe('cpu').sample(torch.ones((1, 4)), step_callback=bad)


def test_hook_refuses_steering(monkeypatch):
    monkeypatch.setenv('BOLTZ_DENOISER_CACHE', 'off')
    with pytest.raises(ValueError, match='steering'):
        Probe('cpu').sample(torch.ones((1, 4)), step_callback=lambda x, *_: x,
                            steering_args={'fk_steering': True})


def test_private_mps_fit_uses_real_fixed_axis_not_padding():
    axis = tuple(AtomIdentity('a', '1037', '', n) for n in ('N', 'CA', 'C', 'CB'))
    target = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    fixed = dict(zip(axis[:3], target[:3]))
    callback = MPSFixedAtomCallback(axis, 8, fixed)
    x = np.zeros((8, 3))
    x[:4] = target @ np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])+5
    x[4:] = 1000
    y = callback(torch.tensor(x[None], device='mps', dtype=torch.float32), 0, 1)
    assert np.allclose(y[0, :4].cpu().numpy(), target, atol=1e-5)
    assert callback.steps[0]['post_overwrite_max_angstrom'] == 0
    assert callback.last_pre_overwrite.shape == (1, 8, 3)
    with pytest.raises(ValueError, match='nonsequential'):
        callback(y, 0, 1)


def test_disabled_callback_and_identity_callback_do_not_change_sampler_output(monkeypatch):
    monkeypatch.setenv('BOLTZ_DENOISER_CACHE', 'off')
    model = Probe('cpu')
    mask = torch.ones((1, 4))
    torch.manual_seed(7)
    original = model.sample(mask)['sample_atom_coords']
    torch.manual_seed(7)
    identity = model.sample(mask, step_callback=lambda x, *_: x)['sample_atom_coords']
    assert torch.equal(original, identity)
