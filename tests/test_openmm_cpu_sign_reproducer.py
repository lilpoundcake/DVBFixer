"""Synthetic pure-math/System engineering controls, no native minimization."""
import random

import numpy as np
import pytest
from scripts.research_openmm_cpu_sign import (
    analytic,
    energy_terms,
    envelope,
    finite_difference,
    make_system,
    read,
    write,
)


@pytest.mark.parametrize('stress', [1000., 1e14])
@pytest.mark.parametrize('constrained', [False, True])
def test_independent_termwise_derivative_and_translation_free_pair(stress, constrained):
    p = np.array([[-0.05, 0.01, -0.02], [0.050005, 0.012, -0.01], [1., 1., 0.]])
    physical, gp = analytic(p, stress, constrained, 1e6)
    derivative = finite_difference(p, stress, constrained, 1e6)
    assert derivative == pytest.approx(physical+gp, abs=1e-5, rel=1e-6)
    assert physical[0]+physical[1] == pytest.approx(np.zeros(3), abs=1e-12)
    assert gp[0]+gp[1] == pytest.approx(np.zeros(3), abs=1e-12)


def test_wrong_sign_is_not_positive_energy_derivative():
    p = np.array([[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [1., 1., 0.]])
    physical, gp = analytic(p, 1e14, True, 1e6)
    fd = finite_difference(p, 1e14, True, 1e6)
    assert fd[:2] == pytest.approx((physical+gp)[:2], abs=1e-5)
    assert np.max(np.abs(fd[:2]-(physical-gp)[:2])) > 9.9
    assert energy_terms(p, 1e14, True, 1e6)[2] > 0


def test_envelope_does_not_hide_source_sign_with_original_gpu_allowance():
    policy = dict(gradient_absolute=1e-3, gradient_relative=2e-5)
    p = np.array([[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [1., 1., 0.]])
    physical, gp = analytic(p, 1e14, True, 1e6)
    native = (physical-gp).astype(np.float32).astype(np.float64)
    allowance = np.zeros((3, 3))
    allowance[:2] = 32*np.finfo(np.float32).eps*1e6*(0.100005+0.1)
    assert not envelope(physical+gp, native, policy, allowance)['allowance_passed']
    assert envelope(physical-gp, native, policy, allowance)['allowance_passed']


def test_system_has_only_declared_synthetic_forces_constraints_and_no_rng():
    import openmm

    before = random.getstate()
    system = make_system(1e14, True)
    assert system.getNumParticles() == 3 and system.getNumConstraints() == 1
    assert [f.__class__.__name__ for f in system.getForces()] == ['HarmonicBondForce', 'CustomExternalForce']
    assert not any(system.isVirtualSite(i) for i in range(3))
    assert make_system(1000., False).getNumConstraints() == 0
    assert openmm.XmlSerializer.serialize(system) == openmm.XmlSerializer.serialize(make_system(1e14, True))
    assert random.getstate() == before


def test_exclusive_private_json_transport_and_nonfinite_refusal(tmp_path):
    p = tmp_path/'record.json.gz'
    write(p, dict(value=[1., 2., 3.]))
    assert read(p) == dict(value=[1., 2., 3.])
    with pytest.raises(FileExistsError):
        write(p, {})
    with pytest.raises(ValueError):
        write(tmp_path/'bad.json', dict(value=float('nan')))
    assert not (tmp_path/'bad.json').exists()
