"""Saved native force attribution is exact readback, never optimizer replay."""
import math

import pytest
from scripts.attribute_carbon_H_whole import finite_gradient_error, readback
from scripts.audit_stereo_force_projection import volume_gradient


def test_gradient_is_regular_at_nearplanar_positive_and_negative_volumes():
    for z in (.001, 0., -.001):
        points = [[0., 0., 0.], [.14, 0., 0.], [0., .14, 0.], [.05, .03, z]]
        assert finite_gradient_error(points) < 1e-10
        volume, _, norm = volume_gradient(points)
        assert norm > 0
        assert volume == pytest.approx(.14*.14*z)


def fixture():
    openmm = pytest.importorskip('openmm')
    rows = [dict(key=['d', '82', 'A', 'THR', name], element=e)
        for name, e in [('CA', 'C'), ('N', 'N'), ('C', 'C'), ('CB', 'C')]]
    points = [[0., 0., 0.], [.14, 0., 0.], [0., .14, 0.], [.05, .03, .001]]
    system = openmm.System()
    for _ in rows:
        system.addParticle(12.)
    force = openmm.CustomExternalForce('0.5*k*(x*x+y*y+z*z)')
    force.addGlobalParameter('k', 100.)
    for i in range(4):
        force.addParticle(i, [])
    force.setForceGroup(9)
    system.addForce(force)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 3, .1, 500.)
    system.addForce(bond)
    centre = dict(label='target_ca', residue=['d', '82', 'A', 'THR'], names=['CA', 'N', 'C', 'CB'])
    energy = .5*100*sum(x*x for p in points for x in p)+.5*500*(math.dist(points[0], points[3])-.1)**2
    return openmm, rows, points, system, centre, energy


def test_original_native_force_energy_and_projection_close_without_minimization(monkeypatch):
    openmm, rows, points, system, centre, energy = fixture()
    original = openmm.XmlSerializer.serialize(system)
    def forbidden(*args, **kwargs):
        raise AssertionError('saved-only attribution must not minimize')
    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', forbidden)
    result, forces = readback(system, points, {'k':100.}, energy, [centre], rows)
    assert result['maximum_native_force_closure_error_kj_mol_nm'] < 1e-10
    assert result['energy_closure_error_kj_mol'] < 1e-10
    assert result['centres'][0]['volume_angstrom3'] == pytest.approx(.0196)
    assert result['centres'][0]['finite_difference_max_error_nm2'] < 1e-10
    assert len(forces['component_forces_kj_mol_nm']) == 2
    assert openmm.XmlSerializer.serialize(system) == original


@pytest.mark.parametrize('fault', ['energy', 'parameters', 'identity'])
def test_saved_energy_global_or_exact_identity_drift_refuses(fault):
    _, rows, points, system, centre, energy = fixture()
    parameters = {'k':100.}
    if fault == 'energy':
        energy += 1.
    elif fault == 'parameters':
        parameters['k'] = 99.
    else:
        centre['residue'][0] = 'D'
    with pytest.raises(ValueError):
        readback(system, points, parameters, energy, [centre], rows)
