"""Native Reference calibration of analytical LJ/reaction-field/exception forces."""
import math

import pytest
from scripts.audit_ha_nonbonded_pairs import pair_force, parameters


@pytest.mark.parametrize('exception,distance', [(False, .25), (False, 1.2), (True, .25), (True, 1.2)])
def test_exact_native_two_particle_force(exception, distance):
    openmm = pytest.importorskip('openmm')
    from openmm.unit import kilojoules_per_mole, nanometer

    system = openmm.System()
    system.addParticle(1)
    system.addParticle(12)
    force = openmm.NonbondedForce()
    force.setNonbondedMethod(force.CutoffNonPeriodic)
    force.setCutoffDistance(1.)
    force.setReactionFieldDielectric(78.3)
    force.addParticle(.2, .1, .2)
    force.addParticle(-.3, .3, .8)
    if exception:
        force.addException(0, 1, -.017, .18, .13)
    system.addForce(force)
    points = [[0., 0., 0.], [distance, .03, -.01]]
    vector = [a-b for a, b in zip(points[0], points[1])]
    charge, sigma, epsilon = (-.017, .18, .13) if exception else (-.06, .2, .4)
    lj, coulomb = pair_force(vector, charge, sigma, epsilon, 1., 78.3, exception)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        actual = list(context.getState(forces=True).getForces().value_in_unit(kilojoules_per_mole/nanometer)[0])
        assert actual == pytest.approx([a+b for a, b in zip(lj, coulomb)], abs=1e-10, rel=1e-12)
    finally:
        del context
        del integrator


def test_zero_exclusion_is_not_a_short_distance_force():
    assert pair_force([0, 0, 0], 0., .1, 0., 1., 78.3, True) == ([0., 0., 0.], [0., 0., 0.])
    with pytest.raises(ValueError, match='coincident'):
        pair_force([0, 0, 0], .1, .1, .1, 1., 78.3)


@pytest.mark.parametrize('fault', ['nan', 'negative_epsilon', 'zero_cutoff'])
def test_invalid_pair_inputs_are_not_guessed(fault):
    vector = [math.nan, 0., 0.] if fault == 'nan' else [.1, 0., 0.]
    with pytest.raises(ValueError):
        pair_force(vector, .1, .1, -.1 if fault == 'negative_epsilon' else .1,
            0. if fault == 'zero_cutoff' else 1., 78.3)


def test_switching_and_parameter_offsets_are_explicitly_outside_scope():
    openmm = pytest.importorskip('openmm')
    system = openmm.System()
    system.addParticle(1)
    force = openmm.NonbondedForce()
    force.addParticle(0., .1, 0.)
    force.setNonbondedMethod(force.CutoffNonPeriodic)
    force.setUseSwitchingFunction(True)
    force.setSwitchingDistance(.8)
    system.addForce(force)
    with pytest.raises(ValueError, match='unswitched'):
        parameters(system)
    force.setUseSwitchingFunction(False)
    force.addGlobalParameter('lambda', 1.)
    force.addParticleParameterOffset('lambda', 0, .1, 0., 0.)
    with pytest.raises(ValueError, match='no-offset'):
        parameters(system)
