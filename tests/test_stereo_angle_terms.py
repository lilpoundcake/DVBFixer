"""Independent energy/force derivative and sparse-projection controls."""
import math

import pytest
from scripts.audit_stereo_angle_terms import angle_energy_forces, hydrogen_class, term_projection


def test_angle_derivatives_match_finite_difference_and_translation_invariance():
    points = [[.12, .03, -.04], [.02, -.01, .02], [-.03, .13, .07]]
    reference, stiffness = 1.9, 500.
    _, energy, forces = angle_energy_forces(points, reference, stiffness)
    assert energy >= 0
    for i in range(3):
        for j in range(3):
            first, second = [p[:] for p in points], [p[:] for p in points]
            first[i][j] += 1e-7
            second[i][j] -= 1e-7
            derivative = (angle_energy_forces(first, reference, stiffness)[1]-angle_energy_forces(second, reference, stiffness)[1])/2e-7
            assert forces[i][j] == pytest.approx(-derivative, abs=2e-6)
    assert all(sum(f[j] for f in forces) == pytest.approx(0, abs=1e-10) for j in range(3))
    moved = [[x+4, y-2, z+7] for x, y, z in points]
    _, moved_energy, moved_forces = angle_energy_forces(moved, reference, stiffness)
    assert moved_energy == pytest.approx(energy)
    for first, second in zip(forces, moved_forces, strict=True):
        assert second == pytest.approx(first)


def test_angle_forces_match_openmm_without_minimization():
    openmm = pytest.importorskip('openmm')
    from openmm.unit import kilojoules_per_mole, nanometer

    points = [[.12, .03, -.04], [.02, -.01, .02], [-.03, .13, .07]]
    _, expected_energy, expected_forces = angle_energy_forces(points, 1.9, 500.)
    system = openmm.System()
    for _ in points:
        system.addParticle(12)
    force = openmm.HarmonicAngleForce()
    force.addAngle(0, 1, 2, 1.9, 500.)
    system.addForce(force)
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions(points*nanometer)
        state = context.getState(energy=True, forces=True)
        assert state.getPotentialEnergy().value_in_unit(kilojoules_per_mole) == pytest.approx(expected_energy)
        for actual, expected in zip(state.getForces().value_in_unit(kilojoules_per_mole/nanometer), expected_forces, strict=True):
            assert tuple(actual) == pytest.approx(expected)
    finally:
        del context
        del integrator


@pytest.mark.parametrize('points', [
    [[0, 0, 0]]*3,
    [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
    [[math.nan, 0, 0], [1, 0, 0], [0, 1, 0]],
])
def test_undefined_derivatives_are_not_guessed(points):
    with pytest.raises(ValueError):
        angle_energy_forces(points, 1.9, 500.)


@pytest.mark.parametrize('elements,expected', [
    (['C', 'N', 'C'], 'heavy_only'), (['C', 'C', 'H'], 'H_containing'),
    (['D', 'N', 'C'], 'H_containing'),
])
def test_angle_classification_uses_elements(elements, expected):
    assert hydrogen_class(elements) == expected


def test_unknown_elements_are_refused():
    with pytest.raises(ValueError):
        hydrogen_class(['C', None, 'H'])


def test_sparse_projection_counts_only_the_declared_focal_atoms():
    forces = [[2., 0, 0], [0, 3., 0], [0, 0, 4.]]
    gradient = [[1., 0, 0], [0, 0, 2.]]
    assert term_projection([7, 8, 9], forces, [7, 9], gradient, math.sqrt(5)) == pytest.approx(10/math.sqrt(5))
