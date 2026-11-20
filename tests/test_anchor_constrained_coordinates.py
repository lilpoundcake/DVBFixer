import numpy as np
import pytest
from openmm import (
    Context,
    CustomExternalForce,
    HarmonicBondForce,
    Platform,
    System,
    VerletIntegrator,
)
from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole, nanometer
from scipy.optimize import minimize
from scripts.anchor_constrained_coordinates import ReducedCoordinates


def model():
    topology, system = Topology(), System()
    chain = topology.addChain("d")
    residues = [topology.addResidue(name, chain, "82", insertionCode=icode)
                for name, icode in (("HOH", ""), ("ALA", "A"))]
    for residue, names in zip(residues, (("O", "H1", "H2"), ("C", "H1", "H2")), strict=True):
        for name in names:
            topology.addAtom(name, element.get_by_symbol(name[0]), residue)
    for mass in (0, 1, 1, 12, 1, 1):
        system.addParticle(mass)
    force = HarmonicBondForce()
    axis = list(topology.atoms())
    for parent, h, length in ((0, 1, .096), (0, 2, .096), (3, 4, .109), (3, 5, .109)):
        topology.addBond(axis[parent], axis[h])
        force.addBond(parent, h, length, 400000)
    system.addForce(force)
    points = np.array([[.2, .1, 0], [.3, .1, 0], [.2, .2, 0],
                       [.8, .7, .1], [.8, .7, .21], [.91, .7, .1]])
    return topology, system, points


def test_exact_fixed_coordinates_and_ff_lengths_at_arbitrary_variables():
    topology, system, points = model()
    reference = points.copy()
    chart = ReducedCoordinates(topology, system, points)
    q = chart.initial+np.random.default_rng(7).normal(0, .07, chart.size)
    final = chart.positions(q)
    assert np.array_equal(final[~chart.mobile], reference[~chart.mobile])
    assert np.array_equal(points, reference)
    assert chart.size == 11  # 3 heavy Cartesian + 4*2 hydrogen directions
    assert np.max(np.abs(np.linalg.norm(final[chart.h]-final[chart.parent], axis=1)-chart.length)) < 1.e-14
    assert system.getNumConstraints() == 0


def test_pullback_matches_finite_difference_including_mobile_parent_translation():
    topology, system, points = model()
    chart = ReducedCoordinates(topology, system, points)
    rng = np.random.default_rng(8)
    q = chart.initial+rng.normal(0, .025, chart.size)
    target = points+rng.normal(0, .15, points.shape)
    def energy(x):
        return .5*np.sum((chart.positions(x)-target)**2)
    analytic = chart.pullback(q, chart.positions(q)-target)
    step = 1.e-6
    finite = np.array([(energy(q+np.eye(chart.size)[i]*step)-energy(q-np.eye(chart.size)[i]*step))/(2*step)
                       for i in range(chart.size)])
    np.testing.assert_allclose(analytic, finite, rtol=1.e-7, atol=1.e-9)


def test_projection_matches_independent_dense_constraint_jacobian():
    topology, system, points = model()
    chart = ReducedCoordinates(topology, system, points)
    q = chart.initial+np.random.default_rng(9).normal(0, .025, chart.size)
    positions = chart.positions(q)
    forces = np.random.default_rng(10).normal(size=points.shape)
    free = np.flatnonzero(np.repeat(chart.mobile, 3))
    jacobian = np.zeros((len(chart.h), points.size))
    for row, (h, parent) in enumerate(zip(chart.h, chart.parent, strict=True)):
        n = positions[h]-positions[parent]
        n /= np.linalg.norm(n)
        jacobian[row, 3*h:3*h+3] = n
        jacobian[row, 3*parent:3*parent+3] = -n
    j = jacobian[:, free]
    f = forces.ravel()[free]
    expected = f-j.T@np.linalg.solve(j@j.T, j@f)
    actual = chart.project_forces(positions, forces)
    np.testing.assert_allclose(actual.ravel()[free], expected, atol=1.e-13)
    np.testing.assert_allclose(j@actual.ravel()[free], 0, atol=1.e-13)
    assert np.array_equal(actual[~chart.mobile], np.zeros((1, 3)))


def test_radial_reactions_are_not_misreported_as_nonstationary():
    topology, system, points = model()
    chart = ReducedCoordinates(topology, system, points)
    positions = chart.positions(chart.initial)
    forces = np.zeros_like(points)
    for h, parent in zip(chart.h, chart.parent, strict=True):
        force = 1.e6*(positions[h]-positions[parent])/np.linalg.norm(positions[h]-positions[parent])
        forces[h] += force
        forces[parent] -= force
    report = chart.audit(positions, forces)
    assert report["stationarity_passed"]
    assert report["projected_cartesian_dof_rms_kj_mol_nm"] < 1.e-9
    assert report["xh_equalities_satisfied"] and report["fixed_exact"]


def test_chart_gradient_is_not_the_stationarity_gate_near_antipode():
    topology, system, points = model()
    chart = ReducedCoordinates(topology, system, points)
    q = chart.initial.copy()
    q[chart.heavy_size:] = 1000.
    final = chart.positions(q)
    forces = np.zeros_like(points)
    forces[1] = 1000.*chart.first[0]
    assert np.linalg.norm(chart.pullback(q, -forces)) < 1.e-3
    assert not chart.audit(final, forces)["stationarity_passed"]


def test_openmm_energy_gradient_and_small_constrained_minimization():
    topology, system, points = model()
    chart = ReducedCoordinates(topology, system, points)
    target = points+np.random.default_rng(11).normal(0, .03, points.shape)
    external = CustomExternalForce("500*((x-a)^2+(y-b)^2+(z-c)^2)")
    for name in ("a", "b", "c"):
        external.addPerParticleParameter(name)
    for i in range(len(points)):
        external.addParticle(i, target[i].tolist())
    system.addForce(external)
    context = Context(system, VerletIntegrator(.001), Platform.getPlatformByName("CPU"), {"Threads": "1"})
    def objective(q):
        positions = chart.positions(q)
        context.setPositions(positions*nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        energy = float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
        forces = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
        return energy, chart.pullback(q, -forces)
    before = objective(chart.initial)[0]
    result = minimize(objective, chart.initial, jac=True, method="L-BFGS-B",
                      options={"maxiter": 150, "gtol": 1.e-8, "ftol": 1.e-14})
    assert result.fun < before
    final = chart.positions(result.x)
    context.setPositions(final*nanometer)
    force = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
    report = chart.audit(final, force)
    assert report["fixed_exact"] and report["xh_equalities_satisfied"]
    assert report["projected_cartesian_dof_rms_kj_mol_nm"] < .001
    assert np.array_equal(points[0], final[0])


def test_fixed_h_and_undefined_h_directions_are_refused():
    topology, system, points = model()
    system.setParticleMass(1, 0)
    with pytest.raises(ValueError, match="fixed/observed hydrogens"):
        ReducedCoordinates(topology, system, points)
    system.setParticleMass(1, 1)
    points[1] = points[0]
    with pytest.raises(ValueError, match="undefined initial X"):
        ReducedCoordinates(topology, system, points)
