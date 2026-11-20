import numpy as np
import pytest
from openmm import HarmonicBondForce, System
from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole, nanometer
from scripts.anchor_constrained_control import Objective
from scripts.anchor_constrained_coordinates import ReducedCoordinates
from scripts.audit_anchor_constrained_control import independent_projection


def chart_model():
    topology, system = Topology(), System()
    residue = topology.addResidue("SER", topology.addChain("B"), "161")
    oxygen = topology.addAtom("OG", element.oxygen, residue)
    hydrogen = topology.addAtom("HG", element.hydrogen, residue)
    topology.addBond(oxygen, hydrogen)
    system.addParticle(0)
    system.addParticle(1)
    force = HarmonicBondForce()
    force.addBond(0, 1, .096, 400000)
    system.addForce(force)
    chart = ReducedCoordinates(topology, system, np.array([[0., 0, 0], [.097, 0, 0]]))
    return chart, system


class State:
    def __init__(self, force, energy):
        self.force, self.energy = force, energy

    def getPotentialEnergy(self):
        return self.energy*kilojoule_per_mole

    def getForces(self, asNumpy):
        return self.force*(kilojoule_per_mole/nanometer)


class FakeContext:
    def __init__(self):
        self.force = np.array([[1.e8, 1.e8, 1.e8], [1., 2., 3.]])
        self.energy = 1.
        self.calls = 0

    def setPositions(self, points):
        self.positions = points.value_in_unit(nanometer)

    def getState(self, getEnergy, getForces, groups):
        self.calls += 1
        self.groups = groups
        return State(self.force, self.energy)


def test_objective_retains_force_group_and_matches_chart_pullback():
    chart, _ = chart_model()
    context = FakeContext()
    objective = Objective(chart, context, (1 << 31) | 1)
    energy, gradient = objective.evaluate(chart.initial)
    assert energy == 1
    np.testing.assert_allclose(gradient, chart.pullback(chart.initial, -context.force))
    assert context.groups == (1 << 31) | 1
    assert np.array_equal(context.positions[0], chart.reference[0])
    assert np.linalg.norm(context.positions[1]-context.positions[0]) == pytest.approx(.096)
    objective.evaluate(chart.initial.copy())
    assert objective.calls == context.calls == 1


@pytest.mark.parametrize("kind", ["force", "energy"])
def test_failed_ff_evaluation_is_retained_without_clipping(kind):
    chart, _ = chart_model()
    context = FakeContext()
    if kind == "force":
        context.force[1, 0] = np.inf
    else:
        context.energy = np.nan
    objective = Objective(chart, context, 1)
    with pytest.raises(FloatingPointError, match="nonfinite"):
        objective.evaluate(chart.initial)
    assert objective.failed_state is not None
    assert np.array_equal(objective.failed_state[0], chart.initial)
    assert objective.last is None


def test_independent_sparse_auditor_excludes_radial_and_fixed_reactions():
    chart, system = chart_model()
    points = chart.positions(chart.initial)
    forces = np.array([[1.e12, 1.e12, 1.e12], [1.e8, 3., 4.]])
    rms, tangent_error = independent_projection(system, chart.bonds, points, forces)
    assert rms == pytest.approx(5/np.sqrt(2))
    assert tangent_error == 0
    assert rms == pytest.approx(chart.audit(points, forces)["projected_cartesian_dof_rms_kj_mol_nm"])
