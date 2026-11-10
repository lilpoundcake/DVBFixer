"""Actual Reference contexts verify force expressions, masks and wall gradients."""
import numpy as np
import pytest
from openmm import (
    Context,
    HarmonicBondForce,
    LocalEnergyMinimizer,
    NonbondedForce,
    Platform,
    System,
    VerletIntegrator,
)
from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole, nanometer

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, AnchorPositionRestraint
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.geometry_conditioning import (
    add_anchor_boundary_walls,
    add_geometry_conditioning,
)


def test_geometry_projection_uses_equilibrium_bonds_and_keeps_fixed_particle():
    system = System()
    system.addParticle(0)
    system.addParticle(12)
    bond = HarmonicBondForce()
    bond.addBond(0, 1, 0.15, 100000)
    system.addForce(bond)
    nonbonded = NonbondedForce()
    for _ in range(2):
        nonbonded.addParticle(0, 0.3, 0.1)
    nonbonded.addException(0, 1, 0, 0.3, 0)
    system.addForce(nonbonded)
    add_geometry_conditioning(system)
    integrator = VerletIntegrator(0.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions([[0, 0, 0], [0.03, 0, 0]])
    integrator.setIntegrationForceGroups(1)
    LocalEnergyMinimizer.minimize(context, 0.01, 100)
    xyz = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer)
    np.testing.assert_array_equal(xyz[0], [0, 0, 0])
    assert np.linalg.norm(xyz[1]) == pytest.approx(0.15, abs=1e-6)
    for mask, alpha in ((1 | (1 << 30), 1.0), (1 | (1 << 30), 0.2), (1 | (1 << 31), 0.2)):
        context.setParameter("soft_scale", 1 if mask & (1 << 30) else 0)
        context.setParameter("soft_alpha", alpha)
        integrator.setIntegrationForceGroups(mask)
        LocalEnergyMinimizer.minimize(context, 0.01, 10)
        assert np.isfinite(context.getState(getEnergy=True, groups=mask).getPotentialEnergy().value_in_unit(kilojoule_per_mole))


def test_anchor_walls_match_grouped_backbone_and_case_sensitive_identity():
    topology, system = Topology(), System()
    restraints = []
    for chain_id in ("D", "d"):
        chain = topology.addChain(chain_id)
        for code in ("", "A"):
            residue = topology.addResidue("ALA", chain, id="82", insertionCode=code)
            for name in ("N", "CA", "C", "O", "CB"):
                topology.addAtom(name, element.carbon, residue)
                system.addParticle(12)
                restraints.append(AnchorPositionRestraint(AtomIdentity(chain_id, "82", code, name), (0, 0, 0), 10000))
    add_anchor_boundary_walls(system, topology, tuple(restraints), AnchorMovementLimits())
    integrator = VerletIntegrator(0.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    xyz = np.zeros((20, 3))
    xyz[0, 0] = 0.11  # N alone: backbone RMSD .55 A, heavy below 1.25 A.
    context.setPositions(xyz)
    state = context.getState(getEnergy=True, getForces=True)
    assert state.getPotentialEnergy().value_in_unit(kilojoule_per_mole) > 0
    forces = state.getForces(asNumpy=True)._value
    assert forces[0, 0] < 0
    np.testing.assert_array_equal(forces[5:], 0)
    xyz[:] = 0
    xyz[1, 0] = 0.06  # CA wall despite backbone RMSD only .3 A.
    context.setPositions(xyz)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) > 0
    xyz[:] = 0
    xyz[4, 0] = 0.13  # Sidechain-heavy ceiling is separate from backbone/CA.
    context.setPositions(xyz)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) > 0
    xyz[:] = 0.01
    context.setPositions(xyz)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == 0


def test_softcore_is_finite_repulsive_and_does_not_replace_final_nonbonded():
    system = System()
    original = NonbondedForce()
    for _ in range(2):
        system.addParticle(12)
        original.addParticle(0, 0.3, 0.4)
    system.addForce(original)
    add_geometry_conditioning(system)
    integrator = VerletIntegrator(0.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions([[0, 0, 0], [0.03, 0, 0]])
    context.setParameter("soft_scale", 1)
    state = context.getState(getEnergy=True, getForces=True, groups=1 << 30)
    energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    assert 0 < energy < 1
    assert state.getForces(asNumpy=True)._value[1, 0] > 0
    context.setParameter("soft_scale", 0)
    assert context.getState(getEnergy=True, groups=1 << 30).getPotentialEnergy()._value == 0
    full_energy = context.getState(getEnergy=True, groups=1 << 31).getPotentialEnergy()._value
    assert full_energy > 1e10  # Original hard Lennard-Jones potential is intact.
    context.setParameter("soft_scale", 1)
    assert context.getState(getEnergy=True, groups=1 << 31).getPotentialEnergy()._value == full_energy
