"""Actual Reference energy minima and observed-only/native-bond safety gates."""
import pytest
from openmm import Context, Platform, System, VerletIntegrator
from openmm.app import Topology, element
from openmm.unit import kilojoule_per_mole

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError
from dvbfixer.model.diffusion.rfd3_disulfide_restraint import (
    build_generated_disulfide_torsion_force,
)


def fixture():
    top = Topology()
    chain = top.addChain("A")
    first, last = top.addResidue("CYS", chain, "7"), top.addResidue("CYS", chain, "16", insertionCode="A")
    a, b = top.addAtom("CB", element.carbon, first), top.addAtom("SG", element.sulfur, first)
    c, d = top.addAtom("SG", element.sulfur, last), top.addAtom("CB", element.carbon, last)
    for pair in ((a, b), (b, c), (c, d)):
        top.addBond(*pair)
    first_id, last_id = AtomIdentity("A", "7", "", "SG"), AtomIdentity("A", "16", "A", "SG")
    return top, (first_id, last_id)


def test_actual_twofold_torsion_minima_are_plus_minus90_not0_or180():
    top, pair = fixture()
    force = build_generated_disulfide_torsion_force(top, pairs=(pair,), mobile=frozenset({pair[0]}), strength_kj_mol=100.)
    assert force.getForceGroup() == 0 and force.getNumTorsions() == 1
    system = System()
    for _ in range(4):
        system.addParticle(12)
    system.addForce(force)
    integrator = VerletIntegrator(.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    for sign in (-1, 1):
        context.setPositions([(0., 1., 0.), (0., 0., 0.), (1., 0., 0.), (1., 0., sign)])
        assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(0.)
    context.setPositions([(0., 1., 0.), (0., 0., 0.), (1., 0., 0.), (1., 1., 0.)])
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(200.)


@pytest.mark.parametrize("problem", ["observed-only", "two-mobile", "duplicate", "missing-icode", "no-bond", "strength"])
def test_undeclared_or_incomplete_disulfide_force_domain_is_refused(problem):
    top, pair = fixture()
    pairs, mobile, strength = (pair,), frozenset({pair[0]}), 100.
    if problem == "observed-only":
        mobile = frozenset({AtomIdentity("B", "7", "", "SG")})
    elif problem == "two-mobile":
        mobile = frozenset(pair)
    elif problem == "duplicate":
        pairs += (tuple(reversed(pair)),)
    elif problem == "missing-icode":
        pairs = ((pair[0], AtomIdentity("A", "16", "", "SG")),)
    elif problem == "no-bond":
        # A fresh topology has matching atoms but no native bond graph.
        fresh = Topology()
        chain = fresh.addChain("A")
        for number, code in (("7", ""), ("16", "A")):
            residue = fresh.addResidue("CYS", chain, number, insertionCode=code)
            fresh.addAtom("CB", element.carbon, residue)
            fresh.addAtom("SG", element.sulfur, residue)
        top = fresh
    else:
        strength = True
    with pytest.raises(CustodyError):
        build_generated_disulfide_torsion_force(top, pairs=pairs, mobile=mobile, strength_kj_mol=strength)
