import numpy as np
import pytest
from openmm import Context, Platform, System, VerletIntegrator
from openmm.unit import kilojoule_per_mole, nanometer
from scripts.shape_1non_control import shape_force


@pytest.mark.parametrize('z', [-.02, .01, .15])
def test_signed_elevation_force_against_independent_energy_and_gradient(z):
    system = System()
    for _ in range(4):
        system.addParticle(12.)
    f = shape_force()
    f.addBond(range(4), [])
    system.addForce(f)
    context = Context(system, VerletIntegrator(.001), Platform.getPlatformByName('Reference'))
    p = np.array([[.15,0,0],[0,0,0],[-.04,.14,0],[.08,.08,z]])
    def energy(p):
        a,b,c = p[[0,2,3]]-p[1]
        s = np.dot(np.cross(a,b),c)/(np.linalg.norm(np.cross(a,b))*np.linalg.norm(c))
        return 5000*max(0,np.sin(np.deg2rad(20))-s)**2
    context.setPositions(p)
    state = context.getState(getEnergy=True,getForces=True)
    assert state.getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(energy(p))
    force = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
    for i in range(12):
        delta = np.eye(12)[i].reshape(4,3)*1.e-7
        assert force.flat[i] == pytest.approx((energy(p-delta)-energy(p+delta))/2.e-7,abs=1.e-4)
    assert (force[3,2] > 0) == (z < .15)
