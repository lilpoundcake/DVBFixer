"""Research beta barrier units and analytic forces, not production qualification."""
import math

import numpy as np
import pytest
from scripts.research_beta_wall import DRIVE_K, FLOOR_A3, SCALE_NM3, WALL_K, admission, beta_force


def test_wrong_source_is_refused_and_val_is_not_a_beta_stereocentre():
    state=dict(target=dict(residue=['A','2','','THR'],target_cip=dict(CA='S',CB='S'),beta_volume_angstrom3=-2.),negative_ca=[])
    assert admission(state)=='refused_initial_non_R_beta'
    state['target']['residue'][-1]='VAL'
    assert admission(state)=='not_applicable_non_THR'


@pytest.mark.parametrize('kind,z_a3',[('wall',.2),('wall',2.),('drive',.2)])
def test_force_energy_units_and_cartesian_derivatives(kind,z_a3):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    system=openmm.System()
    for _ in range(4):
        system.addParticle(12.)
    system.addForce(beta_force([0,1,2,3],kind))
    integrator=openmm.VerletIntegrator(.001)
    context=openmm.Context(system,integrator,openmm.Platform.getPlatformByName('Reference'))
    xyz=np.array([[.1,0.,0.],[0.,0.,0.],[0.,.1,0.],[0.,0.,z_a3/10.]])
    context.setPositions(xyz*nanometer)
    state=context.getState(energy=True,forces=True)
    q=(z_a3/1000)/SCALE_NM3
    expected=.5*WALL_K*max(0.,FLOOR_A3/1000/SCALE_NM3-q)**2 if kind=='wall' else .5*DRIVE_K*(q+1)**2
    assert state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)==pytest.approx(expected)
    analytic=state.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
    for i in range(4):
        for axis in range(3):
            values=[]
            for sign in (-1,1):
                moved=xyz.copy()
                moved[i,axis]+=sign*1e-6
                context.setPositions(moved*nanometer)
                values.append(context.getState(energy=True).getPotentialEnergy().value_in_unit(kilojoules_per_mole))
            numeric=-(values[1]-values[0])/2e-6
            assert math.isclose(numeric,analytic[i,axis],rel_tol=1e-7,abs_tol=1e-4)
    del context,integrator
