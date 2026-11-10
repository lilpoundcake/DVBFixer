"""Saved force comparisons do not infer equality from equal positions."""
import pytest
from scripts.audit_geometry_mediators import compare_systems


def test_physical_parameter_difference_is_reported_even_with_same_atom_map():
    import openmm
    a=openmm.System()
    a.addParticle(12.)
    a.addParticle(12.)
    force=openmm.HarmonicBondForce()
    force.addBond(0,1,.1,100.)
    a.addForce(force)
    b=openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(a))
    assert compare_systems(openmm.XmlSerializer.serialize(a),openmm.XmlSerializer.serialize(b),[])==[
        dict(type='HarmonicBondForce',exact=True)]
    b.getForce(0).setBondParameters(0,0,1,.1,200.)
    assert not compare_systems(openmm.XmlSerializer.serialize(a),openmm.XmlSerializer.serialize(b),[])[0]['exact']


def test_changed_system_inventory_is_not_silently_compared_as_an_intersection():
    import openmm
    a,b=openmm.System(),openmm.System()
    a.addParticle(12.)
    with pytest.raises(ValueError,match='inventory'):
        compare_systems(openmm.XmlSerializer.serialize(a),openmm.XmlSerializer.serialize(b),[])
