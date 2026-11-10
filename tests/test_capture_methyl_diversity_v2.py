"""ForceField observer keeps the signature consulted by native ArgTracker."""
import inspect
from unittest.mock import patch

import pytest
from scripts.capture_methyl_diversity_v2 import forcefield_observer


def test_signature_and_defaults_survive_the_actual_native_argtracker():
    pytest.importorskip('openmm')
    from openmm.app import ForceField, Topology

    forcefield = ForceField('amber19/protein.ff19SB.xml', 'amber19/tip3p.xml')
    topology = Topology()
    original = ForceField.createSystem
    systems = []
    observed = forcefield_observer(original, lambda: True, systems)
    assert inspect.signature(observed) == inspect.signature(original)
    with patch.object(ForceField, 'createSystem', observed):
        system = forcefield.createSystem(topology, flexibleConstraints=False)
    assert system.getNumParticles() == 0 and len(systems) == 1
    assert systems[0]['pointer'] == int(system.this)
    assert ForceField.createSystem is original


def test_unwrapped_observer_reproduces_the_retained_V1_instrumentation_blocker():
    pytest.importorskip('openmm')
    from openmm.app import ForceField, Topology

    original = ForceField.createSystem
    forcefield = ForceField('amber19/protein.ff19SB.xml')
    def unwrapped(self, topology, *args, **kwargs):
        return original(self, topology, *args, **kwargs)
    with patch.object(ForceField, 'createSystem', unwrapped):
        with pytest.raises(ValueError, match='flexibleConstraints'):
            forcefield.createSystem(Topology(), flexibleConstraints=False)
    assert ForceField.createSystem is original
