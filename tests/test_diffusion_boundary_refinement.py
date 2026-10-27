"""Force-field mobility must distinguish source atoms from temporary completion."""

import numpy as np
from openmm import System
from openmm.app import Topology, element
from openmm.unit import dalton

from dvbfixer.model.diffusion.boundary_refinement import _freeze_observed_heavy_atoms
from dvbfixer.model.diffusion.contract import AtomIdentity


def test_only_observed_non_generated_heavy_atoms_are_frozen():
    topology = Topology()
    system = System()
    original = {}
    generated = set()
    expected_fixed = set()
    # Case-sensitive chains and insertion codes must not cross-freeze an atom.
    for chain_id in ("D", "d"):
        chain = topology.addChain(chain_id)
        for code in ("", "A"):
            residue = topology.addResidue("ALA", chain, id="82", insertionCode=code)
            for name, chemical_element in (("CA", element.carbon), ("CB", element.carbon), ("H", element.hydrogen)):
                atom = topology.addAtom(name, chemical_element, residue)
                system.addParticle(chemical_element.mass)
                key = AtomIdentity(chain_id, "82", code, name)
                if name == "CA":
                    original[key] = np.zeros(3)
                    if chain_id == "D" and not code:
                        generated.add(key)
                    else:
                        expected_fixed.add(atom.index)
                # CB and H exist only in the temporary FF topology.
    _freeze_observed_heavy_atoms(system, topology, original, generated)
    fixed = {i for i in range(system.getNumParticles())
        if system.getParticleMass(i).value_in_unit(dalton) == 0}
    assert fixed == expected_fixed
