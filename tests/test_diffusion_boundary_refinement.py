"""Force-field mobility must distinguish source atoms from temporary completion."""

import numpy as np
import pytest
from openmm import System
from openmm.app import Topology, element
from openmm.unit import dalton

from dvbfixer.model.diffusion.boundary_refinement import _freeze_observed_heavy_atoms
from dvbfixer.model.diffusion.contract import AtomIdentity


@pytest.mark.parametrize("precondition", [False, True])
def test_full_score_drives_restarts_and_selected_coordinates(precondition):
    from pathlib import Path

    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    fixture = Path(__file__).parent / "fixtures/8cz8/8cz8_a_u.pdb"
    lines = [line for line in fixture.read_text().splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[21] == "C" and 64 <= int(line[22:26]) <= 70
        and line[76:78].strip() != "H"]
    text = "".join(lines) + "END\n"
    mobile = tuple(AtomIdentity("C", line[22:26].strip(), line[26].strip(), line[12:16].strip())
        for line in lines if 65 <= int(line[22:26]) <= 69)
    trials = []

    def score(coordinates):
        trials.append(coordinates)
        # Stand-in for a final gate omitted by the old local score (Rama/chi).
        return (0, 0, 1 if len(trials) == 1 else 0)

    result = refine_generated_region(text,
        generated_residues=tuple(ResidueIdentity("C", str(i)) for i in range(65, 70)),
        generated_atoms=mobile, platform_name="Reference", max_iterations=10,
        restart_count=3, candidate_score=score, geometry_preconditioning=precondition)
    assert len(trials) == 2
    assert set(result.coordinates_angstrom) == set(mobile)
    for atom in mobile:
        np.testing.assert_array_equal(result.coordinates_angstrom[atom], trials[1][atom])
    phases = [d["phase"] for d in result.search_diagnostics]
    assert phases.count("full-validator-ranking") == 2
    assert phases.count("bonded-geometry") == (2 if precondition else 0)


def test_geometry_first_repairs_compressed_proline_ring_with_fixed_flanks():
    from pathlib import Path

    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    fixture = Path(__file__).parent / "fixtures/8cz8/8cz8_a_u.pdb"
    lines = [line for line in fixture.read_text().splitlines(keepends=True)
        if line.startswith("ATOM  ") and line[21] == "C" and 6 <= int(line[22:26]) <= 8
        and line[76:78].strip() != "H"]
    nitrogen = next(line for line in lines if int(line[22:26]) == 7 and line[12:16].strip() == "N")
    xyz_n = np.array([float(nitrogen[i:i + 8]) for i in (30, 38, 46)])
    compressed = xyz_n + [0.1, 0, 0]
    lines = [line[:30] + "".join(f"{v:8.3f}" for v in compressed) + line[54:]
        if int(line[22:26]) == 7 and line[12:16].strip() == "CD" else line for line in lines]
    mobile = tuple(AtomIdentity("C", "7", "", line[12:16].strip())
        for line in lines if int(line[22:26]) == 7)
    result = refine_generated_region("".join(lines) + "END\n",
        generated_residues=(ResidueIdentity("C", "7"),), generated_atoms=mobile,
        platform_name="Reference", max_iterations=100, restart_count=0,
        geometry_preconditioning=True)
    coordinates = result.coordinates_angstrom
    ring_closure = np.linalg.norm(coordinates[AtomIdentity("C", "7", "", "CD")]
        - coordinates[AtomIdentity("C", "7", "", "N")])
    assert 1.3 < ring_closure < 1.65
    assert all(d.get("finite_all_l", True) for d in result.search_diagnostics)


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
