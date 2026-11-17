import numpy as np
import pytest
from openmm import HarmonicBondForce, System
from scripts.continue_retained_1he8_segment_frames import (
    particle_graph,
    remaining_arms,
    saved_coordinates,
)

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_exact_complement_does_not_repeat_the_two_previously_checked_arms():
    rows = [{"final_path": f"/retained/{i}/final.pdb", "accepted": i >= 10} for i in range(12)]
    selected = remaining_arms(rows, ["/retained/10", "/retained/11"])
    assert selected == rows[:10]
    assert rows[10]["accepted"]


@pytest.mark.parametrize("change", ["duplicate", "unknown", "omitted", "unaccepted", "invalid-flag"])
def test_selector_refuses_ambiguous_history(change):
    rows = [{"final_path": "/retained/a/final.pdb", "accepted": False},
            {"final_path": "/retained/b/final.pdb", "accepted": True}]
    completed = ["/retained/b"]
    if change == "duplicate":
        rows.append(dict(rows[0]))
    elif change == "unknown":
        completed.append("/retained/missing")
    elif change == "omitted":
        completed.clear()
    elif change == "unaccepted":
        completed.append("/retained/a")
    else:
        rows[0]["accepted"] = 0
    with pytest.raises(ValueError):
        remaining_arms(rows, completed)


def graph():
    system = System()
    for _ in range(3):
        system.addParticle(12.)
    force = HarmonicBondForce()
    force.addBond(0, 1, 0.14, 100.)
    system.addForce(force)
    particles = [{"index": index, "chain": chain, "residue_number": "82", "insertion_code": "A", "atom_name": name}
                 for index, (chain, name) in enumerate((("A", "N"), ("A", "CA"), ("a", "C")))]
    return particles, system


def test_graph_uses_only_serialized_bonds_not_adjacent_particle_indices():
    particles, system = graph()
    atoms, bonds, edges = particle_graph(particles, system)
    pair = tuple(sorted((atoms[0], atoms[1])))
    assert bonds == {pair}
    assert edges == (pair,)
    assert all(atoms[2] not in edge for edge in edges)


@pytest.mark.parametrize("change", ["indices", "duplicate-identity", "count", "no-bonds"])
def test_invalid_particle_graph_refused(change):
    particles, system = graph()
    if change == "indices":
        particles[0]["index"] = 1
    elif change == "duplicate-identity":
        particles[1] = {**particles[0], "index": 1}
    elif change == "count":
        system.addParticle(12.)
    else:
        system.removeForce(0)
    with pytest.raises(ValueError):
        particle_graph(particles, system)


def test_saved_frame_units_are_converted_without_moving_atoms(tmp_path):
    path = tmp_path / "frame.npz"
    atoms = [AtomIdentity("A", "82", "A", "CA")]
    np.savez(path, positions_nm=np.array([[0.1, 0.2, 0.3]]))
    assert np.array_equal(saved_coordinates(path, atoms)[atoms[0]], [1., 2., 3.])


@pytest.mark.parametrize("xyz", [np.zeros((2, 3)), np.array([[np.nan, 0., 0.]])])
def test_invalid_saved_frame_fails_closed(tmp_path, xyz):
    path = tmp_path / "frame.npz"
    np.savez(path, positions_nm=xyz)
    with pytest.raises(ValueError):
        saved_coordinates(path, [AtomIdentity("A", "1", "", "N")])
