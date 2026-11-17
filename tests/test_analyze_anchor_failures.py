import math

import numpy as np
from scripts.analyze_anchor_failures import force_summary, pair_terms, topology_only_contacts


def metadata(name, number, residue, index, element="O"):
    return {"identity": {"chain": "B", "residue_number": str(number), "insertion_code": "", "atom_name": name},
            "resname": residue, "residue_index": index, "element": element, "temporary_addition": False}


def test_short_water_contact_is_not_a_covalent_exclusion_in_control():
    axis = [metadata("OG", 161, "SER", 0), metadata("O", 12, "HOH", 1)]
    coordinates = np.array([[0, 0, 0], [.1658, 0, 0]])
    rows = topology_only_contacts(axis, [], {0, 1}, coordinates, set(), set())
    assert len(rows) == 1
    assert rows[0]["hidden_by_distance_inferred_exclusions"]
    assert rows[0]["both_fixed"] and not rows[0]["touches_generated"]
    assert abs(rows[0]["overlap_angstrom"]-1.382) < 1.e-10
    assert topology_only_contacts(axis, [(0, 1)], {0, 1}, coordinates, set(), set()) == []


def test_actual_1_3_and_1_4_exclusions_remain_in_force():
    axis = [metadata("C", i, "ALA", i, "C") for i in range(4)]
    coordinates = np.array([[0, 0, 0], [.1, 0, 0], [.2, 0, 0], [.05, 0, 0]])
    assert topology_only_contacts(axis, [(0, 1), (1, 2), (2, 3)], set(range(4)), coordinates, set(), set()) == []


def test_zero_epsilon_polar_h_has_no_lj_barrier_to_nonexception_acceptor():
    axis = [metadata("HG", 161, "SER", 0, "H"), metadata("O", 12, "HOH", 1)]
    xyz = np.array([[0, 0, 0], [.0002877, 0, 0]])
    charges, sigmas, epsilons = np.array([.4275, -.834]), np.array([1., .315]), np.array([0., .636])
    pair = pair_terms(0, xyz, charges, sigmas, epsilons, {}, axis)[0]
    assert pair["lj_energy_kj_mol"] == 0
    assert pair["coulomb_energy_kj_mol"] < -170000
    assert pair["pair_force_norm_kj_mol_nm"] > 590000000
    assert not pair["nonbonded_exception"]
    assert pair_terms(0, xyz, charges, sigmas, epsilons, {(0, 1): (0., 0., 0.)}, axis) == []


def test_fixed_reaction_forces_do_not_fake_mobile_nonconvergence():
    forces = np.array([[3., 0, 0], [1000., 0, 0], [6., 0, 0]])
    row = force_summary(forces, np.array([1., 0, 1.]), {0})
    assert row["mobile_gradient_rms_over_all_components"] == math.sqrt(5.)
    assert row["max_mobile_force"] == 6.
    assert row["max_generated_heavy_force"] == 3.
