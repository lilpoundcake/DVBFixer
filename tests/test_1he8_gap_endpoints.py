import numpy as np
import pytest
from scripts.analyze_1he8_joint_gap_endpoints import junction, motion_summary

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def test_motion_is_against_the_given_start_not_a_hidden_reference():
    result = motion_summary([0., 3., 4.])
    assert result["atoms"] == 3
    assert result["rms_angstrom"] == pytest.approx(np.sqrt(25/3))
    assert result["max_angstrom"] == 4.
    with pytest.raises(ValueError, match="finite"):
        motion_summary([float("nan")])


def test_junction_retains_case_and_insertion_code_and_undefined_dihedral():
    left, right = ResidueIdentity("a", "82", "A"), ResidueIdentity("a", "82", "B")
    coordinates = {AtomIdentity(r.chain, r.residue_number, r.insertion_code, name): np.array([float(i), 0., 0.])
                   for i, (r, name) in enumerate(((left, "CA"), (left, "C"), (right, "N"), (right, "CA")))}
    result = junction(coordinates, left, right)
    assert result["left"]["chain"] == "a"
    assert result["right"]["insertion_code"] == "B"
    assert result["C_N_distance_angstrom"] == 1.
    assert result["omega_degrees"] is None
