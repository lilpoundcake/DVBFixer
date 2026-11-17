import numpy as np
import pytest
from scripts.preflight_retained_loop_protenix_native import check_template_mask, token_ownership

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_native_token_masks_exclude_generated_and_cover_observed_positions():
    hidden = np.array([False, True, False])
    result = check_template_mask(np.array([[1, 0, 1]]), hidden, 3)
    assert result["all_observed_context_tokens_covered"]
    with pytest.raises(ValueError, match="generated"):
        check_template_mask(np.array([[1, 1, 1]]), hidden, 3)
    with pytest.raises(ValueError, match="lost observed"):
        check_template_mask(np.array([[1, 0, 0]]), hidden, 3)


def test_native_pair_masks_check_both_axes_not_just_one():
    hidden = np.array([False, True, False])
    clean = np.zeros((1, 3, 3), dtype=bool)
    clean[0, 0, 0] = clean[0, 2, 2] = True
    assert check_template_mask(clean, hidden, 3)["generated_tokens_zero"]
    dirty = clean.copy()
    dirty[0, 0, 1] = True
    with pytest.raises(ValueError, match="pair mask"):
        check_template_mask(dirty, hidden, 3)
    dirty = clean.copy()
    dirty[0, 1, 0] = True
    with pytest.raises(ValueError, match="pair mask"):
        check_template_mask(dirty, hidden, 3)


def test_atom_to_token_handles_case_insertion_codes_and_disposable_generated_oxt():
    atoms = [AtomIdentity("a", "82", "A", "CA"), AtomIdentity("a", "82", "B", "CA"),
             AtomIdentity("a", "82", "B", "OXT"), AtomIdentity("A", "82", "B", "CA")]
    hidden = token_ownership(atoms, np.array([0, 1, 1, 2]), {atoms[1]}, {atoms[0], atoms[3]})
    np.testing.assert_array_equal(hidden, [False, True, False])


def test_token_aliases_missing_indices_and_partial_residue_masks_fail_closed():
    first, second = AtomIdentity("A", "1", "", "CA"), AtomIdentity("A", "2", "", "CA")
    with pytest.raises(ValueError, match="mixes"):
        token_ownership([first, second], np.array([0, 0]), {second})
    with pytest.raises(ValueError, match="incomplete"):
        token_ownership([first, second], np.array([0, 2]), {second})
    with pytest.raises(ValueError, match="partial-residue"):
        token_ownership([first, second], np.array([0, 1]), {second}, [AtomIdentity("A", "2", "", "N")])
