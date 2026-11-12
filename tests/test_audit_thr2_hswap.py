"""Independent transplant ownership must reject edits outside the selected mask."""
import pytest
from scripts.audit_thr2_hswap import independent_mask_check


def test_masks_and_unchanged_heavy_are_checked_independently():
    atoms = [dict(key=['A', '2', '', 'THR', name], element=element)
        for name, element in (('CB', 'C'), ('HB', 'H'), ('HA', 'H'))]
    a, b = [[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [3, 0, 0], [4, 0, 0]]
    assert independent_mask_check(atoms, a, b, [a[0], b[1], a[2]], 'HB') == [1]
    assert independent_mask_check(atoms, a, b, [a[0], a[1], b[2]], 'other-H') == [2]
    with pytest.raises(ValueError, match='ownership'):
        independent_mask_check(atoms, a, b, [[0, 1, 0], b[1], a[2]], 'HB')
    with pytest.raises(ValueError, match='ownership'):
        independent_mask_check(atoms, a, b, [a[0], b[1], b[2]], 'HB')
