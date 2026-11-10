"""Object-index transport moves coordinates, not labels or guessed aliases."""
from unittest.mock import patch

import pytest
from scripts.preflight_graph_methyl_transport import mapped_start


def test_exact_bijection_maps_repaired_coordinates_to_full_indices():
    points = [[1., 2., 3.], [4., 5., 6.]]
    rows = [{'key': ['a', '82', 'A', 'LYS', 'HZ3']}, {'key': ['a', '82', 'A', 'LYS', 'HZ2']}]
    with patch('scripts.preflight_graph_methyl_transport.post_h_repair', return_value=(points, 1)) as repair:
        result, count = mapped_start(rows, [], points, [1, 0])
    assert result == points[::-1] and count == 1
    repair.assert_called_once_with(rows, [], points)


def test_no_partial_lineage_transport():
    with pytest.raises(ValueError, match='bijection'):
        mapped_start([{}, {}], [], [[0., 0., 0.]]*2, [0, 0])
