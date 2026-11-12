"""Toy algebra only; no enrollment execution, H sampling or minimization."""
import math

import pytest
from scripts.audit_thr2_h_initialization import initializer
from scripts.research_thr2_initializer_frames import (
    determinant,
    initialize_offset,
    inverse,
    rotate,
    rotations,
)


def test_all_24_are_distinct_proper_orthogonal_rotations():
    matrices = rotations()
    assert len(matrices) == len({tuple(tuple(row) for row in m) for m in matrices}) == 24
    for m in matrices:
        assert determinant(m) == 1
        assert all(sum(row[i]*other[i] for i in range(3)) == (a == b)
                   for a, row in enumerate(m) for b, other in enumerate(m))


def test_native_offset_formula_matches_existing_saved_rng_formula_exactly():
    p, neighbors, draws = [0.3, 0.2, 0.1], [[0.1, 0.2, 0.1], [0.3, 0.4, 0.1], [0.3, 0.2, 0.3]], [0.2, 0.4, 0.8]
    _, expected = initializer(p, neighbors, draws)
    assert initialize_offset(p, neighbors, [0.05*d for d in draws]) == expected


def test_covariant_offset_and_translation_roundtrip_for_every_declared_rotation():
    p, neighbors, offset = [0.3, 0.2, 0.1], [[0.1, 0.2, 0.1], [0.3, 0.4, 0.1], [0.3, 0.2, 0.3]], [0.01, 0.02, 0.04]
    expected = initialize_offset(p, neighbors, offset)
    for matrix in rotations():
        for shift in ([0., 0., 0.], [1.25, -0.75, 0.5]):
            def transform(point):
                return [a+b for a, b in zip(rotate(matrix, point), shift, strict=True)]
            point = initialize_offset(transform(p), [transform(n) for n in neighbors], rotate(matrix, offset))
            assert math.dist(inverse(matrix, point, shift), expected) < 1e-12


def test_fixed_lab_offset_is_not_claimed_covariant():
    p, neighbors, offset = [0., 0., 0.], [[-0.1, 0., 0.], [0., -0.1, 0.], [0., 0., -0.1]], [0.01, 0.02, 0.04]
    matrix = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    point = initialize_offset(rotate(matrix, p), [rotate(matrix, n) for n in neighbors], offset)
    assert math.dist(inverse(matrix, point, [0., 0., 0.]), initialize_offset(p, neighbors, offset)) > 0.01


def test_singular_neighbor_sum_refuses_instead_of_inventing_a_direction():
    with pytest.raises(ValueError, match='Singular'):
        initialize_offset([0., 0., 0.], [[1., 0., 0.], [-1., 0., 0.]], [0.01, 0.02, 0.03])
