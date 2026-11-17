from types import SimpleNamespace

import numpy as np
import pytest
from scripts.check_retained_context_reinjection import represented_fixed

from dvbfixer.model.diffusion.contract import AtomIdentity


def example():
    inside = AtomIdentity("a", "82", "A", "CA")
    outside = AtomIdentity("A", "82", "A", "CA")
    coordinates = {inside: np.array([1., 2., 3.]), outside: np.array([4., 5., 6.])}
    request = SimpleNamespace(fixed_atoms=(inside, outside), protein_context=SimpleNamespace(
        fixed_conditioning=(SimpleNamespace(atom=inside),), preserved_only=(outside,)))
    return request, inside, outside, coordinates


def test_context_callback_receives_only_represented_observations():
    request, inside, outside, coordinates = example()
    selected = represented_fixed(request, (inside,), coordinates)
    assert set(selected) == {inside}
    assert selected[inside] is coordinates[inside]
    assert outside in coordinates and outside in request.fixed_atoms


def test_preserved_only_atoms_cannot_become_conditioning_by_accident():
    request, inside, outside, coordinates = example()
    with pytest.raises(ValueError, match="preserved-only"):
        represented_fixed(request, (inside, outside), coordinates)


def test_observation_partition_and_source_coordinates_are_required():
    request, inside, outside, coordinates = example()
    request.protein_context.preserved_only = (inside, outside)
    with pytest.raises(ValueError, match="partition"):
        represented_fixed(request, (inside,), coordinates)
    request.protein_context.preserved_only = (outside,)
    with pytest.raises(ValueError, match="incomplete"):
        represented_fixed(request, (inside,), {inside: coordinates[inside]})
