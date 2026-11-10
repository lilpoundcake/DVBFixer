"""Native object lineage permits label mutation but never alias-based remapping."""
import pytest
from scripts.capture_methyl_diversity import object_lineage


class Topology:
    def __init__(self, atoms):
        self._atoms = atoms

    def atoms(self):
        return iter(self._atoms)


def test_native_objects_not_names_establish_the_cross_stage_mapping():
    a, b, c = object(), object(), object()
    assert object_lineage([a, b, c], Topology([c, a, b])) == [1, 2, 0]


@pytest.mark.parametrize('fault', ['replacement', 'missing', 'duplicate', 'extra'])
def test_unavailable_object_lineage_refuses_instead_of_guessing_aliases(fault):
    a, b = object(), object()
    final = dict(replacement=[a, object()], missing=[a], duplicate=[a, a], extra=[a, b, object()])[fault]
    with pytest.raises(ValueError, match='lineage'):
        object_lineage([a, b], Topology(final))
