import numpy as np
import pytest

from dvbfixer.model.diffusion.anchor_symmetry import (
    measure_label_sensitivity,
    residue_label_permutations,
)
from dvbfixer.model.diffusion.contract import ResidueIdentity
from dvbfixer.model.diffusion.protein_atoms import canonical_heavy_atom_identities


def axis(resname):
    names = [a.atom_name for a in canonical_heavy_atom_identities(ResidueIdentity("D", "82", "A"), resname)]
    rng = np.random.default_rng(7)
    return {name: rng.normal(size=3) for name in names}


@pytest.mark.parametrize("resname,first,second", [("ASP", "OD1", "OD2"), ("GLU", "OE1", "OE2"),
                                                ("LEU", "CD1", "CD2"), ("ARG", "NH1", "NH2")])
def test_equivalent_labels_are_an_additional_metric_not_coordinate_repair(resname, first, second):
    reference = axis(resname)
    mobile = {name: point.copy() for name, point in reference.items()}
    mobile[first], mobile[second] = mobile[second], mobile[first]
    before = {name: point.copy() for name, point in mobile.items()}
    result = measure_label_sensitivity(resname, reference, mobile)
    assert result.exact_maximum_angstrom > 0
    assert result.symmetry_minimax_angstrom == pytest.approx(0)
    assert result.automorphism_count >= 2
    for name in mobile:
        np.testing.assert_array_equal(before[name], mobile[name])


def test_backbone_names_cannot_be_permuted_to_hide_real_movement():
    reference = axis("ALA")
    mobile = {name: point.copy() for name, point in reference.items()}
    mobile["CA"] += [0.7, 0, 0]
    result = measure_label_sensitivity("ALA", reference, mobile)
    assert result.symmetry_minimax_angstrom == pytest.approx(result.exact_maximum_angstrom)
    assert result.symmetry_minimax_angstrom == pytest.approx(0.7)


def test_ambiguous_histidine_and_partial_axes_are_explicitly_unavailable():
    with pytest.raises(ValueError, match="unambiguous"):
        residue_label_permutations("HIS", tuple(sorted(axis("HIS"))))
    reference = axis("LEU")
    reference.pop("CD1")
    with pytest.raises(ValueError, match="unambiguous"):
        measure_label_sensitivity("LEU", reference, reference)
