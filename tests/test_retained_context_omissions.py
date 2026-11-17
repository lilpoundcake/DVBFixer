from types import SimpleNamespace

import pytest
from scripts.audit_retained_context_omissions import context_facts

from dvbfixer.model.diffusion.contract import ResidueIdentity


def atom(serial, chain, number, name="CA", residue="ALA", x=0., record="ATOM  ", icode=" ", element="C"):
    return (f"{record}{serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}{icode}   "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}\n")


def request():
    return SimpleNamespace(gaps=(SimpleNamespace(chain="a", left_anchor=ResidueIdentity("a", "82", "A"),
                                                right_anchor=None),))


def test_exact_omission_instances_and_case_sensitive_anchors_not_distance_admission():
    protein = atom(1, "a", 82, icode="A")+atom(2, "A", 82, icode="A", x=100.)
    context = atom(3, "a", 90, residue="ZN", name="ZN", record="HETATM", x=2., element="Zn")
    context += atom(4, "A", 90, residue="HOH", name="O", record="HETATM", x=99., element="O")
    facts = context_facts(protein+context, protein, request())
    assert facts["omitted_atom_counts_by_residue_name"] == {"ZN": 1, "HOH": 1}
    assert len(facts["omitted_residue_instances"]) == 2
    nearest = facts["gap_anchor_context"][0]["nearest"]["outside_HOH_WAT"]
    assert nearest["distance_angstrom"] == 2. and nearest["observed_anchor_atom"]["chain"] == "a"
    assert not facts["scientific_admission"] and not facts["missing_region_positions_or_contacts_inferred"]


def test_projection_cannot_hide_a_protein_chain_or_coordinate_change():
    protein = atom(1, "a", 82, icode="A")
    with pytest.raises(ValueError, match="omitted an ATOM"):
        context_facts(protein+atom(2, "A", 82), protein, request())
    with pytest.raises(ValueError, match="changed observations"):
        context_facts(protein, atom(1, "a", 82, icode="A", x=1.), request())


def test_empty_omission_does_not_imply_scientific_admission():
    protein = atom(1, "a", 82, icode="A")
    facts = context_facts(protein, protein, request())
    assert facts["omitted_atoms"] == 0 and facts["gap_anchor_context"][0]["nearest"] == {}
    assert not facts["scientific_admission"]
