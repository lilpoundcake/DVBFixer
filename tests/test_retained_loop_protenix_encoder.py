from copy import deepcopy
from types import SimpleNamespace

import gemmi
import numpy as np
import pytest
from scripts.materialize_retained_loop_protenix_inputs import audit_template

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def example(*, omit_oxt=False, duplicate_atom=False, moved_atom=False, engine_chain="B"):
    observed = (ResidueIdentity("a", "82", "A"), ResidueIdentity("a", "82", "B"))
    entity = SimpleNamespace(engine_chain="B", target=SimpleNamespace(chain="a", sequence="AAAA"),
        placement=SimpleNamespace(observed_target_indices=(0, 3), observed_residues=observed))
    structure = gemmi.Structure()
    model, chain = gemmi.Model("1"), gemmi.Chain(engine_chain)
    coordinates = {}
    for index, identity in enumerate(observed):
        residue = gemmi.Residue()
        residue.name, residue.het_flag = "ALA", "A"
        residue.seqid = gemmi.SeqId(int(identity.residue_number), identity.insertion_code)
        names = ("N", "CA", "C", "O", "CB") + (("OXT",) if index == 1 else ())
        for ordinal, name in enumerate(names):
            atom = gemmi.Atom()
            atom.name, atom.element = name, gemmi.Element(name[0])
            xyz = np.array([index*10.+ordinal, 2., 3.])
            coordinates[AtomIdentity("a", identity.residue_number, identity.insertion_code, name)] = xyz
            if omit_oxt and name == "OXT":
                continue
            atom.pos = gemmi.Position(*(xyz + ([.1, 0., 0.] if moved_atom and name == "CA" else [0., 0., 0.])))
            residue.add_atom(atom)
            if duplicate_atom and name == "N":
                residue.add_atom(atom.clone())
        chain.add_residue(residue)
    # Same author number in a case-distinct source chain must not enter this entity.
    coordinates[AtomIdentity("A", "82", "A", "CA")] = np.array([100., 100., 100.])
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    payload = [{"mmcif": structure.make_mmcif_document().as_string(),
                "queryIndices": [0, 3], "templateIndices": [0, 1]}]
    return payload, entity, coordinates


def test_discontiguous_template_preserves_chain_case_insertion_codes_and_oxt():
    payload, entity, coordinates = example()
    assert audit_template(payload, entity, coordinates) == 11


@pytest.mark.parametrize("query_indices,template_indices", [([0, 1], [0, 1]), ([3, 0], [0, 1]), ([0, 3], [1, 0])])
def test_generated_target_positions_and_reordered_observations_cannot_enter_masks(query_indices, template_indices):
    payload, entity, coordinates = example()
    payload[0]["queryIndices"], payload[0]["templateIndices"] = query_indices, template_indices
    with pytest.raises(ValueError, match="alignment differs"):
        audit_template(payload, entity, coordinates)


def test_terminal_oxt_loss_is_not_silently_accepted():
    payload, entity, coordinates = example(omit_oxt=True)
    with pytest.raises(ValueError, match="lost observed atoms"):
        audit_template(payload, entity, coordinates)


def test_changed_coordinate_is_not_hidden_by_template_reserialization():
    payload, entity, coordinates = example(moved_atom=True)
    with pytest.raises(ValueError, match="changed an observed coordinate"):
        audit_template(payload, entity, coordinates)


def test_duplicate_atom_or_wrong_engine_entity_fails_closed():
    payload, entity, coordinates = example(duplicate_atom=True)
    with pytest.raises(ValueError, match="duplicate"):
        audit_template(payload, entity, coordinates)
    payload, entity, coordinates = example(engine_chain="b")
    with pytest.raises(ValueError, match="engine entity identity"):
        audit_template(payload, entity, coordinates)


def test_repeated_template_record_and_non_observed_atom_are_rejected():
    payload, entity, coordinates = example()
    with pytest.raises(ValueError, match="ambiguous"):
        audit_template(payload + deepcopy(payload), entity, coordinates)
    del coordinates[AtomIdentity("a", "82", "A", "CB")]
    with pytest.raises(ValueError, match="non-observed"):
        audit_template(payload, entity, coordinates)
