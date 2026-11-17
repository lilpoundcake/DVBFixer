import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from dvbfixer.model.diffusion import protein_context
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


@pytest.fixture
def helper():
    path = Path(__file__).resolve().parents[1] / "deploy/protenix-v1/checkpoint_gap_smoke.py"
    spec = importlib.util.spec_from_file_location("retained_context_axis_test", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def example(monkeypatch):
    residues = (ResidueIdentity("a", "82", "A"), ResidueIdentity("A", "82", "A"))
    atoms = tuple(AtomIdentity(r.chain, r.residue_number, r.insertion_code, "CA") for r in residues)
    entities = tuple(SimpleNamespace(engine_chain=chain, residues=(residue,), target=SimpleNamespace(sequence="A"))
                     for chain, residue in zip(("a", "A"), residues, strict=True))
    monkeypatch.setattr(protein_context, "sampling_entities", lambda request: entities)
    request = SimpleNamespace(protein_context=SimpleNamespace(fixed_conditioning=(SimpleNamespace(atom=atoms[0]),)),
                              generated_atoms=(atoms[1],))
    array = SimpleNamespace(label_entity_id=np.array(["1", "2"]), chain_id=np.array(["a", "A"]),
        res_id=np.array([1, 1]), atom_name=np.array(["CA", "CA"]),
        res_name=np.array(["ALA", "ALA"]), element=np.array(["C", "C"]))
    return request, array, atoms


def test_native_context_axis_preserves_case_distinct_sources_and_insertion_codes(helper, monkeypatch):
    request, array, atoms = example(monkeypatch)
    axis, names, elements = helper._context_atom_axis(array, request)
    assert axis == atoms
    assert names == ("ALA", "ALA") and elements == ("C", "C")


def test_same_sequence_entities_cannot_be_silently_coalesced_or_relabelled(helper, monkeypatch):
    request, array, _ = example(monkeypatch)
    array.label_entity_id[:] = "1"
    with pytest.raises(ValueError, match="entity/ordinal mapping changed"):
        helper._context_atom_axis(array, request)


def test_native_missing_observed_oxt_is_a_failure_not_silent_atom_loss(helper, monkeypatch):
    request, array, atoms = example(monkeypatch)
    oxt = AtomIdentity(atoms[0].chain, atoms[0].residue_number, atoms[0].insertion_code, "OXT")
    request.protein_context.fixed_conditioning += (SimpleNamespace(atom=oxt),)
    with pytest.raises(ValueError, match="incomplete"):
        helper._context_atom_axis(array, request)


def test_native_wrong_sequence_is_not_fixed_by_identity_projection(helper, monkeypatch):
    request, array, _ = example(monkeypatch)
    array.res_name[1] = "GLY"
    with pytest.raises(ValueError, match="sequence changed"):
        helper._context_atom_axis(array, request)
