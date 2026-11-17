import sys
from types import SimpleNamespace

import pytest
from scripts.run_retained_modeller_fixation_controls import selected_generated_initial

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def test_native_initial_selection_must_cover_all_and_only_generated_atoms():
    left, middle, right = (ResidueIdentity("a", str(n), "A") for n in (82, 83, 84))
    generated = tuple(AtomIdentity("a", "83", "A", name) for name in ("N", "CA"))
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="AAA"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=(0, 2), observed_residues=(left, right)),),
        gaps=(SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=1, stop=2), generated_residues=(middle,)),),
        generated_atoms=generated)
    residues = tuple(SimpleNamespace(index=i) for i in (1, 2, 3))
    selection = [SimpleNamespace(residue=residues[1], name=name) for name in ("N", "CA")]
    model = SimpleNamespace(residues=residues, select_loop_atoms=lambda: selection)
    selected, identities = selected_generated_initial(model, request)
    assert selected is selection and identities == generated
    selection.append(SimpleNamespace(residue=residues[0], name="CA"))
    with pytest.raises(ValueError, match="exactly all generated"):
        selected_generated_initial(model, request)


def test_terminal_engine_oxt_is_not_added_to_frozen_mobility_mask(monkeypatch):
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    residue = ResidueIdentity("a", "82", "A")
    generated = tuple(AtomIdentity("a", "82", "A", name) for name in ("N", "CA"))
    request = SimpleNamespace(target_sequences=(SimpleNamespace(chain="a", sequence="A"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=(), observed_residues=()),),
        gaps=(SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=0, stop=1), generated_residues=(residue,)),),
        generated_atoms=generated)
    native_residue = SimpleNamespace(index=1)
    atoms = [SimpleNamespace(residue=native_residue, name=name) for name in ("N", "CA", "OXT")]
    model = SimpleNamespace(residues=(native_residue,), select_loop_atoms=lambda: atoms)
    monkeypatch.setitem(sys.modules, "modeller", SimpleNamespace(Selection=lambda *args: list(args)))
    selection, identities = selected_generated_initial(model, request)
    assert identities == generated
    assert [atom.name for atom in selection] == ["N", "CA"]
    atoms.append(SimpleNamespace(residue=native_residue, name="CB"))
    with pytest.raises(ValueError, match="exactly all generated"):
        selected_generated_initial(model, request)
    atoms.pop()
    atoms.append(atoms[-1])
    with pytest.raises(ValueError, match="exactly all generated"):
        selected_generated_initial(model, request)
