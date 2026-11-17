from types import SimpleNamespace

import numpy as np
from scripts.audit_retained_loop_context_sampling import fixed_records
from scripts.audit_retained_loop_modeller_controls import restore_generated_only

from dvbfixer.model.diffusion.contract import AtomIdentity, GapKind, ResidueIdentity


def test_classical_restoration_does_not_copy_native_observation_drift():
    left, generated, right = (ResidueIdentity("a", str(n), "A") for n in (82, 83, 84))
    atoms = tuple(AtomIdentity(r.chain, r.residue_number, r.insertion_code, "CA") for r in (left, generated, right))
    line = "ATOM  {serial:5d}  CA  ALA a{number:4d}A      1.000   2.000   3.000  1.00 20.00           C\n"
    source = line.format(serial=1, number=82)+line.format(serial=2, number=84)+"END\n"
    request = SimpleNamespace(fixed_atoms=(atoms[0], atoms[2]), generated_atoms=(atoms[1],),
        gaps=(SimpleNamespace(gap_kind=GapKind.INTERNAL, left_anchor=left, right_anchor=right, generated_residues=(generated,)),),
        target_sequences=(SimpleNamespace(chain="a", sequence="AAA"),))
    drifted = np.array([[100., 200., 300.], [4., 5., 6.], [-100., -200., -300.]])
    native = {a: ("ALA", "C", xyz) for a, xyz in zip(atoms, drifted, strict=True)}
    restored = restore_generated_only(source, atoms, drifted, native, request)
    assert fixed_records(restored, set(request.fixed_atoms)) == fixed_records(source, set(request.fixed_atoms))
    assert "   4.000   5.000   6.000" in restored
