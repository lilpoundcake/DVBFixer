"""Same-case proposal transfer cannot mutate observations or collapse identities."""
import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_assembly_copy import transfer_missing_assembly_copy
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError


def fixture():
    donor = AtomIdentity("A", "7", "B", "SG")
    recipient = AtomIdentity("a", "7", "B", "SG")
    known = AtomIdentity("a", "16", "", "SG")
    return dict(coordinates={donor: (1., 2., 3.), recipient: (0., 0., 0.), known: (5., 5., 5.)},
                donor_chain="A", recipient_chain="a", donor_atoms=frozenset({donor}), recipient_atoms=frozenset({recipient}),
                observed=frozenset({known}), rotation=((-1., 0., 0.), (0., -1., 0.), (0., 0., 1.)), translation_angstrom=(0., 0., 5.))


def test_transfer_uses_exact_chain_and_insertion_code_only_missing_recipient_changes():
    values = fixture()
    after = transfer_missing_assembly_copy(**values)
    recipient = next(iter(values["recipient_atoms"]))
    assert after[recipient] == (-1., -2., 8.)
    assert all(after[a] == xyz for a, xyz in values["coordinates"].items() if a != recipient)
    assert after.keys() == values["coordinates"].keys()


@pytest.mark.parametrize("problem", ["observed-recipient", "observed-donor", "missing-donor", "icode", "wrong-case", "reflection", "scale", "nan", "wrong-shape"])
def test_unapproved_or_chemically_mismatched_copy_is_refused(problem):
    values = fixture()
    if problem == "observed-recipient":
        values["observed"] |= values["recipient_atoms"]
    elif problem == "observed-donor":
        values["observed"] |= values["donor_atoms"]
    elif problem == "missing-donor":
        values["coordinates"].pop(next(iter(values["donor_atoms"])))
    elif problem == "icode":
        values["donor_atoms"] = frozenset({AtomIdentity("A", "7", "", "SG")})
    elif problem == "wrong-case":
        values["recipient_chain"] = "A"
    elif problem == "reflection":
        values["rotation"] = np.diag([-1., 1., 1.])
    elif problem == "scale":
        values["rotation"] = np.eye(3) * 2.
    elif problem == "nan":
        values["translation_angstrom"] = (0., float("nan"), 0.)
    else:
        values["rotation"] = ((1., 0.), (0., 1.))
    with pytest.raises(CustodyError):
        transfer_missing_assembly_copy(**values)


def test_decimal_rounding_of_real_biomt_is_admitted_without_reflection():
    values = fixture()
    values["rotation"] = ((-.5, .866025, 0.), (.866025, .5, 0.), (0., 0., -1.))
    after = transfer_missing_assembly_copy(**values)
    assert all(np.all(np.isfinite(xyz)) for xyz in after.values())
