"""One missing atom can invert a centre without touching observations or angles."""
import math

import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError
from dvbfixer.model.diffusion.rfd3_missing_stereocentre import reflect_reviewed_missing_atom


def fixture():
    center, first, second, atom, fixed = [AtomIdentity("C", "60", "", name) for name in ("CB", "CA", "OG1", "CG2", "N")]
    return dict(coordinates={center: (0., 0., 0.), first: (1., 0., 0.), second: (0., 1., 0.), atom: (.3, .4, -1.), fixed: (1., 2., 3.)},
                atom=atom, plane=(center, first, second), allowed_missing=frozenset({atom}), observed=frozenset({fixed}),
                max_displacement_angstrom=3.)


def test_reviewed_missing_reflection_keeps_plane_bond_distances_and_every_other_atom():
    values = fixture()
    before = values["coordinates"]
    after = reflect_reviewed_missing_atom(**values)
    assert after[values["atom"]] == (.3, .4, 1.)
    for anchor in values["plane"]:
        assert math.dist(after[values["atom"]], after[anchor]) == pytest.approx(math.dist(before[values["atom"]], before[anchor]))
    assert all(after[a] == xyz for a, xyz in before.items() if a != values["atom"])


@pytest.mark.parametrize("problem", ["observed", "no-permission", "degenerate", "bound", "nonfinite"])
def test_unapproved_or_degenerate_missing_stereocentre_repairs_are_refused(problem):
    values = fixture()
    if problem == "observed":
        values["observed"] |= {values["atom"]}
    elif problem == "no-permission":
        values["allowed_missing"] = frozenset()
    elif problem == "degenerate":
        values["coordinates"][values["plane"][2]] = (2., 0., 0.)
    elif problem == "bound":
        values["max_displacement_angstrom"] = 1.
    else:
        values["coordinates"][values["atom"]] = (.3, .4, float("nan"))
    with pytest.raises(CustodyError):
        reflect_reviewed_missing_atom(**values)
