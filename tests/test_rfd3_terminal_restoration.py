"""Native fidelity is a pre-materialization gate; inactive proposals are fixed."""
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, digest, read_records
from dvbfixer.model.diffusion.rfd3_terminal_restoration import (
    audit_native_terminal_transport,
    audit_terminal_derivative,
)
from tests.test_rfd3_refinement import case, move


def coordinates():
    active = AtomIdentity("A", "1", "", "SG")
    inactive = AtomIdentity("a", "1", "A", "SG")
    water = AtomIdentity("W", "1", "", "O")
    return {active: (0., 0., 0.), inactive: (2., 0., 0.), water: (3., 0., 0.)}, frozenset({active})


def test_native_transport_accepts_explicitly_larger_missing_atom_motion_only():
    before, active = coordinates()
    after = {**before, next(iter(active)): (16., 0., 0.)}
    report = audit_native_terminal_transport(before, after, active=active, max_displacement_angstrom=60.)
    assert report.max_active_displacement_angstrom == 16.
    assert report.nonactive_atoms == 2 and report.max_nonactive_residual_angstrom == 0.
    assert not report.scientific_acceptance


@pytest.mark.parametrize("mutation", ["water", "inactive", "missing", "extra", "nan", "bound", "mask"])
def test_native_nonactive_drift_or_inventory_failure_cannot_be_hidden_by_copyback(mutation):
    before, active = coordinates()
    after, bound = dict(before), 60.
    if mutation in {"water", "inactive"}:
        atom = next(a for a in before if (a.chain == "W") == (mutation == "water") and a not in active)
        after[atom] = (4., 0., 0.)
    elif mutation == "missing":
        after.pop(next(iter(before)))
    elif mutation == "extra":
        after[AtomIdentity("X", "1", "", "O")] = (0., 0., 0.)
    elif mutation == "nan":
        after[next(iter(active))] = (float("nan"), 0., 0.)
    elif mutation == "bound":
        bound = True
    else:
        active = frozenset()
    with pytest.raises(CustodyError):
        audit_native_terminal_transport(before, after, active=active, max_displacement_angstrom=bound)


def test_terminal_derivative_cannot_move_a_different_missing_atom():
    request, source, context, raw = case()
    atoms, _, _, _ = read_records(raw.decode())
    active = frozenset([next(a.identity for a in reversed(atoms) if a.identity in request.generated_atoms)])
    # move() edits the first atom of residue2, not the last missing atom.
    with pytest.raises(CustodyError, match="inactive"):
        audit_terminal_derivative(request, source, context, raw, move(raw), active=active,
            raw_sha256=digest(raw), context_sha256=digest(context), max_displacement_angstrom=60.)


def test_terminal_derivative_accepts_unchanged_inactive_missing_atoms():
    request, source, context, raw = case()
    final = move(raw)
    before, _, _, _ = read_records(raw.decode())
    after, _, _, _ = read_records(final.decode())
    active = frozenset(a.identity for a, b in zip(before, after, strict=True) if a.record != b.record)
    report = audit_terminal_derivative(request, source, context, raw, final, active=active,
        raw_sha256=digest(raw), context_sha256=digest(context), max_displacement_angstrom=60.)
    assert report.changed_generated_atoms == 1 and not report.scientific_acceptance
