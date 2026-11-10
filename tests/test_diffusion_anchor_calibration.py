from dataclasses import replace

import numpy as np
import pytest

from dvbfixer.model.diffusion.anchor_calibration import compare_anchor
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def coordinates():
    rng = np.random.default_rng(7)
    atoms = tuple(AtomIdentity("d", "82", "A", name) for name in ("N", "CA", "C", "O", "CB"))
    fit = tuple(AtomIdentity("D", str(i), "", "CA") for i in range(1, 9))
    ref = {a: tuple(rng.normal(size=3)) for a in (*atoms, *fit)}
    return ref, atoms, fit


def test_rigid_frame_motion_is_not_counted_as_local_anchor_displacement():
    ref, _, fit = coordinates()
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    mobile = {a: tuple(np.array(v) @ rotation + [10, -8, 2]) for a, v in ref.items()}
    result = compare_anchor(ref, mobile, ResidueIdentity("d", "82", "A"), fit, contact_atoms=fit)
    assert result.maximum_heavy_displacement_angstrom < 1e-10
    assert result.contact_pairs_lost == result.contact_pairs_gained == 0
    assert "current" in result.limits_passed


def test_anchor_motion_is_not_fitted_away_and_proposed_bounds_are_only_measurements():
    ref, atoms, fit = coordinates()
    mobile = {a: tuple(np.array(v) + ([0.7, 0, 0] if a in atoms else [0, 0, 0])) for a, v in ref.items()}
    result = compare_anchor(ref, mobile, ResidueIdentity("d", "82", "A"), fit)
    assert result.backbone_rmsd_angstrom == pytest.approx(0.7)
    assert result.ca_displacement_angstrom == pytest.approx(0.7)
    assert "current" not in result.limits_passed and "candidate-075" in result.limits_passed
    with pytest.raises(ValueError, match="outside the anchor"):
        compare_anchor(ref, mobile, result.residue, (*fit, atoms[0]))


def test_incomplete_axis_nonfinite_or_degenerate_fits_are_not_calibration_evidence():
    ref, atoms, fit = coordinates()
    mobile = dict(ref)
    mobile.pop(atoms[-1])
    with pytest.raises(ValueError, match="identical complete"):
        compare_anchor(ref, mobile, ResidueIdentity("d", "82", "A"), fit)
    mobile = dict(ref)
    mobile[fit[0]] = (np.nan, 0, 0)
    with pytest.raises(ValueError, match="nonfinite"):
        compare_anchor(ref, mobile, ResidueIdentity("d", "82", "A"), fit)
    ref = {a: (float(i), 0.0, 0.0) for i, a in enumerate(ref)}
    with pytest.raises(ValueError, match="degenerate"):
        compare_anchor(ref, ref, replace(ResidueIdentity("d", "82"), insertion_code="A"), fit)
