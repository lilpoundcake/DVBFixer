import json
from dataclasses import asdict

import pytest
from scripts import run_1he8_buffered_walls as experiment

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits


def test_guard_band_tightens_wall_activation_without_relaxing_validator():
    original = AnchorMovementLimits()
    inner = experiment.buffered_limits(original)
    assert inner.backbone_rmsd_angstrom_max == .45
    assert inner.ca_displacement_angstrom_max == .45
    assert inner.heavy_displacement_angstrom_max == 1.20
    assert asdict(original) == asdict(AnchorMovementLimits())
    with pytest.raises(ValueError, match="positive validation"):
        experiment.buffered_limits(AnchorMovementLimits(.01, .02, .03))


def test_wrapper_only_changes_inner_wall_argument_and_restores_refiner(tmp_path, monkeypatch):
    original_limits = AnchorMovementLimits()
    actual = []

    def refiner(*args, **kwargs):
        actual.append((args, kwargs))
        return "same underlying refiner"

    def worker(index, arm, output):
        assert experiment.boundary.refine_generated_region("unchanged proposal", anchor_movement_limits=original_limits,
            generated_atoms=("all generated",)) == "same underlying refiner"
        assert original_limits.backbone_rmsd_angstrom_max == .5

    monkeypatch.setattr(experiment.boundary, "refine_generated_region", refiner)
    monkeypatch.setattr(experiment.traced, "worker", worker)
    experiment.worker(1, "direct", tmp_path)
    assert experiment.boundary.refine_generated_region is refiner
    assert actual[0][0] == ("unchanged proposal",)
    assert actual[0][1]["generated_atoms"] == ("all generated",)
    assert actual[0][1]["anchor_movement_limits"].backbone_rmsd_angstrom_max == .45
    policy = json.loads((tmp_path / "wall-activation-policy.json").read_text())
    assert policy["validation_limits_unchanged"]["backbone_rmsd_angstrom_max"] == .5
