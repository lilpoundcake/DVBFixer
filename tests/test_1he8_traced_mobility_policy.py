from types import SimpleNamespace

import pytest
from scripts import run_1he8_traced_mobility as study

from dvbfixer.model.diffusion.contract import ResidueIdentity


def test_sequence_window_is_frozen_nonrecursive_and_excludes_generated_positions(monkeypatch):
    axis = [ResidueIdentity("a", str(i)) for i in range(1, 11)]
    observed = tuple(axis[i] for i in (0, 1, 2, 6, 7, 8, 9))
    request = SimpleNamespace(gaps=[SimpleNamespace(left_anchor=axis[2], right_anchor=axis[6])],
                              sequence_placements=[SimpleNamespace(observed_residues=observed)])
    monkeypatch.setattr(study, "target_axis", lambda _: {"a": axis})
    assert study.anchor_residues(request, "source", "proposal", "fixed") == ()
    assert set(study.anchor_residues(request, "source", "proposal", "direct")) == {axis[2], axis[6]}
    assert set(study.anchor_residues(request, "source", "proposal", "sequence2")) == {axis[i] for i in (0, 1, 2, 6, 7, 8)}
    with pytest.raises(ValueError, match="undeclared"):
        study.anchor_residues(request, "source", "proposal", "unrestrained-all-observations")


def test_spatial_arm_reuses_shared_frozen_shell_policy(monkeypatch):
    selected = (ResidueIdentity("a", "82", "A"), ResidueIdentity("A", "82"))
    request = SimpleNamespace(gaps=[])
    calls = []

    def mobility(req, source, proposal):
        calls.append((req, source, proposal))
        return SimpleNamespace(observed_residues=selected)

    monkeypatch.setattr(study, "spatial_mobility", mobility)
    assert study.anchor_residues(request, "original", "initial-proposal", "spatial5-padding1") == selected
    assert calls == [(request, "original", "initial-proposal")]
