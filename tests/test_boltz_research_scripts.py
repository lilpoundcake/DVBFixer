"""Static invariants for the private maintained Boltz-2 hook spike."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PATCH = REPO_ROOT / "deploy/boltz-2/per-step-callback.patch"
SMOKE = REPO_ROOT / "deploy/boltz-2/hook-smoke.py"


def test_callback_runs_after_predictor_update() -> None:
    patch = PATCH.read_text(encoding="utf-8")
    update = patch.index("atom_coords = atom_coords_next")
    callback = patch.index("updated = step_callback")

    assert callback > update
    assert 'getattr(self, "dvbfixer_step_callback", None)' in patch
    assert "step_callback must preserve tensor shape" in patch
    assert "step_callback must preserve tensor device" in patch
    assert "step_callback must preserve tensor dtype" in patch


def test_smoke_uses_pinned_real_sampler_and_exact_projection() -> None:
    smoke = SMOKE.read_text(encoding="utf-8")

    assert "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc" in smoke
    assert "from boltz.model.modules.diffusionv2 import AtomDiffusion" in smoke
    assert "num_sampling_steps=2" in smoke
    assert 'step_callback=callback' in smoke
    assert '"maximum_post_projection_error"' in smoke
