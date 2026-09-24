"""Invariants for the private maintained Boltz-2 hook spike."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

REPO_ROOT = Path(__file__).resolve().parent.parent
PATCH = REPO_ROOT / "deploy/boltz-2/per-step-callback.patch"
SMOKE = REPO_ROOT / "deploy/boltz-2/hook-smoke.py"
MAPPING = REPO_ROOT / "deploy/boltz-2/atom_mapping.py"
CHECKPOINT_SMOKE = REPO_ROOT / "deploy/boltz-2/checkpoint_gap_smoke.py"
INPUT_BUILDER = REPO_ROOT / "deploy/boltz-2/build-template-input.py"
REFINEMENT = REPO_ROOT / "deploy/boltz-2/refine_candidate.py"


def _load_mapping_module():
    spec = importlib.util.spec_from_file_location("boltz_atom_mapping", MAPPING)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_checkpoint_smoke_uses_strict_model_and_direct_identity_output() -> None:
    smoke = CHECKPOINT_SMOKE.read_text(encoding="utf-8")

    assert "Boltz2.load_from_checkpoint(" in smoke
    assert "strict=True" in smoke
    assert "model_module.dvbfixer_step_callback = callback" in smoke
    assert "return_predictions=True" in smoke
    assert "BoltzWriter" not in smoke
    assert "validate_runner_result(request, runner_result" in smoke
    assert "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1" in smoke


def test_checkpoint_input_and_refinement_remain_fail_closed() -> None:
    builder = INPUT_BUILDER.read_text(encoding="utf-8")
    refinement = REFINEMENT.read_text(encoding="utf-8")

    assert 'f"SEQRES {row:3d}' in builder
    assert '"template_id": f"{target.chain}1"' in builder
    assert 'if set(request.fixed_atoms) != source_atoms:' in builder
    assert "refine_generated_region(" in refinement
    assert "generated_atoms=request.generated_atoms" in refinement
    assert "validate_runner_result(request, runner_result" in refinement
    assert 'backend="boltz-2-hook-spike+boundary-refinement"' in refinement


def test_atom_mapping_validates_token_and_feature_axes() -> None:
    mapping = _load_mapping_module()
    residues = (
        ResidueIdentity("A", "10", ""),
        ResidueIdentity("A", "10", "A"),
    )
    request = SimpleNamespace(
        target_sequences=(SimpleNamespace(chain="A", sequence="AG"),),
        sequence_placements=(
            SimpleNamespace(
                chain="A",
                target_length=2,
                observed_target_indices=(0,),
                observed_residues=(residues[0],),
            ),
        ),
        gaps=(
            SimpleNamespace(
                chain="A",
                target_interval=SimpleNamespace(start=1, stop=2),
                generated_residues=(residues[1],),
            ),
        ),
        fixed_atoms=(AtomIdentity("A", "10", "", "CA"),),
        generated_atoms=(AtomIdentity("A", "10", "A", "N"),),
    )
    token_dtype = np.dtype(
        [
            ("token_idx", "i4"),
            ("atom_idx", "i4"),
            ("atom_num", "i4"),
            ("res_idx", "i4"),
            ("asym_id", "i4"),
            ("res_name", "U4"),
        ]
    )
    atom_dtype = np.dtype([("name", "U4")])
    tokenized = SimpleNamespace(
        tokens=np.array(
            [(0, 0, 2, 0, 0, "ALA"), (1, 2, 1, 1, 0, "GLY")],
            dtype=token_dtype,
        ),
        structure=SimpleNamespace(
            atoms=np.array([("N",), ("CA",), ("N",)], dtype=atom_dtype)
        ),
    )

    axis, residue_names, token_indices = mapping.model_atom_axis(
        tokenized,
        request,
        chain_name_by_asym_id={0: "A"},
    )

    assert axis == (
        AtomIdentity("A", "10", "", "N"),
        AtomIdentity("A", "10", "", "CA"),
        AtomIdentity("A", "10", "A", "N"),
    )
    assert residue_names == ("ALA", "ALA", "GLY")
    assert token_indices == (0, 0, 1)

    atom_names = np.zeros((4, 4, 64), dtype=np.int64)
    for atom_index, atom_name in enumerate(("N", "CA", "N")):
        for char_index, char in enumerate(atom_name):
            atom_names[atom_index, char_index, ord(char) - 32] = 1
    atom_to_token = np.zeros((4, 2), dtype=np.int64)
    atom_to_token[0:2, 0] = 1
    atom_to_token[2, 1] = 1
    elements = np.zeros((4, 10), dtype=np.int64)
    elements[0, 7] = elements[2, 7] = 1
    elements[1, 6] = 1
    features = {
        "atom_pad_mask": np.array([1, 1, 1, 0]),
        "atom_to_token": atom_to_token,
        "ref_atom_name_chars": atom_names,
        "ref_element": elements,
    }

    assert mapping.validate_feature_axis(axis, token_indices, features) == (7, 6, 7)
