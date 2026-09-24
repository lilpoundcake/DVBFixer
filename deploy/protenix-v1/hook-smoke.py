"""CPU smoke for the maintained Protenix v1 per-step callback patch."""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
from reinjection import FixedAtomReinjector

from dvbfixer.model.diffusion.contract import AtomIdentity

REVISION = "85767b811c40ed46e73a9b39519cf6bfca8701ba"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, help="Patched pinned Protenix checkout")
    args = parser.parse_args()
    source = args.source.resolve()
    revision = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if revision != REVISION:
        raise SystemExit(f"unexpected Protenix revision: {revision}")

    import torch

    for name in ("protenix", "protenix.model", "protenix.tfg", "protenix.utils"):
        sys.modules[name] = types.ModuleType(name)
    model_utils = types.ModuleType("protenix.model.utils")
    model_utils.centre_random_augmentation = (
        lambda *, x_input_coords, N_sample: x_input_coords.unsqueeze(-3)
    )
    sys.modules["protenix.model.utils"] = model_utils
    sys.modules["protenix.tfg"].parse_tfg_config = lambda _cfg: types.SimpleNamespace(
        enable=False
    )
    sys.modules["protenix.tfg"].TFGEngine = object
    logger_module = types.ModuleType("protenix.utils.logger")
    logger_module.get_logger = lambda _name: types.SimpleNamespace(
        info=lambda *_args: None
    )
    sys.modules["protenix.utils.logger"] = logger_module

    generator_path = source / "protenix/model/generator.py"
    spec = importlib.util.spec_from_file_location("audited_generator", generator_path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load {generator_path}")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)

    identities = tuple(AtomIdentity("A", str(index + 1), "", "CA") for index in range(5))
    fixed_array = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    fixed_coordinates = {
        identities[index]: fixed_array[index].copy() for index in range(3)
    }
    reinjector = FixedAtomReinjector(identities, fixed_coordinates)
    fixed = torch.as_tensor(fixed_array)
    observed: list[tuple[int, object]] = []

    def reinject(state: object, step_index: int, step_count: int) -> object:
        updated = reinjector(state, step_index, step_count)
        observed.append((step_index, updated.clone()))
        return updated

    result = generator.sample_diffusion(
        denoise_net=lambda **kwargs: torch.zeros_like(kwargs["x_noisy"]),
        input_feature_dict={"atom_to_token_idx": torch.arange(5)},
        s_inputs=torch.zeros((1, 1)),
        s_trunk=torch.zeros((1, 1)),
        z_trunk=torch.zeros((1, 1)),
        pair_z=torch.zeros((1, 1)),
        p_lm=torch.zeros((1, 1)),
        c_l=torch.zeros((1, 1)),
        noise_schedule=torch.tensor([1.0, 0.5, 0.0]),
        step_callback=reinject,
    )
    if [index for index, _state in observed] != [0, 1]:
        raise AssertionError("callback did not run after every denoising update")
    if not all(torch.equal(state[0, :3], fixed) for _index, state in observed):
        raise AssertionError("fixed coordinates were not reinserted exactly")
    if not torch.equal(result[0, :3], fixed):
        raise AssertionError("final fixed coordinates changed")
    print("patched Protenix callback: 2/2 exact reinjections")


if __name__ == "__main__":
    main()
