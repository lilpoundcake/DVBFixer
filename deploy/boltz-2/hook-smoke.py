"""Exercise the maintained Boltz-2 post-update callback on its real sampler."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import torch
from torch import nn

PINNED_REVISION = "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"


def _revision(checkout: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=checkout,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def run(checkout: Path, *, device: str) -> dict[str, object]:
    checkout = checkout.resolve()
    if _revision(checkout) != PINNED_REVISION:
        raise ValueError("Boltz checkout revision does not match the maintained patch")

    from boltz.model.modules.diffusionv2 import AtomDiffusion

    class SmokeDiffusion(AtomDiffusion):
        def __init__(self) -> None:
            nn.Module.__init__(self)
            self.probe = nn.Parameter(torch.zeros(()))
            self.sigma_min = 0.0004
            self.sigma_max = 1.0
            self.sigma_data = 1.0
            self.rho = 7.0
            self.gamma_0 = 0.0
            self.gamma_min = 1.0
            self.noise_scale = 1.0
            self.step_scale = 1.0
            self.step_scale_random = None
            self.num_sampling_steps = 2
            self.alignment_reverse_diff = False

        @property
        def device(self) -> torch.device:
            return self.probe.device

        def preconditioned_network_forward(
            self,
            noised_atom_coords: torch.Tensor,
            sigma: float,
            network_condition_kwargs: dict[str, object],
        ) -> torch.Tensor:
            del sigma, network_condition_kwargs
            return torch.zeros_like(noised_atom_coords)

    target = torch.tensor([[1.0, 2.0, 3.0], [-2.0, 0.5, 4.0]])
    steps: list[dict[str, object]] = []

    def callback(state: torch.Tensor, step_index: int, step_count: int) -> torch.Tensor:
        updated = state.clone()
        projected = target.to(device=state.device, dtype=state.dtype)
        updated[:, :2] = projected
        error = torch.max(torch.abs(updated[:, :2] - projected)).item()
        steps.append(
            {
                "step_index": step_index,
                "step_count": step_count,
                "post_projection_max_error": error,
            }
        )
        return updated

    steering_args = {
        "fk_steering": False,
        "physical_guidance_update": False,
        "contact_guidance_update": False,
    }
    diffusion = SmokeDiffusion().to(device)
    result = diffusion.sample(
        atom_mask=torch.ones((1, 3)),
        num_sampling_steps=2,
        multiplicity=1,
        steering_args=steering_args,
        step_callback=callback,
    )["sample_atom_coords"]
    if [step["step_index"] for step in steps] != [0, 1]:
        raise RuntimeError("callback was not invoked after every update")
    if not torch.equal(result[:, :2].cpu(), target.unsqueeze(0)):
        raise RuntimeError("final fixed coordinates were not preserved exactly")

    return {
        "boltz_revision": PINNED_REVISION,
        "torch_version": torch.__version__,
        "device": str(diffusion.device),
        "callback_count": len(steps),
        "expected_callback_count": 2,
        "maximum_post_projection_error": max(
            float(step["post_projection_max_error"]) for step in steps
        ),
        "final_fixed_coordinates_exact": True,
        "steps": steps,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkout", type=Path)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    print(json.dumps(run(args.checkout, device=args.device), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
