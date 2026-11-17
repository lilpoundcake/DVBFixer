"""Private explicit MPS cache-release hooks; never change live model tensors."""
from __future__ import annotations

import json

import torch


def _mps_tensors(value):
    if isinstance(value, torch.Tensor):
        yield value.device.type == "mps"
    elif isinstance(value, dict):
        for item in value.values():
            yield from _mps_tensors(item)
    elif isinstance(value, (tuple, list)):
        for item in value:
            yield from _mps_tensors(item)


def install_boundary(module, log_path, *, label):
    """Attach no-output-changing hooks to one explicitly selected eval module."""
    def release(phase, values):
        if not any(_mps_tensors(values)):
            return
        if module.training or torch.is_grad_enabled():
            raise ValueError("private MPS cache release requires eval/no_grad inference")
        before = {"driver": torch.mps.driver_allocated_memory(), "live": torch.mps.current_allocated_memory()}
        torch.mps.synchronize()
        torch.mps.empty_cache()
        after = {"driver": torch.mps.driver_allocated_memory(), "live": torch.mps.current_allocated_memory()}
        with log_path.open("a") as stream:
            stream.write(json.dumps({"module": label, "phase": phase, "before": before, "after": after}) + "\n")

    def pre(_module, args, kwargs):
        release("pre", (args, kwargs))

    def post(_module, args, kwargs, result):
        release("post", result)

    return (module.register_forward_pre_hook(pre, with_kwargs=True),
            module.register_forward_hook(post, with_kwargs=True))
