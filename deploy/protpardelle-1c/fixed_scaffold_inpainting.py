"""Private VE/ODE known-atom inpainting; no shared engine files are modified.

Known x_sigma = clean + sigma*epsilon, with epsilon reused from initial native
Gaussian noise. Known x0 is clean. Only known entries are projected at initial,
denoiser-input, clean-prediction and accepted-step boundaries. This is NOT the
upstream replacement_guidance option or end-of-run source restoration.
"""
from __future__ import annotations

import ast
import hashlib
import inspect
import json
import re
import textwrap
import types
from contextlib import contextmanager

import numpy as np

REVISION = "private-cc94-fixed-scaffold-ve-inpainting-v1"


class FixedScaffold:
    def __init__(self, known_mask, generated_mask, output, *, expected_steps=500):
        self.known_numpy = np.asarray(known_mask, dtype=bool)
        self.generated_numpy = np.asarray(generated_mask, dtype=bool)
        if (self.known_numpy.ndim != 3 or self.known_numpy.shape[0] != 1
                or self.known_numpy.shape[-1] != 37
                or self.generated_numpy.shape != self.known_numpy.shape
                or (self.known_numpy & self.generated_numpy).any()
                or not self.known_numpy.any() or not self.generated_numpy.any()):
            raise ValueError("inpainting requires disjoint explicit known/generated atom37 masks")
        self.output = output
        output.mkdir(parents=True, exist_ok=False)
        self.expected_steps = expected_steps
        self.steps = 0
        self.inputs = 0
        self.predictions = 0
        self.initialized = False
        self.maximum_post_error = 0.

    def validate_mode(self, *, num_steps, step_scale, s_churn, jump_steps, uniform_steps,
                      gt_aatype, cc, pd, sidechain_mode):
        if (num_steps != self.expected_steps or step_scale != 1. or s_churn != 0.
                or jump_steps or not uniform_steps or gt_aatype is None or not cc.enabled
                or cc.num_recurrence_steps != 1 or cc.reconstruction_guidance.enabled
                or cc.replacement_guidance.enabled or pd.enabled or sidechain_mode):
            raise ValueError("private inpainting supports fixed-sequence uniform VE Euler only; no fallback")

    def initialize(self, state, motif_indices, motif_coordinates, motif_mask, sigma):
        import torch

        if self.initialized or tuple(state.shape[:-1]) != self.known_numpy.shape or state.shape[-1] != 3:
            raise ValueError("invalid/repeated inpainting initialization")
        if motif_indices is None or len(motif_indices) != 1 or motif_mask is None:
            raise ValueError("inpainting requires native atom-level motif evidence")
        self.known = torch.as_tensor(self.known_numpy, device=state.device)
        self.generated = torch.as_tensor(self.generated_numpy, device=state.device)
        self.clean = torch.zeros_like(state)
        represented = torch.zeros_like(self.known)
        positions = motif_indices[0]
        coords = motif_coordinates[0].to(device=state.device, dtype=state.dtype)
        mask = motif_mask[0].to(device=state.device, dtype=torch.bool)
        if len(set(positions)) != len(positions) or len(positions) != len(coords) or mask.shape != coords.shape[:-1]:
            raise ValueError("native motif placement/masks are ambiguous")
        self.clean[0, positions] = coords
        represented[0, positions] = mask
        if not torch.equal(represented, self.known):
            raise ValueError("native motif mask does not equal the declared fixed observations")
        sigma = self._sigma(sigma, state)
        if not bool((sigma > 0).all()):
            raise ValueError("initial noise sigma must be positive")
        # No new RNG draws: reuse the exact native initialized unit Gaussian.
        self.epsilon = state.detach().clone()/sigma
        if not bool(torch.isfinite(self.clean[self.known]).all()) or not bool(torch.isfinite(self.epsilon).all()):
            raise FloatingPointError("nonfinite inpainting reference/noise")
        self.initialized = True
        np.savez_compressed(self.output / "reference.npz", clean=self.clean.detach().cpu().numpy(),
            epsilon=self.epsilon.detach().cpu().numpy(), known_mask=self.known_numpy,
            generated_mask=self.generated_numpy)

    @staticmethod
    def _sigma(value, state):
        import torch

        sigma = torch.as_tensor(value, device=state.device, dtype=state.dtype).reshape(-1, 1, 1, 1)
        if sigma.shape[0] not in {1, state.shape[0]} or not bool(torch.isfinite(sigma).all()) or bool((sigma < 0).any()):
            raise ValueError("invalid inpainting noise level")
        return sigma

    def _project(self, state, target, *, stage, step, sigma=None, save_state=False):
        import torch

        if not self.initialized or state.shape != self.clean.shape or state.device != self.clean.device or state.dtype != self.clean.dtype:
            raise ValueError("inpainting state changed axis/device/dtype")
        if not bool(torch.isfinite(state).all()):
            raise FloatingPointError("nonfinite state before inpainting projection")
        result = state.clone()
        before = float(torch.linalg.vector_norm(state[self.known]-target[self.known], dim=-1).max().item())
        result[self.known] = target[self.known]
        unchanged = torch.equal(result[~self.known], state[~self.known])
        if not unchanged:
            raise RuntimeError("inpainting modified a generated/unowned atom")
        error = float(torch.linalg.vector_norm(result[self.known]-target[self.known], dim=-1).max().item())
        self.maximum_post_error = max(self.maximum_post_error, error)
        event = {"stage": stage, "step": int(step), "pre_projection_error_angstrom": before,
                 "post_projection_error_angstrom": error, "nonfixed_bits_unchanged": unchanged,
                 "known_atoms": int(self.known_numpy.sum()), "generated_atoms": int(self.generated_numpy.sum())}
        if sigma is not None:
            event["sigma"] = float(sigma.reshape(-1)[0].item())
        with (self.output / "events.jsonl").open("a") as stream:
            stream.write(json.dumps(event, allow_nan=False)+"\n")
        if stage != "self-conditioning":
            np.savez_compressed(self.output / f"projection-{stage}-{step:04d}.npz",
                known_before=state[self.known].detach().cpu().numpy(),
                known_after=result[self.known].detach().cpu().numpy(),
                generated_before=state[self.generated].detach().cpu().numpy(),
                generated_after=result[self.generated].detach().cpu().numpy(),
                sigma=0. if sigma is None else float(sigma.reshape(-1)[0].item()))
        if save_state:
            np.savez_compressed(self.output / f"{stage}-{step:04d}.npz",
                before=state.detach().cpu().numpy(), after=result.detach().cpu().numpy(),
                target=target.detach().cpu().numpy(), sigma=0. if sigma is None else float(sigma.reshape(-1)[0].item()))
        return result

    def noisy(self, state, sigma, step, stage="denoiser-input"):
        sigma = self._sigma(sigma, state)
        result = self._project(state, self.clean+sigma*self.epsilon, stage=stage, step=step, sigma=sigma,
                               save_state=stage == "initial")
        if stage == "denoiser-input":
            if step != self.inputs:
                raise ValueError("unexpected denoiser inpainting order")
            self.inputs += 1
        return result

    def prediction(self, state, step):
        if step != self.predictions:
            raise ValueError("unexpected clean prediction inpainting order")
        self.predictions += 1
        if getattr(self, 'placement_chains', None) is not None and step >= 250:
            import torch
            from scipy.spatial import cKDTree

            gr = np.where(self.generated_numpy)[1]
            if not hasattr(self, 'placement_tree'):
                self.fixed_rows = np.where(self.known_numpy)[1]
                self.placement_tree = cKDTree(self.clean[self.known].detach().cpu().numpy())
            before = state[self.generated].detach().cpu().numpy().copy()
            shifts = residue_shifts(before, gr, self.placement_tree, self.fixed_rows, self.placement_chains)
            self.placement_delta = torch.as_tensor(shifts, dtype=state.dtype, device=state.device)
            state = state.clone()
            state[self.generated] += self.placement_delta
            np.savez_compressed(self.output / f'placement-{step:04d}.npz', before=before,
                                after=state[self.generated].detach().cpu().numpy(), shifts=shifts)
        return self._project(state, self.clean, stage="clean-prediction", step=step)

    def self_condition(self, state, sigma_data, step):
        if state is None:
            return None
        if sigma_data <= 0:
            raise ValueError("invalid model self-conditioning scale")
        if hasattr(self, 'placement_delta'):
            state = state.clone()
            state[self.generated] += self.placement_delta/sigma_data
        return self._project(state, self.clean/sigma_data, stage="self-conditioning", step=step)

    def after_step(self, state, sigma, step, step_count):
        if step != self.steps or step_count != self.expected_steps:
            raise ValueError("unexpected accepted-step inpainting order/count")
        sigma = self._sigma(sigma, state)
        if step == step_count-1 and bool((sigma != 0).any()):
            raise ValueError("final inpainting noise is not zero")
        result = self._project(state, self.clean+sigma*self.epsilon, stage="accepted-step", step=step, sigma=sigma,
                               save_state=step < 3 or (step+1) % 50 == 0 or step == step_count-1)
        self.steps += 1
        return result

    def output_mask(self, mask):
        if self.steps != self.expected_steps or self.inputs != self.expected_steps or self.predictions != self.expected_steps:
            raise ValueError("incomplete native inpainting trajectory")
        if mask.shape != self.known.shape:
            raise ValueError("native output mask axis differs")
        # Input-only OXT is carried by the native process, NOT model-generated.
        return mask.bool() | self.known

    def summary(self):
        return {"revision": REVISION, "accepted_steps": self.steps, "denoiser_inputs": self.inputs,
                "clean_predictions": self.predictions, "expected_steps": self.expected_steps,
                "maximum_post_projection_error_angstrom": self.maximum_post_error,
                "known_atoms": int(self.known_numpy.sum()), "generated_atoms": int(self.generated_numpy.sum()),
                "noise_policy": "VE fixed native initial epsilon; x_known=clean+sigma*epsilon; x0_known=clean",
                "source_restoration_after_sampling": False}


def residue_shifts(points, rows, tree, fixed_rows, chains):
    """Translate each generated residue away from its closest eligible scaffold contact."""
    shifts = np.zeros_like(points)
    closest = {}
    for i, hits in enumerate(tree.query_ball_point(points, 3.0)):
        row = rows[i]
        for j in hits:
            other = fixed_rows[j]
            if chains[row] == chains[other] and abs(row-other) <= 1:
                continue
            delta = points[i]-tree.data[j]
            distance = np.linalg.norm(delta)
            if row not in closest or distance < closest[row][0]:
                closest[row] = (distance, delta)
    for row, (distance, delta) in closest.items():
        if distance > 1.e-8:
            shifts[rows == row] = min(.35, .5*(3.-distance))*delta/distance
    return shifts


def patched_source(original):
    """Exact, bounded in-memory edits to the verified native sampling method."""
    source = textwrap.dedent(inspect.getsource(original))
    marker = "    # Sampling trajectory\n"
    if source.count(marker) != 1:
        raise ValueError("unsupported sampler initialization location")
    setup = """    _dvbfixer_inpaint.validate_mode(num_steps=num_steps, step_scale=step_scale,
        s_churn=s_churn, jump_steps=jump_steps, uniform_steps=uniform_steps,
        gt_aatype=gt_aatype, cc=cc, pd=pd, sidechain_mode=sidechain_mode)
    _dvbfixer_inpaint.initialize(xt, motif_idx, motif_all_atom, motif_atom_mask, sigma)
    xt = _dvbfixer_inpaint.noisy(xt, sigma, -1, stage="initial")
"""
    source = source.replace(marker, setup+marker)
    pattern = r"(?m)^( +)xt_rep = xt_hat\.clone\(\)$"
    if len(re.findall(pattern, source)) != 2:
        raise ValueError("unsupported sampler denoiser input locations")
    source = re.sub(pattern, lambda m: m[1]+"xt_hat = _dvbfixer_inpaint.noisy(xt_hat, sigma_hat if gamma > 0 else sigma, i)\n"+m[0], source)
    tree = ast.parse(source)
    predictions = [n for n in ast.walk(tree) if isinstance(n, ast.Assign)
                   and isinstance(n.value, ast.Call) and isinstance(n.value.func, ast.Attribute)
                   and isinstance(n.value.func.value, ast.Name) and n.value.func.value.id == "self"
                   and n.value.func.attr == "forward" and isinstance(n.targets[0], ast.Tuple)
                   and [a.id for a in n.targets[0].elts if isinstance(a, ast.Name)]
                   == ["x0", "s_logprobs", "x_self_cond", "s_self_cond"]]
    if len(predictions) != 2:
        raise ValueError("unsupported sampler prediction locations")
    lines = source.splitlines(keepends=True)
    for node in sorted(predictions, key=lambda n: n.end_lineno, reverse=True):
        indent = " "*node.col_offset
        lines.insert(node.end_lineno, indent+"x0 = _dvbfixer_inpaint.prediction(x0, i)\n"+
                     indent+"x_self_cond = _dvbfixer_inpaint.self_condition(x_self_cond, self.sigma_data, i)\n")
    source = "".join(lines)
    marker = "        sigma = sigma_next\n"
    if source.count(marker) != 1:
        raise ValueError("unsupported sampler accepted-step location")
    source = source.replace(marker, "        xt = _dvbfixer_inpaint.after_step(xt, sigma_next, i, num_steps)\n"+marker)
    marker = r"(?m)^    atom_mask = atom37_mask_from_aatype\(s_hat, seq_mask\)$"
    if len(re.findall(marker, source)) != 1:
        raise ValueError("unsupported sampler output-mask location")
    return re.sub(marker, lambda m: m[0]+"\n    atom_mask = _dvbfixer_inpaint.output_mask(atom_mask)", source)


@contextmanager
def installed(model, controller):
    original = model.sample
    source = patched_source(original)
    namespace = dict(original.__func__.__globals__)
    namespace["_dvbfixer_inpaint"] = controller
    exec(compile(source, "<private-cc94-fixed-scaffold-inpainting>", "exec"), namespace)
    (controller.output / "executed-sampler.py").write_text(source)
    (controller.output / "patch.json").write_text(json.dumps({"revision": REVISION,
        "original_method_sha256": hashlib.sha256(inspect.getsource(original).encode()).hexdigest(),
        "patched_method_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "shared_engine_files_modified": False}, indent=2)+"\n")
    had_override = "sample" in model.__dict__
    previous_override = model.__dict__.get("sample")
    model.sample = types.MethodType(namespace["sample"], model)
    try:
        yield
    finally:
        if had_override:
            model.sample = previous_override
        else:
            del model.sample
