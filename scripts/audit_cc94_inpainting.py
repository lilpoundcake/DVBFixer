"""Independent NumPy audit of every retained known-atom projection boundary."""
from __future__ import annotations

import json
from collections import Counter

import numpy as np
from scripts.mac_reconstruction_comparison import digest

REVISION = "private-cc94-fixed-scaffold-ve-inpainting-v1"


def audit_trajectory(directory, expected_known, expected_generated, native, *, steps=500, tolerance=1.e-4):
    work = directory / "inpainting"
    summary = json.loads((directory / "inpainting-summary.json").read_text())
    if (summary["revision"] != REVISION or summary["accepted_steps"] != steps
            or summary["denoiser_inputs"] != steps or summary["clean_predictions"] != steps
            or summary["maximum_post_projection_error_angstrom"] != 0):
        raise ValueError("inpainting process did not retain a complete exact projection trace")
    with np.load(work / "reference.npz", allow_pickle=False) as data:
        clean, epsilon = data["clean"].copy(), data["epsilon"].copy()
        known, generated = data["known_mask"].copy(), data["generated_mask"].copy()
    if (not np.array_equal(known, expected_known) or not np.array_equal(generated, expected_generated)
            or (known & generated).any() or clean.shape != (*known.shape, 3)
            or epsilon.shape != clean.shape or not np.isfinite(clean).all() or not np.isfinite(epsilon).all()):
        raise ValueError("trajectory mask/reference inventory differs from source-bound ownership")
    events = [json.loads(line) for line in (work / "events.jsonl").read_text().splitlines()]
    counts = Counter(event["stage"] for event in events)
    for stage in ("denoiser-input", "clean-prediction", "accepted-step"):
        selected = [e for e in events if e["stage"] == stage]
        if len(selected) != steps or [e["step"] for e in selected] != list(range(steps)):
            raise ValueError("missing/reordered native inpainting events")
    if any(not e["nonfixed_bits_unchanged"] or e["post_projection_error_angstrom"] != 0 for e in events):
        raise ValueError("projection changed generated/unowned atoms or failed exact assignment")
    maximum = 0.
    sigma_after = []
    checked = []
    for stage in ("initial", "denoiser-input", "clean-prediction", "accepted-step"):
        indices = [-1] if stage == "initial" else range(steps)
        for index in indices:
            path = work / f"projection-{stage}-{index:04d}.npz"
            with np.load(path, allow_pickle=False) as data:
                before, after = data["known_before"], data["known_after"]
                first, second = data["generated_before"], data["generated_after"]
                sigma = float(data["sigma"])
            if before.shape != (int(known.sum()), 3) or after.shape != before.shape or first.shape != (int(generated.sum()), 3):
                raise ValueError("retained projection arrays have an inconsistent ownership axis")
            if not np.array_equal(first, second):
                raise ValueError("projection changed a generated coordinate")
            if stage == "clean-prediction":
                expected = clean[known]
            else:
                expected = (clean+np.asarray(sigma, dtype=clean.dtype)*epsilon)[known]
            error = float(np.linalg.norm(after-expected, axis=-1).max())
            maximum = max(maximum, error)
            if not np.isfinite(error) or error > tolerance:
                raise ValueError("retained known coordinates violate the declared VE noise path")
            if stage == "initial":
                if sigma <= 0 or not np.allclose(before/np.float32(sigma), epsilon[known], rtol=1.e-6, atol=1.e-6):
                    raise ValueError("known-noise epsilon was not the original native initialized noise")
            if stage == "accepted-step":
                sigma_after.append(sigma)
            checked.append((str(path), digest(path)))
    if sigma_after[-1] != 0 or any(a < b for a, b in zip(sigma_after, sigma_after[1:])):
        raise ValueError("accepted known-noise schedule is not decreasing to zero")
    if not np.array_equal(native[expected_known[0]], clean[expected_known]):
        raise ValueError("returned native endpoint differs from the zero-noise fixed reference")
    return {"revision": REVISION, "steps": steps, "event_counts": dict(counts),
        "projection_payloads_checked": len(checked), "known_atoms": int(known.sum()),
        "generated_atoms": int(generated.sum()), "generated_projection_bits_exact": True,
        "numpy_noise_path_max_error_angstrom": maximum, "noise_path_tolerance_angstrom": tolerance,
        "native_endpoint_known_bits_exact": True, "source_copy_after_sampling": False,
        "reference_sha256": digest(work / "reference.npz"),
        "payload_sha256": dict(checked)}
