"""Private, passive sampler-return evidence; never bypass downstream model work."""
from __future__ import annotations

import json
import os
from contextlib import contextmanager
from dataclasses import asdict

import numpy as np
from scripts.mac_reconstruction_comparison import digest

from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb


def write_json(path, value):
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
        stream.flush()
        os.fsync(stream.fileno())


def capture_return(output, tensor, *, request, request_path, source, axis, names, elements, steps):
    """Commit complete raw evidence before returning the exact original tensor."""
    if (len(steps) != 200 or [s["step_index"] for s in steps] != list(range(200))
            or any(s["step_count"] != 200 or s["post_overwrite_max_angstrom"] != 0 for s in steps)):
        raise ValueError("sampler return requires all200 exact reinjection callbacks")
    if tensor.device.type != "mps" or str(tensor.dtype) != "torch.float32":
        raise ValueError("capture requires native MPS float32 sampler return")
    if len(set(axis)) != len(axis) or not len(axis) == len(names) == len(elements):
        raise ValueError("native identity axis is not unique/complete")
    xyz = tensor.detach().cpu().numpy()
    if xyz.shape != (1, len(axis), 3) or not np.isfinite(xyz).all():
        raise ValueError("sampler coordinates have unexpected shape/nonfinite values")
    if digest(source) != request.normalized_pdb.sha256:
        raise ValueError("source changed before capture")
    selected = {item.atom for item in request.protein_context.fixed_conditioning} | set(request.generated_atoms)
    if not selected <= set(axis):
        raise ValueError("native sampler axis omitted requested atoms")
    indices = tuple(i for i, atom in enumerate(axis) if atom in selected)
    native = xyz[0].astype(np.float64)
    text = materialize_candidate_pdb(source.read_text(), tuple(axis[i] for i in indices), native[list(indices)],
        tuple(names[i] for i in indices), request, elements=tuple(elements[i] for i in indices))
    output.mkdir(exist_ok=False)
    with (output / "native-arrays.npz").open("xb") as stream:
        np.savez_compressed(stream, coordinates=native)
        stream.flush()
        os.fsync(stream.fileno())
    write_json(output / "native-axis.json", {"atoms": [asdict(atom) for atom in axis],
        "residue_names": list(names), "elements": list(elements), "materialized_indices": list(indices)})
    write_json(output / "reinjection-steps.json", steps)
    with (output / "raw-materialized.pdb").open("x") as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    payloads = sorted(output.iterdir())
    write_json(output / "capture.json", {"kind": "private-protenix-sampler-return-v1",
        "stage": "sample_diffusion returned; downstream heads not yet executed",
        "request_path": str(request_path), "request_sha256": digest(request_path),
        "source_path": str(source), "source_sha256": digest(source), "seed": request.seeds[0],
        "native_atoms": len(axis), "materialized_atoms": len(indices), "updates": len(steps),
        "files": {p.name: digest(p) for p in payloads}, "downstream_work_skipped": False,
        "scientific_acceptance": False, "native_full_inference_success": False,
        "coordinate_tensor_returned_unchanged": True})
    return tensor


@contextmanager
def passive_capture(legacy, model_class, *, output, request, request_path, source):
    original_axis, original_sample = legacy._model_atom_axis, model_class.sample_diffusion
    state = {}

    def record_axis(array, current_request):
        if current_request != request:
            raise ValueError("captured request differs from native request")
        result = original_axis(array, current_request)
        if state:
            raise ValueError("unexpected second native identity axis")
        state["axis"], state["names"], state["elements"] = result
        return result

    def sample(model, **kwargs):
        result = original_sample(model, **kwargs)
        return capture_return(output, result, request=request, request_path=request_path, source=source,
            steps=model.dvbfixer_step_callback.steps, **state)

    legacy._model_atom_axis, model_class.sample_diffusion = record_axis, sample
    try:
        yield
    finally:
        legacy._model_atom_axis, model_class.sample_diffusion = original_axis, original_sample
