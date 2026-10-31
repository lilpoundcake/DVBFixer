"""Independent bounded audit of complete private completion sampling evidence."""

from __future__ import annotations

import hashlib
import json
import math
import os
import zipfile
from dataclasses import asdict, dataclass
from io import BytesIO
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from dvbfixer.model.diffusion import cc94
from dvbfixer.model.diffusion.atom_completion_contract import AtomCompletionRequest
from dvbfixer.model.diffusion.atom_completion_refinement import COMPLETION_REFINEMENT_REVISION
from dvbfixer.model.diffusion.atom_completion_sampling import (
    assess_completion_sampling,
    completion_engine_view,
)
from dvbfixer.model.diffusion.atom_completion_trace import (
    CompletionSamplerTrace,
    completion_conditioning,
)
from dvbfixer.model.diffusion.atom_completion_validation import (
    CompletionValidation,
    _coordinates,
    audit_completion_candidate,
    validate_completion_candidate,
)
from dvbfixer.model.diffusion.boundary_refinement import BOUNDARY_REFINEMENT_REVISION
from dvbfixer.model.diffusion.contract import ArtifactReference, DiffusionContractError
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.preflight import RunnerPreflightReport
from dvbfixer.model.diffusion.private_artifact import (
    read_private_artifact_bytes,
    read_private_file_bytes,
)
from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER


@dataclass(frozen=True, slots=True)
class CompletionEvidenceAudit:
    request: AtomCompletionRequest
    payloads: dict[str, bytes]
    raw_validation: CompletionValidation
    refined_validation: CompletionValidation


def attest_completion_run(workspace: Path, adapter_source: bytes) -> CompletionSamplerTrace:
    """Supplement retained legacy/private evidence without rewriting any input.

    The provided adapter snapshot must hash to the source recorded at sampling.
    New attestations never adopt or overwrite an existing attempt/trace.
    """
    if workspace.is_symlink() or not workspace.is_dir() or len(adapter_source) > 128_000:
        raise ValueError("completion attestation requires a private run and bounded adapter snapshot")
    root = workspace.resolve(strict=True)
    output = root / "completion-sampler-trace.json"
    adapter_path = root / "completion-adapter-source.py"
    if output.exists() or output.is_symlink() or adapter_path.exists() or adapter_path.is_symlink():
        raise FileExistsError("completion attestation never overwrites retained evidence")
    names = {"request": "request.json", "source": "source.pdb", "motif": "motif.pdb", "reader": "reader.npz",
             "native": "native.npz", "raw": "raw.pdb", "refined": "refined.pdb", "engine_view": "engine-view.json",
             "preflight": "preflight.json", "refinement": "refinement.json", "progress": "trace.json"}
    data = {role: read_private_file_bytes(root, name, max_bytes=2_000_000) for role, name in names.items()}
    progress = json.loads(data["progress"])
    adapter_digest = hashlib.sha256(adapter_source).hexdigest()
    if progress.get("adapter_sha256") != adapter_digest:
        raise DiffusionContractError("completion adapter snapshot does not match the executed source digest")
    references = {role: ArtifactReference(names[role], hashlib.sha256(payload).hexdigest()) for role, payload in data.items()}
    trace = CompletionSamplerTrace.from_dict({"schema_version": 1,
        **{role: asdict(reference) for role, reference in references.items()},
        "adapter_source": asdict(ArtifactReference(adapter_path.name, adapter_digest)),
        "seed": progress["seed"], "device": progress["device"]})
    with adapter_path.open("xb") as stream:
        stream.write(adapter_source)
        stream.flush()
        os.fsync(stream.fileno())
    audit_completion_evidence(trace, root)
    with output.open("xb") as stream:
        stream.write(trace.to_json().encode("utf-8"))
        stream.flush()
        os.fsync(stream.fileno())
    return trace


def audit_completion_evidence(trace: CompletionSamplerTrace, workspace: Path) -> CompletionEvidenceAudit:
    """Recompute transport, source/native/raw linkage and final geometry.

    Digest binding is not cryptographic proof of execution on a hostile host.
    It does prevent accidental substitution or accepting a report without its
    bounded, independently checked source/native/raw/refined evidence.
    """
    if workspace.is_symlink() or not workspace.is_dir():
        raise DiffusionContractError("completion evidence workspace must be a regular directory")
    root = workspace.resolve(strict=True)
    payloads = {role: read_private_artifact_bytes(root, artifact, max_bytes=2_000_000)
                for role, artifact in trace.artifacts().items()}
    request = AtomCompletionRequest.from_json(payloads["request"].decode("utf-8"))
    source, raw, refined = (payloads[role] for role in ("source", "raw", "refined"))
    if request.seeds != (trace.seed,) or request.normalized_pdb.sha256.lower() != trace.source.sha256.lower():
        raise DiffusionContractError("completion trace seed/source differs from frozen request")
    admission = assess_completion_sampling(request, source)
    if not admission.supported:
        raise DiffusionContractError("completion evidence source is outside admitted private scope: " + ", ".join(admission.reasons))
    view = completion_engine_view(request, source)
    if payloads["motif"] != view.motif_text.encode("utf-8"):
        raise DiffusionContractError("completion motif differs from independent source transport")
    if json.loads(payloads["engine_view"]) != json.loads(json.dumps(asdict(view))):
        raise DiffusionContractError("completion engine view differs from independent source mapping")
    preflight = RunnerPreflightReport.from_json(payloads["preflight"].decode("utf-8"))
    mps = trace.device == "mps"
    facts = preflight.facts
    if (not preflight.passed or preflight.profile != (cc94.MPS_PROFILE if mps else cc94.CUDA_PROFILE)
            or preflight.runner_protocol_version != 4 or preflight.contract_schema_version != 4
            or facts.checkpoint_sha256 != trace.checkpoint_sha256 or facts.engine_revision != trace.engine_revision
            or facts.effective_device != trace.device or not facts.accelerator_available or not facts.fallback_disabled
            or facts.environment_identity != (cc94.MPS_ENVIRONMENT if mps else cc94.CUDA_ENVIRONMENT)
            or facts.framework_version != "2.6.0" or facts.python_version != "3.12"
            or facts.platform != ("darwin" if mps else "linux")
            or facts.architecture not in ({"arm64"} if mps else {"x86_64", "amd64"})
            or facts.patch_identity != "cc94-native-relative-v1;source-physical-v1;" + ("apple-patch-v1" if mps else "callback-patch-v1")
            or facts.refinement_platform != ("OpenMM CPU" if mps else "OpenMM Reference")):
        raise DiffusionContractError("completion preflight does not attest the pinned native runtime")
    progress = json.loads(payloads["progress"])
    required = {"schema_version": 1, "kind": "private-cc94-atom-completion-trace", "status": "accepted-private",
                "seed": trace.seed, "device": trace.device, "ode_updates": 500, "step_scale": 1.0,
                "s_churn": 0.0, "chain_phase": "native", "checkpoint_sha256": trace.checkpoint_sha256,
                "config_sha256": trace.config_sha256, "engine_revision": trace.engine_revision,
                "request_sha256": trace.request.sha256, "source_sha256": trace.source.sha256,
                "motif_sha256": trace.motif.sha256, "native_sha256": trace.native.sha256,
                "raw_sha256": trace.raw.sha256, "adapter_sha256": trace.adapter_source.sha256,
                "conditioning": completion_conditioning(view.fixed_positions)}
    if any(progress.get(key) != value for key, value in required.items()):
        raise DiffusionContractError("completion progress trace disagrees with bound evidence or frozen schedule")
    if any(type(progress.get(key)) is not int for key in ("schema_version", "seed", "ode_updates")):
        raise DiffusionContractError("completion progress counters require integer types")
    if (progress.get("checkpoint_sampling_executed") is not True or progress.get("checkpoint_sampling_started") is not True
            or any(progress.get(key) is not False for key in ("published", "scientific_qualification", "public_profile_enabled"))):
        raise DiffusionContractError("completion progress has invalid execution/qualification flags")
    reader = _arrays(payloads["reader"], {"atom_mask", "atom_positions", "residue_index"})
    count = len(view.fixed_positions)
    if (reader["atom_mask"].shape != (count, 37) or reader["atom_positions"].shape != (count, 37, 3)
            or not np.array_equal(reader["residue_index"], np.asarray(view.fixed_positions) + 1)):
        raise DiffusionContractError("completion reader axes do not match observed target positions")
    expected_mask = np.zeros((count, 37))
    target_index = dict(zip(view.atom_axis, view.atom_target_indices, strict=True))
    row_index = {index: row for row, index in enumerate(view.fixed_positions)}
    authoritative = _coordinates(source.decode("utf-8").splitlines())
    for atom in view.represented_fixed_atoms:
        row, column = row_index[target_index[atom]], ATOM37_ORDER.index(atom.atom_name)
        expected_mask[row, column] = 1
        if np.linalg.norm(reader["atom_positions"][row, column] - authoritative[atom]) > 1e-5:
            raise DiffusionContractError("completion reader moved an observed source coordinate")
    if not np.array_equal(expected_mask, reader["atom_mask"]):
        raise DiffusionContractError("completion reader changed fixed/partial/gap masks")
    native = _arrays(payloads["native"], {"atom37"})["atom37"]
    if native.dtype.kind != "f" or native.shape != (len(view.target_residues), 37, 3):
        raise DiffusionContractError("completion native tensor shape differs from the full target axis")
    coordinates = np.asarray([native[index, ATOM37_ORDER.index(atom.atom_name)]
                              for atom, index in zip(view.atom_axis, view.atom_target_indices, strict=True)], dtype=np.float64)
    by_atom = dict(zip(view.atom_axis, coordinates, strict=True))
    transform = weighted_kabsch(np.asarray([by_atom[atom] for atom in view.represented_fixed_atoms]),
                                np.asarray([authoritative[atom] for atom in view.represented_fixed_atoms]))
    fitted = dict(zip(view.atom_axis, transform.apply(coordinates), strict=True))
    raw_audit = audit_completion_candidate(source, raw, request)
    expected = np.asarray([[float(f"{value:.3f}") for value in fitted[atom]] for atom in request.generated_atoms])
    if not np.array_equal(expected, np.asarray(raw_audit.generated_coordinates)):
        raise DiffusionContractError("completion raw coordinates are not materialized from the bound native sample")
    refinement = json.loads(payloads["refinement"])
    expected_refinement = {"status": "accepted-private", "revision": COMPLETION_REFINEMENT_REVISION,
                           "kernel_revision": BOUNDARY_REFINEMENT_REVISION, "seed": trace.seed,
                           "platform": "CPU" if mps else "Reference", "source_sha256": trace.source.sha256,
                           "raw_sha256": trace.raw.sha256, "refined_sha256": trace.refined.sha256}
    if (any(refinement.get(key) != value for key, value in expected_refinement.items())
            or any(refinement.get(key) is not False for key in ("published", "scientific_qualification"))):
        raise DiffusionContractError("completion refinement evidence differs from bound source/raw/refined artifacts")
    for key in ("initial_energy_kj_mol", "final_energy_kj_mol"):
        if type(refinement.get(key)) not in {int, float} or not math.isfinite(refinement[key]):
            raise DiffusionContractError("completion refinement evidence requires finite energies")
    raw_validation = validate_completion_candidate(source, raw, request)
    final_validation = validate_completion_candidate(source, refined, request)
    if not final_validation.summary.passed:
        raise DiffusionContractError("completion final geometry failed independent validation: "
                                     + ", ".join(final_validation.summary.hard_gate_failures))
    return CompletionEvidenceAudit(request, payloads, raw_validation, final_validation)


def _arrays(data: bytes, names: set[str]) -> dict[str, NDArray[np.generic]]:
    """Reject compressed bombs, object arrays, extra members and nonfinite data."""
    with zipfile.ZipFile(BytesIO(data)) as archive:
        members = archive.infolist()
        if ({member.filename for member in members} != {name + ".npy" for name in names}
                or len(members) != len(names) or sum(member.file_size for member in members) > 1_000_000):
            raise DiffusionContractError("completion tensor archive has unexpected members or exceeds its bound")
        for member in members:
            with archive.open(member) as stream:
                version = np.lib.format.read_magic(stream)
                if version not in {(1, 0), (2, 0)}:
                    raise DiffusionContractError("unsupported completion array header version")
                shape, _, dtype = (np.lib.format.read_array_header_1_0(stream) if version == (1, 0)
                                   else np.lib.format.read_array_header_2_0(stream))
                count = math.prod(shape)
                if (len(shape) > 3 or any(size < 0 for size in shape) or count > 128 * 37 * 3
                        or dtype.hasobject or dtype.kind not in "fiu" or dtype.itemsize > 8
                        or stream.tell() + count * dtype.itemsize != member.file_size):
                    raise DiffusionContractError("completion array header exceeds its numeric/shape/byte bounds")
    with np.load(BytesIO(data), allow_pickle=False) as archive:
        result = {name: archive[name] for name in names}
    if any(array.dtype.kind not in "fiu" or not np.isfinite(array).all() for array in result.values()):
        raise DiffusionContractError("completion tensor evidence must be finite numeric arrays")
    return result
