"""Deterministic CPU-only external runner used by diffusion protocol tests."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any

from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    BackendProvenance,
    DiffusionRequest,
    DiffusionStatus,
    RunnerCandidate,
    RunnerDiagnostics,
    RunnerResult,
)
from dvbfixer.model.diffusion.runner import (
    DIFFUSION_RUNNER_PROTOCOL_VERSION,
    REQUEST_MANIFEST,
    RESULT_MANIFEST,
)


def main() -> int:
    request_name = os.environ.get("DVBFIXER_DIFFUSION_REQUEST", REQUEST_MANIFEST)
    result_name = os.environ.get("DVBFIXER_DIFFUSION_RESULT", RESULT_MANIFEST)
    request = DiffusionRequest.from_json(Path(request_name).read_text(encoding="utf-8"))

    input_path = Path(request.normalized_pdb.path)
    candidates: list[RunnerCandidate] = []
    for index, seed in enumerate(request.seeds, start=1):
        candidate_id = f"candidate-{index:04d}"
        output_path = Path("candidates") / f"{candidate_id}.pdb"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input_path, output_path)
        candidates.append(
            RunnerCandidate(
                candidate_id=candidate_id,
                coordinate_artifact=ArtifactReference(
                    path=output_path.as_posix(),
                    sha256=_sha256(output_path),
                ),
                generated_atoms=request.generated_atoms,
                generated_residues=tuple(
                    residue
                    for gap in request.gaps
                    for residue in gap.generated_residues
                ),
                raw_backend_score=float(seed),
                score_provenance="dvbfixer-fake-runner:seed-v1",
            )
        )

    result = RunnerResult(
        schema_version=DIFFUSION_SCHEMA_VERSION,
        status=DiffusionStatus.SUCCESS,
        candidates=tuple(candidates),
        runner_diagnostics=RunnerDiagnostics(exit_code=None, timed_out=False),
        backend_provenance=BackendProvenance(
            backend="dvbfixer-fake",
            runner_protocol_version=DIFFUSION_RUNNER_PROTOCOL_VERSION,
            engine_repository="builtin://dvbfixer",
            engine_revision="fake-runner-v1",
            environment_hash=_environment_hash(),
        ),
    )
    Path(result_name).write_text(result.to_json(), encoding="utf-8")
    return 0


def _environment_hash() -> str:
    payload: dict[str, Any] = {
        "contract_schema_version": DIFFUSION_SCHEMA_VERSION,
        "runner_protocol_version": DIFFUSION_RUNNER_PROTOCOL_VERSION,
        "implementation": "dvbfixer-fake-runner-v1",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
