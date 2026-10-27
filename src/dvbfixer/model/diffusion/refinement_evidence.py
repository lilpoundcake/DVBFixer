"""Bounded artifact custody for staged refinement, independent of runner claims.

This record grants no geometry or mobility acceptance. Consumers must verify
every artifact digest and independently recompute eligibility and geometry.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, replace
from typing import Literal

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    BackendOption,
    DiffusionContractError,
    DiffusionRequest,
    SamplerTrace,
)
from dvbfixer.model.diffusion.refinement_policy import read_refinement_policy
from dvbfixer.model.diffusion.refinement_stages import NumericalFailureEvidence

MAX_REFINEMENT_EVIDENCE_BYTES = 16_384
MAX_REFINEMENT_REPORT_BYTES = 16_000_000
EVIDENCE_REVISION = "staged-refinement-evidence-v1"
_PATH_KEY = "staged_evidence_path"
_SHA_KEY = "staged_evidence_sha256"


def evidence_reference(trace: SamplerTrace) -> ArtifactReference | None:
    selected = [option for option in trace.refinement_parameters if option.name in {_PATH_KEY, _SHA_KEY}]
    values = {option.name: option.value for option in selected}
    if not selected:
        return None
    if len(selected) != 2 or set(values) != {_PATH_KEY, _SHA_KEY}:
        raise DiffusionContractError("incomplete or duplicate staged evidence reference")
    return ArtifactReference(values[_PATH_KEY], values[_SHA_KEY])


def bind_evidence(trace: SamplerTrace, reference: ArtifactReference) -> SamplerTrace:
    return replace(trace, refinement_parameters=(
        *(option for option in trace.refinement_parameters if option.name not in {_PATH_KEY, _SHA_KEY}),
        BackendOption(_PATH_KEY, reference.path), BackendOption(_SHA_KEY, reference.sha256),
    ))


def _digest(value: str) -> None:
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise DiffusionContractError("refinement evidence requires lowercase SHA-256 digests")


@dataclass(frozen=True)
class StagedRefinementEvidence:
    request_sha256: str
    raw: ArtifactReference
    selected: ArtifactReference
    selected_stage: Literal["raw", "generated-only", "generated-and-anchors"]
    stage1: ArtifactReference | None = None
    stage1_input_sha256: str | None = None
    stage2_input_sha256: str | None = None
    stage1_failure: NumericalFailureEvidence | None = None
    report: ArtifactReference | None = None
    revision: str = EVIDENCE_REVISION

    def __post_init__(self) -> None:
        if self.revision != EVIDENCE_REVISION:
            raise DiffusionContractError("unsupported staged refinement evidence revision")
        _digest(self.request_sha256)
        if not isinstance(self.raw, ArtifactReference) or not isinstance(self.selected, ArtifactReference):
            raise DiffusionContractError("staged refinement requires raw and selected artifacts")
        if self.stage1 is not None and not isinstance(self.stage1, ArtifactReference):
            raise DiffusionContractError("invalid stage-1 artifact")
        if self.report is not None and not isinstance(self.report, ArtifactReference):
            raise DiffusionContractError("invalid refinement report artifact")
        if self.stage1_failure is not None and not isinstance(self.stage1_failure, NumericalFailureEvidence):
            raise DiffusionContractError("invalid stage-1 numerical evidence")
        if self.stage1 is not None and self.stage1_failure is not None:
            raise DiffusionContractError("stage 1 cannot have coordinates and failure evidence")
        for digest in (self.stage1_input_sha256, self.stage2_input_sha256):
            if digest is not None:
                _digest(digest)
                if digest != self.raw.sha256:
                    raise DiffusionContractError("refinement stages must start from raw")
        if self.stage1_failure is not None and self.stage1_failure.input_sha256 != self.raw.sha256:
            raise DiffusionContractError("stage-1 numerical failure is not raw-bound")
        if self.selected_stage == "raw":
            if self.selected.sha256 != self.raw.sha256 or any(value is not None for value in (
                self.stage1, self.stage1_failure, self.stage1_input_sha256, self.stage2_input_sha256,
            )):
                raise DiffusionContractError("raw no-op evidence must not contain refinement attempts")
        elif self.selected_stage == "generated-only":
            if (self.stage1 is None or self.stage1.sha256 != self.selected.sha256
                    or self.stage1_input_sha256 is None or self.stage2_input_sha256 is not None):
                raise DiffusionContractError("generated-only evidence must select stage 1")
        elif self.selected_stage == "generated-and-anchors":
            if (self.stage1_input_sha256 is None or self.stage2_input_sha256 is None
                    or (self.stage1 is None and self.stage1_failure is None)):
                raise DiffusionContractError("anchor refinement requires raw-bound stage-1 evidence")
        else:
            raise DiffusionContractError("unknown selected refinement stage")

    def validate_request(self, request: DiffusionRequest, *, seed: int | None = None) -> None:
        digest = hashlib.sha256(request.to_json().encode("utf-8")).hexdigest()
        valid_digests = {digest}
        if seed is not None:
            if seed not in request.seeds:
                raise DiffusionContractError("staged candidate seed is absent from the request")
            # Production runners sample/refine each seed separately. Bind that
            # exact subrequest, never a runner-selected change to masks/policy.
            single = replace(request, candidate_count=1, seeds=(seed,))
            valid_digests.add(hashlib.sha256(single.to_json().encode("utf-8")).hexdigest())
        if self.request_sha256 not in valid_digests:
            raise DiffusionContractError("staged refinement request digest mismatch")
        policy = read_refinement_policy(request.backend_options)
        if policy is None:
            raise DiffusionContractError("historical requests cannot authorize staged refinement")
        if self.selected_stage == "generated-and-anchors" and not policy.relax_anchors:
            raise DiffusionContractError("request disables observed-anchor relaxation")

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"

    @classmethod
    def from_json(cls, text: str) -> StagedRefinementEvidence:
        if len(text.encode("utf-8")) > MAX_REFINEMENT_EVIDENCE_BYTES:
            raise DiffusionContractError("staged refinement evidence exceeds size limit")

        def unique(pairs):
            result = {}
            for key, value in pairs:
                if key in result:
                    raise DiffusionContractError("duplicate refinement evidence key")
                result[key] = value
            return result

        try:
            data = json.loads(text, object_pairs_hook=unique)
            if not isinstance(data, dict) or data.get("revision") != EVIDENCE_REVISION:
                raise DiffusionContractError("missing or unsupported refinement evidence revision")
            for key in ("raw", "selected", "stage1", "report"):
                if data.get(key) is not None:
                    data[key] = ArtifactReference(**data[key])
            if data.get("stage1_failure") is not None:
                data["stage1_failure"] = NumericalFailureEvidence(**data["stage1_failure"])
            return cls(**data)
        except (AttributeError, TypeError, ValueError) as exc:
            raise DiffusionContractError(f"invalid staged refinement evidence: {exc}") from exc
