"""Artifact bindings do not imply scientific acceptance."""

import hashlib
import json
from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.contract import ArtifactReference, DiffusionContractError
from dvbfixer.model.diffusion.refinement_evidence import StagedRefinementEvidence
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
from dvbfixer.model.diffusion.refinement_stages import NumericalFailureEvidence
from tests.test_diffusion_contract import _request

RAW = ArtifactReference("evidence/raw.pdb", "a" * 64)
FIRST = ArtifactReference("evidence/stage1.pdb", "b" * 64)
FINAL = ArtifactReference("candidate.pdb", "c" * 64)


def evidence(*, relax_anchors=True):
    request = replace(_request(), backend_options=RefinementPolicy(relax_anchors=relax_anchors).backend_options())
    record = StagedRefinementEvidence(
        hashlib.sha256(request.to_json().encode()).hexdigest(), RAW, FINAL,
        "generated-and-anchors", FIRST, RAW.sha256, RAW.sha256,
    )
    return request, record


def test_roundtrip_and_request_binding():
    request, record = evidence()
    record.validate_request(request)
    assert StagedRefinementEvidence.from_json(record.to_json()) == record
    with pytest.raises(DiffusionContractError, match="digest mismatch"):
        record.validate_request(replace(request, seeds=tuple(seed + 1 for seed in request.seeds)))


def test_no_anchor_authority_for_disabled_or_historical_requests():
    request, record = evidence(relax_anchors=False)
    with pytest.raises(DiffusionContractError, match="disables"):
        record.validate_request(request)
    historical = _request()
    record = replace(record, request_sha256=hashlib.sha256(historical.to_json().encode()).hexdigest())
    with pytest.raises(DiffusionContractError, match="historical"):
        record.validate_request(historical)


def test_numerical_failure_can_replace_but_not_accompany_stage_one_coordinates():
    _, record = evidence()
    failure = NumericalFailureEvidence(RAW.sha256, "BoundaryNumericalError", "No finite minimum")
    recovered = replace(record, stage1=None, stage1_failure=failure)
    assert StagedRefinementEvidence.from_json(recovered.to_json()) == recovered
    with pytest.raises(DiffusionContractError, match="coordinates and failure"):
        replace(record, stage1_failure=failure)


@pytest.mark.parametrize("mutation", [
    {"stage2_input_sha256": FIRST.sha256}, {"stage1_input_sha256": FINAL.sha256},
    {"stage1": None}, {"selected_stage": "future"}, {"revision": "future"},
    {"stage2_input_sha256": None},
])
def test_incomplete_or_sequential_restart_evidence_is_rejected(mutation):
    _, record = evidence()
    with pytest.raises(DiffusionContractError):
        replace(record, **mutation)


def test_raw_noop_and_generated_only_custody():
    request, record = evidence()
    for selected in (
        StagedRefinementEvidence(record.request_sha256, RAW, RAW, "raw"),
        StagedRefinementEvidence(record.request_sha256, RAW, FIRST, "generated-only", FIRST, RAW.sha256),
    ):
        selected.validate_request(request)
        assert StagedRefinementEvidence.from_json(selected.to_json()) == selected
    with pytest.raises(DiffusionContractError):
        StagedRefinementEvidence(record.request_sha256, RAW, FINAL, "raw")


@pytest.mark.parametrize("text", ["[]", "null", "{}", '{"raw":{},"raw":{}}', "x" * 16385])
def test_bounded_codec_rejects_malformed_evidence(text):
    with pytest.raises(DiffusionContractError):
        StagedRefinementEvidence.from_json(text)


def test_artifact_traversal_and_unknown_fields_are_rejected():
    _, record = evidence()
    data = json.loads(record.to_json())
    data["raw"]["path"] = "../source.pdb"
    with pytest.raises(DiffusionContractError):
        StagedRefinementEvidence.from_json(json.dumps(data))
    data = json.loads(record.to_json()) | {"authorized_anchors": []}
    with pytest.raises(DiffusionContractError):
        StagedRefinementEvidence.from_json(json.dumps(data))


def test_missing_revision_and_misbound_numerical_failure_are_rejected():
    _, record = evidence()
    data = json.loads(record.to_json())
    data.pop("revision")
    with pytest.raises(DiffusionContractError, match="revision"):
        StagedRefinementEvidence.from_json(json.dumps(data))
    failure = NumericalFailureEvidence(FIRST.sha256, "BoundaryNumericalError", "No finite minimum")
    with pytest.raises(DiffusionContractError, match="not raw-bound"):
        replace(record, stage1=None, stage1_failure=failure)


def test_seed_scoped_request_digest_cannot_authorize_different_masks_or_seeds():
    parent, record = evidence()
    seed = parent.seeds[0]
    single = replace(parent, candidate_count=1, seeds=(seed,))
    scoped = replace(record, request_sha256=hashlib.sha256(single.to_json().encode()).hexdigest())
    scoped.validate_request(parent, seed=seed)
    with pytest.raises(DiffusionContractError, match="seed is absent"):
        scoped.validate_request(parent, seed=max(parent.seeds) + 1)
    changed = replace(parent, backend_options=RefinementPolicy(bfactor_weighting=False).backend_options())
    with pytest.raises(DiffusionContractError, match="digest mismatch"):
        scoped.validate_request(changed, seed=seed)
