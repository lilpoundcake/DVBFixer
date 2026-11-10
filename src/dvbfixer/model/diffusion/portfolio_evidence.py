"""Independent bounded custody and geometry audit of a multi-parent artifact."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import (
    DiffusionContractError,
    DiffusionRequest,
    RunnerResult,
    SamplerTrace,
)
from dvbfixer.model.diffusion.gap_portfolio import (
    PORTFOLIO_KIND,
    PORTFOLIO_REVISION,
    GapParent,
    _coordinates,
    audit_portfolio_materialization,
    digest,
    plan_gap_repairs,
)
from dvbfixer.model.diffusion.input_preparation import (
    audit_preparation_payloads,
    preparation_options,
)
from dvbfixer.model.diffusion.portfolio_contract import PortfolioEvidence
from dvbfixer.model.diffusion.private_artifact import read_private_artifact_bytes
from dvbfixer.model.diffusion.refinement_evidence import (
    StagedRefinementEvidence,
    evidence_reference,
)
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter
from dvbfixer.model.diffusion.trace import validate_sampler_trace_context

MAX_PORTFOLIO_PAYLOAD_BYTES = 16_000_000
MAX_PORTFOLIO_TOTAL_BYTES = 128_000_000


def json_bytes(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode("utf-8")


def strict_json(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise DiffusionContractError("duplicate key in portfolio JSON")
            result[key] = value
        return result

    def nonfinite(value):
        raise DiffusionContractError("nonfinite value in portfolio JSON")

    return json.loads(data, object_pairs_hook=unique, parse_constant=nonfinite)


def portfolio_resume_key(evidence: PortfolioEvidence) -> str:
    """Paths are relocatable; bytes, parent identities, choices and sources are not."""
    return hashlib.sha256(json_bytes({"revision": 1,
        "artifacts": {role: a.sha256 for role, a in evidence.artifacts().items()},
        "parents": [(p.candidate_id, p.seed) for p in evidence.parents],
        "choices": [asdict(c) for c in evidence.choices]})).hexdigest()


@dataclass(frozen=True)
class PortfolioEvidenceAudit:
    request: DiffusionRequest
    payloads: dict[str, bytes]
    selected_assessment: RepairAssessment
    final_assessment: RepairAssessment
    preparation_verified: bool

    def validation_record(self):
        return {"kind": "independent-private-portfolio-validation", "schema_version": 1,
            "request_sha256": digest(self.request.to_json()),
            "source_sha256": hashlib.sha256(self.payloads["source"]).hexdigest(),
            "candidate_sha256": hashlib.sha256(self.payloads["candidate"]).hexdigest(),
            "selected_assessment": asdict(self.selected_assessment),
            "final_assessment": asdict(self.final_assessment),
            "preparation_verified": self.preparation_verified,
            "passed": self.final_assessment.geometry_passed and not self.final_assessment.admission_failures,
            "scientific_qualification": False, "public_profile_enabled": False}


def audit_portfolio_evidence(evidence: PortfolioEvidence, workspace: Path) -> PortfolioEvidenceAudit:
    if workspace.is_symlink() or not workspace.is_dir():
        raise DiffusionContractError("portfolio workspace must be a regular directory")
    root = workspace.resolve(strict=True)
    payloads = {}
    total = 0
    for role, artifact in evidence.artifacts().items():
        limit = (256_000 if role in {"selection_source", "repair_selector_source", "repair_kernel_source",
                                    "repair_conditioning_source"} else MAX_PORTFOLIO_PAYLOAD_BYTES)
        data = read_private_artifact_bytes(root, artifact, max_bytes=min(limit, MAX_PORTFOLIO_TOTAL_BYTES - total))
        payloads[role] = data
        total += len(data)
    request = DiffusionRequest.from_dict(strict_json(payloads["request"]))
    source = payloads["source"].decode("utf-8")
    if hashlib.sha256(payloads["source"]).hexdigest() != request.normalized_pdb.sha256:
        raise DiffusionContractError("portfolio normalized source differs from original request")
    prepared = preparation_options(request) is not None
    if prepared != (evidence.preparation is not None):
        raise DiffusionContractError("portfolio must retain exactly the request's source preparation")
    if prepared:
        preparation = {role: payloads[role] for role in evidence.preparation.artifacts()}
        audit_preparation_payloads(request, preparation)
        if preparation["input-prepared"] != payloads["source"]:
            raise DiffusionContractError("portfolio prepared bytes differ from normalized source")
    parents = []
    base_generated_records = None
    for i, entry in enumerate(evidence.parents):
        data = {name: payloads[f"parent-{i:02d}-{name}"] for name in entry.artifacts()}
        trace = SamplerTrace.from_dict(strict_json(data["sampler_trace"]))
        staged = StagedRefinementEvidence.from_json(data["staged_evidence"].decode("utf-8"))
        staged.validate_request(request, seed=entry.seed)
        bound = evidence_reference(trace)
        result = RunnerResult.from_dict(strict_json(data["runner_result"]))
        if len(result.candidates) != 1:
            raise DiffusionContractError("portfolio parent result must identify one original sampling candidate")
        candidate = result.candidates[0]
        if (candidate.seed != entry.seed or candidate.candidate_id != entry.candidate_id
                or tuple(candidate.generated_atoms) != tuple(request.generated_atoms)):
            raise DiffusionContractError("portfolio parent identity/seed/generated mask differs")
        validate_sampler_trace_context(trace, request, result.backend_provenance)
        if (bound is None or bound.sha256 != entry.staged_evidence.sha256 or staged.stage1 is None
                or staged.stage1.sha256 != entry.coordinates.sha256
                or staged.raw.sha256 != entry.raw.sha256 or candidate.coordinate_artifact.sha256 != entry.raw.sha256
                or staged.stage1_input_sha256 != entry.raw.sha256):
            raise DiffusionContractError("portfolio trace/raw/stage-1 custody differs")
        # Non-geometric donor failures cannot become eligible by selecting a gap.
        raw_text = data["raw"].decode("utf-8")
        raw_assessment = assessment_adapter(request, source)(raw_text)
        if raw_assessment.admission_failures:
            raise DiffusionContractError("portfolio parent raw proposal has admission failures")
        stage_text = data["coordinates"].decode("utf-8")
        coordinates, records = _coordinates(stage_text)
        if not set(request.generated_atoms) <= coordinates.keys():
            raise DiffusionContractError("portfolio parent stage-1 omits generated atoms")
        projection = rewrite_generated_coordinates(raw_text, {a: coordinates[a] for a in request.generated_atoms})
        if projection != stage_text:
            raise DiffusionContractError("portfolio parent stage-1 differs from a raw coordinate-only projection")
        structure = _structure(stage_text, "portfolio-donor-stage1.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        if base_generated_records is not None:
            for atom in request.generated_atoms:
                first, current = base_generated_records[atom], records[atom]
                if first[11:30] != current[11:30] or first[54:] != current[54:]:
                    raise DiffusionContractError("portfolio donor generated metadata differs")
        else:
            base_generated_records = records
        parents.append(GapParent(entry.candidate_id, entry.seed, digest(request.to_json()), digest(source),
            entry.coordinates.sha256, entry.sampler_trace.sha256, entry.staged_evidence.sha256,
            stage_text))
    selected, final = (payloads[role].decode("utf-8") for role in ("selected", "candidate"))
    if len(evidence.choices) != len(request.gaps):
        raise DiffusionContractError("portfolio lineage differs from request whole-gap inventory")
    protocol = strict_json(payloads["search_protocol"])
    report = strict_json(payloads["search_report"])
    if (protocol.get("kind") != PORTFOLIO_KIND or protocol.get("stage") != "generated-only"
            or protocol.get("anchor_mobility") is not False
            or protocol.get("seeds") != [p.seed for p in evidence.parents]
            or type(protocol.get("max_evaluations")) is not int
            or not len(parents) <= protocol["max_evaluations"] <= 4096
            or type(protocol.get("beam_width")) is not int or not 1 <= protocol["beam_width"] <= 32):
        raise DiffusionContractError("portfolio search protocol/budget differs")
    if (report.get("kind") != PORTFOLIO_KIND or report.get("revision") != PORTFOLIO_REVISION
            or report.get("publication_qualified") is not False
            or report.get("selected", {}).get("choices") != [asdict(c) for c in evidence.choices]
            or report.get("selected", {}).get("coordinate_sha256") != evidence.selected.sha256
            or not isinstance(report.get("trials"), list)
            or not len(parents) <= len(report["trials"]) <= protocol["max_evaluations"]):
        raise DiffusionContractError("portfolio selection/report differs from bound lineage")
    before = audit_portfolio_materialization(request, source, tuple(parents), evidence.choices, selected)
    if evidence.repair_report is None:
        if payloads["selected"] != payloads["candidate"]:
            raise DiffusionContractError("portfolio changed coordinates without an explicit repair operation")
    else:
        operation = strict_json(payloads["repair_report"])
        scope = plan_gap_repairs(request, before)
        residues = {r for i in scope.gap_indices for r in request.gaps[i].generated_residues}
        mobile = [asdict(a) for a in request.generated_atoms
            if any((a.chain, a.residue_number, a.insertion_code) ==
                   (r.chain, r.residue_number, r.insertion_code) for r in residues)]
        if (operation.get("operation") != "generated-only-issue-local-refinement-v1"
                or operation.get("publication_qualified") is not False
                or operation.get("input_sha256") != evidence.selected.sha256
                or operation.get("output_sha256") != evidence.candidate.sha256
                or json_bytes(operation.get("scope")) != json_bytes(asdict(scope))
                or operation.get("mobile_atoms") != mobile
                or type(operation.get("refinement_seed")) is not int or operation["refinement_seed"] < 0
                or type(operation.get("max_iterations")) is not int or not 1 <= operation["max_iterations"] <= 2000
                or type(operation.get("restart_count")) is not int or not 0 <= operation["restart_count"] <= 8):
            raise DiffusionContractError("portfolio repair operation differs from independent mobility/budget")
    after = audit_portfolio_materialization(request, source, tuple(parents), evidence.choices,
        selected, repaired_text=final if evidence.repair_report is not None else None)
    if not after.geometry_passed or after.admission_failures:
        raise DiffusionContractError("portfolio final whole-complex geometry/admission failed")
    return PortfolioEvidenceAudit(request, payloads, before, after, prepared)
