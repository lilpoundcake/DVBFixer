"""Reporting for every seed, separate from the accepted-candidate contract."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict
from typing import TYPE_CHECKING

from dvbfixer.model.diffusion.contract import (
    ArtifactReference,
    DiffusionCandidate,
    DiffusionRequest,
    DiffusionResult,
    RunnerCandidate,
    ValidationSummary,
)

if TYPE_CHECKING:
    from dvbfixer.model.diffusion.validate import CandidateValidation

MAX_SEED_OUTCOME_BYTES = 10_000_000
SEED_OUTCOMES_FILE = "seed-outcomes.json"


def stage_seed_outcomes(
    request: DiffusionRequest,
    result: DiffusionResult,
    validations: Sequence[CandidateValidation],
    *,
    read_coordinates: Callable[[ArtifactReference], bytes],
    write_artifact: Callable[[str, bytes], None],
) -> ArtifactReference:
    """Keep final summaries and rejected coordinates without changing selection.

    Missing runner candidates are explicitly *not returned*, not geometry-rejected.
    Rejected coordinates are inspection-only, outside the bundle's candidate list.
    """
    by_seed: dict[int, tuple[RunnerCandidate | DiffusionCandidate, ValidationSummary]] = {
        item.candidate.seed: (item.candidate, item.summary) for item in validations}
    if len(by_seed) != len(validations) or set(by_seed) - set(request.seeds):
        raise ValueError("seed outcome validations repeat or introduce a seed")
    # Older callers/mocked validators may supply only the selected result.
    for selected_candidate, selected_summary in zip(result.candidates, result.validation_summaries, strict=True):
        if selected_candidate.seed in by_seed:
            observed, final_summary = by_seed[selected_candidate.seed]
            if observed.candidate_id != selected_candidate.candidate_id or final_summary != selected_summary:
                raise ValueError("seed report differs from the selected validation result")
        by_seed.setdefault(selected_candidate.seed, (selected_candidate, selected_summary))
    published = {candidate.seed: candidate for candidate in result.candidates}
    if len(published) != len(result.candidates) or set(published) - set(request.seeds):
        raise ValueError("selected seed inventory disagrees with the request")
    rows: list[dict[str, object]] = []
    for seed in request.seeds:
        item = by_seed.get(seed)
        if item is None:
            rows.append({"seed": seed, "candidate_id": None, "status": "not_returned",
                         "validation": None, "coordinate": None})
            continue
        candidate, summary = item
        if seed in published:
            # Published PDB coordinates undergo metadata reconciliation; the
            # writer supplies their digest in the bundle index. No
            # duplicate coordinate reference is asserted for these rows.
            coordinate_dict = None
        else:
            if summary.passed:
                raise ValueError("passing seed is absent from the selected result")
            data = read_coordinates(candidate.coordinate_artifact)
            name = f"rejected/seed-{seed}.pdb"
            write_artifact(name, data)
            coordinate = ArtifactReference(name, hashlib.sha256(data).hexdigest())
            coordinate_dict = asdict(coordinate)
        rows.append({"seed": seed, "candidate_id": candidate.candidate_id,
                     "status": "accepted" if summary.passed else "rejected",
                     "validation": asdict(summary), "coordinate": coordinate_dict})
    data = (json.dumps({"version": 1, "request_sha256": hashlib.sha256(request.to_json().encode()).hexdigest(),
                        "requested_seeds": list(request.seeds),
                        "outcomes": rows}, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode()
    if len(data) > MAX_SEED_OUTCOME_BYTES:
        raise ValueError("seed outcome report exceeds the size limit")
    write_artifact(SEED_OUTCOMES_FILE, data)
    return ArtifactReference(SEED_OUTCOMES_FILE, hashlib.sha256(data).hexdigest())


def verify_seed_outcomes(
    data: bytes, *, request_sha256: str, published: Mapping[int, tuple[str, bool]],
    read_artifact: Callable[[ArtifactReference], bytes],
) -> dict:
    """Verify reporting custody/consistency, never grant scientific acceptance."""
    if len(data) > MAX_SEED_OUTCOME_BYTES:
        raise ValueError("seed outcome report exceeds the size limit")
    report = json.loads(data)
    if (not isinstance(report, dict) or type(report.get("version")) is not int
            or report["version"] != 1 or report.get("request_sha256") != request_sha256):
        raise ValueError("invalid seed outcome report binding")
    seeds, rows = report.get("requested_seeds"), report.get("outcomes")
    if (not isinstance(seeds, list) or not 1 <= len(seeds) <= 1000
            or any(type(seed) is not int or not 0 <= seed <= 2**53 - 1 for seed in seeds)
            or len(set(seeds)) != len(seeds) or not isinstance(rows, list) or len(rows) != len(seeds)):
        raise ValueError("invalid seed outcome inventory")
    for seed, row in zip(seeds, rows, strict=True):
        if not isinstance(row, dict) or type(row.get("seed")) is not int or row["seed"] != seed:
            raise ValueError("seed outcome order or identity mismatch")
        status = row.get("status")
        if status == "not_returned":
            if seed in published or any(row.get(key) is not None for key in ("validation", "coordinate", "candidate_id")):
                raise ValueError("not-returned seed claims a candidate")
            continue
        summary = row.get("validation")
        if (status not in {"accepted", "rejected"} or not isinstance(row.get("candidate_id"), str)
                or not row["candidate_id"] or not isinstance(summary, dict)
                or summary.get("passed") is not (status == "accepted")
                or not isinstance(summary.get("hard_gate_failures"), list)
                or not isinstance(summary.get("warnings"), list)
                or any(not isinstance(text, str) for text in (*summary["hard_gate_failures"], *summary["warnings"]))
                or bool(summary["hard_gate_failures"]) != (status == "rejected")):
            raise ValueError("inconsistent seed validation outcome")
        if seed in published:
            if published[seed] != (row["candidate_id"], status == "accepted") or row.get("coordinate") is not None:
                raise ValueError("seed outcome differs from published candidate")
        elif status != "rejected" or row.get("coordinate") is None:
            raise ValueError("unpublished seed must retain rejected coordinates")
        else:
            reference = ArtifactReference(**row["coordinate"])
            if reference.path != f"rejected/seed-{seed}.pdb":
                raise ValueError("rejected coordinate has an unexpected path")
            read_artifact(reference)
    if set(published) - set(seeds):
        raise ValueError("seed report omits a published candidate")
    return report


def seed_outcome_lines(report: dict, bundle: str) -> tuple[str, ...]:
    """Render final-only Python-owned findings, including mixed-success runs."""
    lines = []
    for row in report["outcomes"]:
        seed, status = row["seed"], row["status"]
        if status == "accepted":
            lines.append(f"Diffusion seed {seed}: ACCEPTED ({row['candidate_id']}).")
        elif status == "not_returned":
            lines.append(f"WARNING: Diffusion seed {seed}: NOT RETURNED for independent final validation; "
                         "no scientific acceptance or geometry rejection is inferred.")
        else:
            lines.append(f"ERROR: Diffusion candidate REJECTED: {row['candidate_id']} "
                         f"(seed {seed}; final selected coordinates).")
            summary = row["validation"]
            findings = tuple(dict.fromkeys(warning.removeprefix("Final validation: ")
                for warning in summary["warnings"] if warning.startswith("Final validation: ")))
            lines.extend(f"ERROR: {finding}" for finding in findings)
            lines.append("ERROR: Failed gates: " + ", ".join(summary["hard_gate_failures"]))
            if not findings:
                lines.append("No residue-level detail is available for these failures; inspect the validation report.")
            lines.append("Rejected coordinates are saved for inspection only, not accepted for use.")
    lines.append(f"Per-seed validation report: {bundle}/{SEED_OUTCOMES_FILE}")
    return tuple(lines)
