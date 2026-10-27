"""Audit saved staged 9CQD outputs against EXPLICIT research movement bounds."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from dataclasses import asdict
from pathlib import Path

from dvbfixer.model.diffusion.anchor_acceptance import audit_staged_anchor_repair
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.refinement_stages import NumericalFailureEvidence


def checked_text(root: Path, attempt: dict) -> str:
    path = (root / attempt["coordinate_file"]).resolve(strict=True)
    if not path.is_relative_to(root.resolve()):
        raise ValueError("evidence path escapes pilot")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != attempt["coordinate_sha256"]:
        raise ValueError("coordinate evidence digest mismatch")
    return data.decode("ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--backbone-rmsd-max", required=True, type=float)
    parser.add_argument("--heavy-displacement-max", required=True, type=float)
    args = parser.parse_args()
    limits = AnchorMovementLimits(args.backbone_rmsd_max, args.heavy_displacement_max)
    request = DiffusionRequest.from_json((args.pilot / "request.json").read_text())
    source = (args.pilot / "source.pdb").read_bytes().decode("ascii")
    rows = json.loads((args.pilot / "comparison.json").read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    from dvbfixer.model.diffusion import anchor_acceptance
    shutil.copy2(anchor_acceptance.__file__, args.output / "anchor_acceptance.py")
    shutil.copy2(__file__, args.output / Path(__file__).name)
    manifest = {"movement_limits": asdict(limits), "publication_qualified": False,
        "purpose": "sensitivity audit with explicit research bounds, not scientific calibration",
        "pilot_report_sha256": hashlib.sha256((args.pilot / "comparison.json").read_bytes()).hexdigest(),
        "auditor_sha256": hashlib.sha256(Path(anchor_acceptance.__file__).read_bytes()).hexdigest()}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    report = []
    for row in rows:
        attempts = {attempt["stage"]: attempt for attempt in row["attempts"]}
        if set(attempts) != {"raw", "generated-only", "generated-and-anchors"}:
            raise ValueError("pilot requires all stages")
        raw = checked_text(args.pilot, attempts["raw"])
        first_attempt = attempts["generated-only"]
        first = checked_text(args.pilot, first_attempt) if "coordinate_file" in first_attempt else None
        failure = NumericalFailureEvidence(**first_attempt["numerical_failure"]) if first_attempt.get("numerical_failure") else None
        final = checked_text(args.pilot, attempts["generated-and-anchors"])
        audit = audit_staged_anchor_repair(request, source_text=source, raw_text=raw,
            stage1_text=first, candidate_text=final,
            stage1_input_sha256=attempts["generated-only"]["input_sha256"],
            stage2_input_sha256=attempts["generated-and-anchors"]["input_sha256"], limits=limits,
            stage1_failure=failure)
        result = {"arm": row["arm"], "bfactor_policy": row["bfactor_policy"], **asdict(audit)}
        report.append(result)
        (args.output / "audits.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
        print(row["arm"], row["bfactor_policy"], "AUDIT", audit.summary.passed,
              "FAILURES", audit.summary.hard_gate_failures, "PUBLICATION", audit.publication_qualified, flush=True)


if __name__ == "__main__":
    main()
