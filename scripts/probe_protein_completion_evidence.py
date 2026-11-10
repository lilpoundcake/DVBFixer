"""Publish/read back the six retained deterministic cohort controls, without new masks.

Original artifacts are checksum-verified and never edited. Replay is evidence
validation of an unchanged deterministic control, not new scientific seeds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
from dvbfixer.protein_completion.cohort import complete_cohort, plan_cohort
from dvbfixer.protein_completion.evidence import (
    _targets,
    read_cohort_evidence,
    write_cohort_evidence,
)

CASES = ("ubq-two-ala", "ubq-three-ser", "ubq-four-mixed", "synthetic-ser-joint-options",
         "synthetic-ala-no-joint-solution", "synthetic-rounding-refusal")
_FILES = {"input.pdb", "request.json", "attempts.json", "completed.pdb", "completed.dat", "refusal.json"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    retained = args.retained.resolve(strict=True)
    raw_record = read_private_file_bytes(retained, "evidence.json", max_bytes=2_000_000)
    original = json.loads(raw_record)
    if set(original["results"]) != set(CASES) or original["all_expected_outcomes_matched"] is not True:
        raise ValueError("expected the complete retained six-case control record")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"plan.json").write_text(json.dumps({"retained_record_sha256": sha(raw_record),
        "cases": list(CASES), "purpose": "Retained deterministic control publication/readback, not new scientific seeds"}, indent=2)+"\n")
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = {"schema": "dvbfixer.private-protein-cohort-evidence-controls.v1", "retained_directory": str(args.retained),
              "retained_record_sha256": sha(raw_record), "output_directory": str(args.output),
              "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources},
              "public_profile_enabled": False, "scientific_qualification": False, "results": {}}
    for name in CASES:
        previous = original["results"][name]
        folder = retained/name
        payloads = {}
        for filename, digest in previous["retained_files"].items():
            if filename not in _FILES:
                raise ValueError("unexpected retained payload name")
            data = read_private_file_bytes(folder, filename, max_bytes=64*1024*1024)
            if sha(data) != digest:
                raise ValueError("retained control artifact changed")
            payloads[filename] = data
        source = payloads["input.pdb"]
        saved_plan = json.loads(payloads["request.json"])
        targets = _targets(canonical([r["residue"] for r in saved_plan["requests"]]))
        plan = plan_cohort(source, targets)
        if canonical(asdict(plan)) != canonical(saved_plan):
            raise ValueError("retained request differs from source replay")
        result = complete_cohort(source, plan)
        if canonical(asdict(result)) != canonical(json.loads(payloads["attempts.json"])):
            raise ValueError("retained proposals/attempts differ from deterministic replay")
        destination = write_cohort_evidence(source, plan, result, args.output/name)
        evidence = read_cohort_evidence(destination, expected_source_sha256=sha(source))
        expected_status = "accepted" if previous["outcome"] == "completed" else "rejected"
        if evidence.status != expected_status:
            raise ValueError("evidence publication changed retained scientific outcome")
        if evidence.materialized is not None:
            if evidence.materialized.pdb != payloads["completed.pdb"]:
                raise ValueError("publication changed accepted coordinates/records")
            dat = DatRecord.from_json(payloads["completed.dat"].decode())
            if dat.to_dict() != evidence.materialized.dat_record().to_dict():
                raise ValueError("publication changed generated-only sidecar")
        report["results"][name] = {"status": evidence.status, "stage": evidence.stage,
                                   "plan_digest": plan.digest, "independent_readback_passed": True,
                                   "retained_attempts_identical": True,
                                   "bundle_files": {p.name: sha(p.read_bytes()) for p in sorted(destination.iterdir())}}
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: v["status"] for k, v in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
