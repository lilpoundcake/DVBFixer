"""Commit and replay eight retained mixed controls; no new scientific masks/seeds."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
from dvbfixer.protein_completion.mixed_cohort import complete_mixed_cohort, plan_mixed_cohort
from dvbfixer.protein_completion.mixed_evidence import read_mixed_evidence, write_mixed_evidence

CASES = {"ubq-four-mixed", "ubq-two-thr", "fr2-paired-anchored", "a3n-whole-val-anchored-thr",
         "a3n-occupancy-refusal", "synthetic-mixed-choice", "synthetic-two-anchors-refusal", "synthetic-rounding-refusal"}
_FILES = {"input.pdb", "source-plan.json", "attempts.json", "completed.pdb", "completed.dat", "audit.json"}


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
    if args.retained.is_symlink():
        raise ValueError("retained controls must be a regular directory")
    retained = args.retained.resolve(strict=True)
    raw_record = read_private_file_bytes(retained, "evidence.json", max_bytes=2_000_000)
    original = json.loads(raw_record)
    if original.get("schema") != "dvbfixer.private-mixed-cohort-controls.v1" or set(original["results"]) != CASES:
        raise ValueError("complete retained eight-case mixed control record required")
    declarations = {row["id"]: row for row in original["cases"]}
    if set(declarations) != CASES or len(original["cases"]) != len(CASES):
        raise ValueError("exact eight source/target declarations required")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"plan.json").write_text(json.dumps({"retained_record_sha256": sha(raw_record), "cases": sorted(CASES),
        "purpose": "Commit/readback of retained deterministic mixed controls, not new scientific seeds"}, indent=2)+"\n")
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = {"schema": "dvbfixer.private-mixed-evidence-controls.v1", "retained_directory": str(args.retained),
              "retained_record_sha256": sha(raw_record), "output_directory": str(args.output),
              "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources},
              "public_profile_enabled": False, "scientific_qualification": False, "results": {}}
    for name in sorted(CASES):
        previous = original["results"][name]
        payloads = {}
        for filename, digest in previous["retained_files"].items():
            if filename not in _FILES:
                raise ValueError("unexpected retained file")
            data = read_private_file_bytes(retained/name, filename, max_bytes=64*1024*1024)
            if sha(data) != digest:
                raise ValueError("retained mixed artifact changed")
            payloads[filename] = data
        source = payloads["input.pdb"]
        if sha(source) != declarations[name]["source_sha256"]:
            raise ValueError("retained source differs from declaration")
        targets = tuple(ComponentInstanceRef(**t) for t in declarations[name]["targets"])
        plan = None
        result = None
        if previous["status"] != "admission_refused":
            plan = plan_mixed_cohort(source, targets)
            result = complete_mixed_cohort(source, plan)
            if (canonical(asdict(plan)) != canonical(json.loads(payloads["source-plan.json"]))
                    or canonical(asdict(result)) != canonical(json.loads(payloads["attempts.json"]))):
                raise ValueError("retained mixed plan/attempts differ from current deterministic replay")
        destination = write_mixed_evidence(source, targets, args.output/name, plan=plan, result=result)
        evidence = read_mixed_evidence(destination, expected_source_sha256=sha(source))
        expected_status = "accepted" if previous["status"] == "accepted" else "rejected"
        expected_stage = {"accepted": "rounded", "admission_refused": "admission",
                          "joint_rejected": "raw", "rounded_refused": "rounded"}[previous["status"]]
        if evidence.status != expected_status or evidence.stage != expected_stage:
            raise ValueError("evidence commit changed the retained outcome")
        if evidence.materialized is not None:
            if evidence.materialized.pdb != payloads["completed.pdb"]:
                raise ValueError("evidence commit changed accepted coordinates/records")
            dat = DatRecord.from_json(payloads["completed.dat"].decode())
            if dat.to_dict() != evidence.materialized.dat_record().to_dict():
                raise ValueError("evidence commit changed generated ownership")
        elif previous["status"] == "admission_refused" and evidence.refusal != previous["reason"]:
            raise ValueError("replayed source admission refusal differs")
        report["results"][name] = {"status": evidence.status, "stage": evidence.stage,
            "plan_digest": evidence.plan.digest if evidence.plan else None,
            "retained_outcome_unchanged": True, "independent_replay_passed": True,
            "bundle_files": {p.name: sha(p.read_bytes()) for p in sorted(destination.iterdir())}}
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({n: {"status": r["status"], "stage": r["stage"]} for n, r in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
