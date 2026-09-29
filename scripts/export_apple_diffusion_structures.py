#!/usr/bin/env python3
"""Export Apple diffusion input/output PDB pairs into case-ID directories."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_summary(directory: Path) -> dict[str, Any] | None:
    summary_path = directory / "summary.json"
    candidate_path = directory / "candidate.pdb"
    if not summary_path.is_file() or not candidate_path.is_file():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("candidate_sha256") != _sha256(candidate_path):
        raise ValueError(f"candidate digest mismatch: {candidate_path}")
    return summary


def _load_refinement_failure(directory: Path) -> dict[str, Any] | None:
    path = directory / "failure.json"
    if not path.is_file():
        return None
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") != "scientific-failure":
        raise ValueError(f"invalid refinement failure record: {path}")
    return value


def _atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copyfile(source, temporary)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def export_structures(
    manifest_path: Path,
    cohort_root: Path,
    output_root: Path,
    *,
    require_complete: bool = False,
) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = manifest.get("cases", [])
    prefix_count = int(manifest.get("reused_prefix", {}).get("case_count", 0))
    prefix_raw = str(manifest.get("reused_prefix", {}).get("raw_output_subdir", ""))
    prefix_refined = str(
        manifest.get("reused_prefix", {}).get("refined_output_subdir", "")
    )
    remaining_raw = str(manifest["sampling"]["remaining_output_subdir"])
    remaining_refined = str(manifest["refinement"]["remaining_output_subdir"])
    exported: list[dict[str, Any]] = []
    pending: list[str] = []

    for index, case in enumerate(cases):
        case_id = str(case["case_id"])
        if not case_id or Path(case_id).name != case_id:
            raise ValueError(f"unsafe case ID: {case_id!r}")
        workspace = cohort_root / "cases" / case_id / "workspace"
        if index < prefix_count:
            raw_dir = workspace / prefix_raw
            refined_dir = workspace / prefix_refined
        else:
            raw_dir = workspace / remaining_raw
            refined_dir = workspace / remaining_refined
        raw = _load_summary(raw_dir)
        refined = _load_summary(refined_dir)
        refinement_failure = _load_refinement_failure(refined_dir)
        if raw is None or (refined is None and refinement_failure is None):
            pending.append(case_id)
            continue

        if refined is not None and bool(refined["validation"]["passed"]):
            source = refined_dir / "candidate.pdb"
            selection = "refined"
            selected = refined
        elif bool(raw["validation"]["passed"]):
            source = raw_dir / "candidate.pdb"
            selection = "raw"
            selected = raw
        elif refined is not None:
            source = refined_dir / "candidate.pdb"
            selection = "failed-refined"
            selected = refined
        else:
            source = raw_dir / "candidate.pdb"
            selection = "failed-raw"
            selected = raw

        input_path = workspace / "input" / "normalized.pdb"
        if not input_path.is_file():
            raise FileNotFoundError(f"missing normalized input: {input_path}")
        reference_path = workspace / "reference.pdb"
        if not reference_path.is_file():
            raise FileNotFoundError(f"missing reference structure: {reference_path}")
        case_dir = output_root / case_id
        exported_input = case_dir / f"{case_id}_input.pdb"
        exported_reference = case_dir / f"{case_id}_reference.pdb"
        exported_output = case_dir / f"{case_id}_output.pdb"
        _atomic_copy(input_path, exported_input)
        _atomic_copy(reference_path, exported_reference)
        _atomic_copy(source, exported_output)
        metadata = {
            "schema_version": 1,
            "case_id": case_id,
            "screening_index": case["screening_index"],
            "input": {
                "file": exported_input.name,
                "sha256": _sha256(exported_input),
                "source": input_path.relative_to(cohort_root).as_posix(),
            },
            "reference": {
                "file": exported_reference.name,
                "sha256": _sha256(exported_reference),
                "source": reference_path.relative_to(cohort_root).as_posix(),
            },
            "output": {
                "file": exported_output.name,
                "sha256": _sha256(exported_output),
                "source": source.relative_to(cohort_root).as_posix(),
                "selection": selection,
                "validation_passed": bool(selected["validation"]["passed"]),
                "hard_gate_failures": selected["validation"]["hard_gate_failures"],
            },
        }
        if refinement_failure is not None:
            metadata["refinement_failure"] = refinement_failure
        _atomic_json(case_dir / f"{case_id}_metadata.json", metadata)
        exported.append(metadata)

    result = {
        "schema_version": 1,
        "manifest": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "case_count": len(cases),
        "exported_count": len(exported),
        "pending_count": len(pending),
        "pending_case_ids": pending,
        "cases": exported,
    }
    _atomic_json(output_root / "export-manifest.json", result)
    if require_complete and pending:
        raise RuntimeError(f"{len(pending)} case(s) are not complete")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    result = export_structures(
        args.manifest.resolve(),
        args.cohort_root.resolve(),
        args.output_root.resolve(),
        require_complete=args.require_complete,
    )
    print(
        f"exported={result['exported_count']} pending={result['pending_count']} "
        f"root={args.output_root}"
    )


if __name__ == "__main__":
    main()
