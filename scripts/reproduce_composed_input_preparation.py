#!/usr/bin/env python3
"""Retain an actual new-input completion + observed-geometry repair control."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path


def main():
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.model.diffusion.input_geometry import assess_input_geometry
    from dvbfixer.model.diffusion.input_preparation import audit_preparation_payloads, prepare_input
    from dvbfixer.model.diffusion.protein_context import prepare_cc94_request
    from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy
    from dvbfixer.model.diffusion.scope import assess_diffusion_scope
    from dvbfixer.model.diffusion_cli import _target_sequences, build_cli_diffusion_request

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(exist_ok=False)
    source = args.source.read_bytes()
    path = args.work / "source.pdb"
    path.write_bytes(source)
    repo = Path(__file__).resolve().parents[1]
    snapshots = {}
    for name in ("input_preparation.py", "observed_completion.py", "boundary_refinement.py"):
        data = (repo / "src/dvbfixer/model/diffusion" / name).read_bytes()
        (args.work / name).write_bytes(data)
        snapshots[name] = hashlib.sha256(data).hexdigest()
    (args.work / "source-code.json").write_text(json.dumps(snapshots, indent=2) + "\n")
    request = build_cli_diffusion_request(path, _target_sequences(path, None), seeds=(7, 101, 102, 103, 104),
        profile="protpardelle-1c-cc94-mps-v1", refinement_options=RefinementPolicy().backend_options())
    (args.work / "request-before.json").write_text(request.to_json())
    result = {"passed": False, "scope": "input preparation, NOT diffusion sampling"}
    try:
        before = assess_input_geometry(request, source)
        prepared, output = prepare_input(request, source, args.work)
        prepared = prepare_cc94_request(prepared, output.read_bytes(), "auto")
        payloads = {"input-original": source, "input-prepared": output.read_bytes(),
            "input-repair-report": (output.parent / "input-repair.json").read_bytes(),
            "input-repair-request": prepared.to_json().encode()}
        audit_preparation_payloads(prepared, payloads)
        scope = assess_diffusion_scope(prepared, output.read_bytes())
        after = assess_input_geometry(prepared, output.read_bytes())
        report = json.loads(payloads["input-repair-report"])
        result.update(passed=scope.supported, admission=asdict(scope), prepared_path=str(output),
            source_preserved=args.source.read_bytes() == source,
            original_severe_findings=sum(f.severity is Severity.ERROR for f in before),
            final_severe_findings=sum(f.severity is Severity.ERROR for f in after),
            final_warnings=sum(f.severity is Severity.WARNING for f in after),
            movements=report["movements"], revision=report["revision"],
            gap_ownership_preserved=prepared.gaps == request.gaps and prepared.generated_atoms == request.generated_atoms)
        (args.work / "request-final.json").write_text(prepared.to_json())
    except Exception as exc:
        result.update(error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        (args.work / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
