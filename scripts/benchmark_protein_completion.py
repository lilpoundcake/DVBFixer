"""Private first-slice eight-mask control; not a broad benchmark or public backend.

Fixed declared source: tracked 1UBQ chain-A 1-72 protein-only projection. The
service sees the complete retained projection minus withheld atoms, never their
coordinates. Scoring consults withheld coordinates only after generation; it
does not choose candidates. All failed attempts and masks are retained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate

CASES = ((28, "ALA", ("CB",)), (46, "ALA", ("CB",)),
         (20, "SER", ("OG",)), (20, "SER", ("CB", "OG")),
         (57, "SER", ("OG",)), (57, "SER", ("CB", "OG")),
         (65, "SER", ("OG",)), (65, "SER", ("CB", "OG")))


def run_controls(root: Path, *, cases: tuple = CASES) -> dict:
    source = (root/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    lines = [line for line in source.decode().splitlines(keepends=True)
             if line.startswith("ATOM  ") and line[21] == "A" and int(line[22:26]) <= 72]
    projection = "".join(lines)+"TER\nEND\n"
    result = {"schema": "dvbfixer.private-protein-completion-control.v1",
              "source_sha256": hashlib.sha256(source).hexdigest(),
              "projection": "1UBQ chain A residues 1-72, protein only; not intact deposited complex",
              "projection_sha256": hashlib.sha256(projection.encode()).hexdigest(),
              "cases": {}, "empirical_rotamer_library": False,
              "scripts_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]}}
    for resid, name, mask in cases:
        target = ComponentInstanceRef(1, "A", str(resid))
        input_lines = [line for line in lines
                       if not (int(line[22:26]) == resid and line[12:16].strip() in mask)]
        input_bytes = ("".join(input_lines)+"TER\nEND\n").encode()
        observed = []
        for line in input_lines:
            component = ComponentInstanceRef(1, line[21], line[22:26].strip(), line[26].strip())
            observed.append(Observation(ExactAtomRef(component, line[12:16].strip()),
                                        line[76:78].strip(), tuple(float(line[i:i+8]) for i in (30, 38, 46)),
                                        float(line[54:60])))
        # This is a scoped fixture adapter, not a second production PDB parser or
        # peptide inference algorithm. All 1-72 peptide neighbours are declared;
        # the service independently rejects malformed target C-N link lengths.
        by_identity = {o.atom: o for o in observed}
        links = tuple((ExactAtomRef(ComponentInstanceRef(1, "A", str(i)), "C"),
                       ExactAtomRef(ComponentInstanceRef(1, "A", str(i+1)), "N")) for i in range(1, 72))
        if any(a not in by_identity or b not in by_identity for a, b in links):
            raise ValueError("fixture projection no longer contains declared peptide backbone")
        request = CompletionRequest(hashlib.sha256(input_bytes).hexdigest(), target, name,
                                    tuple(observed), mask, links)
        start = time.monotonic()
        completion = complete_sidechain(request)
        elapsed = time.monotonic()-start
        # Only now read the held-out coordinates for separate reporting.
        truth = {line[12:16].strip(): tuple(float(line[i:i+8]) for i in (30, 38, 46))
                 for line in lines if int(line[22:26]) == resid and line[12:16].strip() in mask}
        scores = []
        for attempt in completion.attempts:
            audit = audit_candidate(request, attempt.candidate)
            if audit != attempt.audit:
                raise ValueError("independent readback disagrees with returned geometry audit")
            distances = [math.dist(a.xyz, truth[a.atom.name]) for a in attempt.candidate.additions]
            scores.append({"heldout_rmsd_angstrom": math.sqrt(sum(d*d for d in distances)/len(distances)),
                           "maximum_error_angstrom": max(distances)})
        result["cases"][f"{name}{resid}-{'-'.join(mask)}"] = {
            "request": asdict(request), "result": asdict(completion), "scores_not_used_for_selection": scores,
            "elapsed_seconds": elapsed, "max_observed_movement_angstrom": 0.0,
            "preservation_basis": "generated-only result; original immutable observations never rewritten"}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    result = run_controls(Path(__file__).resolve().parents[1])
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    if args.summary:
        summary = {k: v for k, v in result.items() if k != "cases"}
        summary["raw_evidence_sha256"] = hashlib.sha256(args.output.read_bytes()).hexdigest()
        summary["cases"] = {}
        for name, case in result["cases"].items():
            completion = case["result"]
            selected = completion["selected_index"]
            summary["cases"][name] = {
                "status": completion["status"], "reason": completion["reason"],
                "request_digest": completion["request_digest"],
                "requested_missing": case["request"]["missing_atoms"],
                "observation_count": len(case["request"]["observations"]),
                "attempts": completion["attempts"], "selected_index": selected,
                "scores_not_used_for_selection": case["scores_not_used_for_selection"],
                "selected_heldout_rmsd_angstrom": None if selected is None else case["scores_not_used_for_selection"][selected]["heldout_rmsd_angstrom"],
                "max_observed_movement_angstrom": case["max_observed_movement_angstrom"],
                "preservation_basis": case["preservation_basis"],
                "elapsed_seconds": case["elapsed_seconds"],
            }
        with args.summary.open("x") as handle:
            json.dump(summary, handle, indent=2, allow_nan=False)
            handle.write("\n")
    print(json.dumps({k: {"status": v["result"]["status"],
                          "attempts": len(v["result"]["attempts"]),
                          "selected_index": v["result"]["selected_index"]}
                      for k, v in result["cases"].items()}, indent=2))


if __name__ == "__main__":
    main()
