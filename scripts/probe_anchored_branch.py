"""Predeclared typed single-branch controls using fixed observed partners.

No PDB materialization, empirical rotamer statistics or public capability. The
complete declared projection minus one atom is passed to the typed kernel.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from dataclasses import asdict
from pathlib import Path

from branched_stereo_oracle import rdkit_labels
from rdkit import rdBase

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.pdbutils.sequence_evidence import observed_peptide_continuity
from dvbfixer.protein_completion.anchored_branch import (
    ANCHORED_POLICY,
    ANCHORED_REVISION,
    audit_anchored_branch,
    complete_anchored_branch,
)
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.source import _layout

FIXTURES = {"1UBQ": "tests/fixtures/protein_context/1UBQ.pdb",
            "1FR2": "tests/fixtures/regressions/1FR2.pdb", "1A3N": "tests/fixtures/protein_context/1A3N.pdb"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    originals = {name: (root/path).read_bytes() for name, path in FIXTURES.items()}
    for name, data in originals.items():
        (args.output/f"original-{name}.pdb").write_bytes(data)
    cases = []
    for fixture, stop, residues in (("1UBQ", 72, (("VAL", 5), ("VAL", 26), ("THR", 7), ("THR", 55))),
                                     ("1FR2", None, (("VAL", 19), ("THR", 11))),
                                     ("1A3N", 59, (("VAL", 1), ("THR", 8)))):
        for name, number in residues:
            branch = BRANCHES[name]
            for missing in (branch.first, branch.second):
                cases.append({"id": f"{fixture}-{name}{number}-{missing}", "fixture": fixture,
                              "stop": stop, "residue_name": name, "number": str(number), "missing": missing})
    plan = {"schema": "dvbfixer.private-observed-branch-controls.v1", "revision": ANCHORED_REVISION,
            "policy": ANCHORED_POLICY, "cases": cases, "fixture_sha256": {n: sha(b) for n, b in originals.items()},
            "projection": "Unmodified ATOM chain A within declared bounds minus one atom; bare TER/END. Not an intact deposited complex.",
            "source_writer_enabled": False, "production_enabled": False, "empirical_statistics_loaded": False}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    sources = [Path(__file__).resolve(), root/"scripts/branched_stereo_oracle.py",
               *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan | {"retained_directory": str(args.output), "python": platform.python_version(),
                     "rdkit": rdBase.rdkitVersion,
                     "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources}, "results": {}}
    for case in cases:
        folder = args.output/case["id"]
        folder.mkdir()
        lines = [l for l in originals[case["fixture"]].splitlines(keepends=True) if l.startswith(b"ATOM  ")
                 and l[21:22] == b"A" and (case["stop"] is None or int(l[22:26]) <= case["stop"])]
        source = b"".join(l for l in lines if not (l[22:26].strip().decode() == case["number"]
                          and l[12:16].strip().decode() == case["missing"]))+b"TER\nEND\n"
        (folder/"input.pdb").write_bytes(source)
        text_lines, atoms, _ = _layout(source)
        target = ComponentInstanceRef(1, "A", case["number"])
        observed = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position) for a in atoms)
        groups = list(dict.fromkeys(_component_ref(a) for a in atoms))
        links = []
        for left, right in zip(groups, groups[1:]):
            if left.chain_id == right.chain_id and observed_peptide_continuity(text_lines, left.chain_id,
                    ((left.sequence_number, left.insertion_code), (right.sequence_number, right.insertion_code))) == (True,):
                links.append((ExactAtomRef(left, "C"), ExactAtomRef(right, "N")))
        request = CompletionRequest(sha(source), target, case["residue_name"], observed, (case["missing"],), tuple(links))
        (folder/"request.json").write_text(json.dumps(asdict(request), indent=2)+"\n")
        before = request.digest
        result = complete_anchored_branch(request)
        (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
        if request.digest != before:
            raise ValueError("observed-anchor request mutated")
        checks = []
        for attempt in result.attempts:
            audit = audit_anchored_branch(request, attempt.candidate)
            if audit != attempt.audit:
                raise ValueError("independent anchored measurement audit disagrees")
            points = {o.atom.name: o.xyz for o in observed if o.atom.component == target}
            points.update({a.atom.name: a.xyz for a in attempt.candidate.additions})
            cip = rdkit_labels(points, case["residue_name"])
            expected = {"CA": "S", "CB": "R"} if case["residue_name"] == "THR" else {"CA": "S"}
            if cip != expected:
                raise ValueError("independent anchored proposal CIP disagrees")
            # Withheld coordinates are consulted only after proposal and selection.
            truth = next(a.position for a in _parse_atom_lines([l.decode() for l in lines])
                         if _component_ref(a) == target and a.name == case["missing"])
            checks.append({"rdkit_cip": cip, "heldout_exact_name_error_not_used_for_selection":
                           math.dist(attempt.candidate.additions[0].xyz, truth)})
        row = {"status": result.status, "reason": result.reason, "request_digest": request.digest,
               "attempts": asdict(result)["attempts"], "selected_index": result.selected_index,
               "independent_checks": checks, "observations_digest_unchanged": request.digest == before}
        (folder/"audit.json").write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
        row["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][case["id"]] = row
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({n: {"status": r["status"], "error": r["independent_checks"][0]["heldout_exact_name_error_not_used_for_selection"]
                         if r["independent_checks"] else None} for n, r in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
