"""Declared VAL/THR paired-branch controls, not rotamer/native-quality qualification."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from benchmark_protein_completion import run_controls
from branched_stereo_oracle import rdkit_labels

from dvbfixer.protein_completion.branched import BRANCHES, named_beta_volume

CASES = ((5, "VAL", ("CG1", "CG2")), (5, "VAL", ("CB", "CG1", "CG2")),
         (26, "VAL", ("CG1", "CG2")), (26, "VAL", ("CB", "CG1", "CG2")),
         (7, "THR", ("OG1", "CG2")), (7, "THR", ("CB", "OG1", "CG2")),
         (55, "THR", ("OG1", "CG2")), (55, "THR", ("CB", "OG1", "CG2")))


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    original = (root/"tests/fixtures/protein_context/1UBQ.pdb").read_bytes()
    (args.output/"original.pdb").write_bytes(original)
    plan = {"schema": "dvbfixer.private-branched-controls.v1", "source_sha256": sha(original),
            "cases": CASES, "source_scope": "Explicit 1UBQ protein-only A1--72 projection",
            "admitted_masks": "paired distal branches, optionally including CB; no retained distal branch closure",
            "sampling": "three-staggered-chi control, no empirical statistics",
            "source_cohort_publication_capability_enabled": False}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    raw = run_controls(root, cases=CASES)
    raw_bytes = (json.dumps(raw, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"raw-controls.json").write_bytes(raw_bytes)
    report = plan | {"raw_record_sha256": sha(raw_bytes), "retained_directory": str(args.output), "results": {},
                     "source_hashes": raw["scripts_sha256"] | {
                         str(p.relative_to(root)): sha(p.read_bytes()) for p in (Path(__file__).resolve(),
                             root/"scripts/branched_stereo_oracle.py")}}
    for case_name, case in raw["cases"].items():
        request = case["request"]
        residue_name = request["residue_name"]
        fixed = {o["atom"]["name"]: tuple(o["xyz"]) for o in request["observations"]
                 if o["atom"]["component"] == request["residue"]}
        branch = BRANCHES[residue_name]
        stereo = []
        for attempt in case["result"]["attempts"]:
            points = fixed | {a["atom"]["name"]: tuple(a["xyz"]) for a in attempt["candidate"]["additions"]}
            labels = rdkit_labels(points, residue_name)
            expected = {"CA": "S", "CB": "R"} if residue_name == "THR" else {"CA": "S"}
            if labels != expected:
                raise ValueError("independent CIP audit disagrees; raw attempts retained")
            stereo.append({"rdkit_cip": labels, "named_beta_volume_angstrom3": named_beta_volume(
                points["CA"], points["CB"], points[branch.first], points[branch.second]),
                "branch_policy": asdict(branch)})
        selected = case["result"]["selected_index"]
        report["results"][case_name] = {"status": case["result"]["status"], "reason": case["result"]["reason"],
            "request_digest": case["result"]["request_digest"], "revision": case["result"]["revision"],
            "attempts": case["result"]["attempts"], "independent_stereo": stereo, "selected_index": selected,
            "heldout_scores_not_used_for_selection": case["scores_not_used_for_selection"],
            "selected_exact_name_rmsd_angstrom": None if selected is None else
                case["scores_not_used_for_selection"][selected]["heldout_rmsd_angstrom"],
            "max_observed_movement_angstrom": case["max_observed_movement_angstrom"]}
    data = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(data)
    with args.summary.open("x") as handle:
        json.dump(report | {"retained_evidence_sha256": sha(data)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({name: {"status": row["status"], "selected_rmsd": row["selected_exact_name_rmsd_angstrom"]}
                      for name, row in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
