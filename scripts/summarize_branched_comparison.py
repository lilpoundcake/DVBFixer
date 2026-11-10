"""Report paired metrics without hiding failed proposals or mixing accepted subsets."""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from collections import Counter
from pathlib import Path


def summarize(record: dict) -> dict:
    result = {}
    for mode in ("paired", "anchored"):
        rows = [r for name, r in record["results"].items() if name.startswith(mode+"-")]
        if not rows:
            raise ValueError("both declared comparison modes required")
        both = [r for r in rows if all(r[a]["common_audit"]["accepted"] for a in ("own", "template"))]
        arms = {}
        for arm in ("own", "template"):
            raw = [r[arm]["scores"]["exact_name_rmsd_angstrom"] for r in rows]
            passed = [r[arm]["scores"]["exact_name_rmsd_angstrom"] for r in rows if r[arm]["common_audit"]["accepted"]]
            same_subset = [r[arm]["scores"]["exact_name_rmsd_angstrom"] for r in both]
            arms[arm] = {"common_passes": len(passed), "total": len(rows),
                         "raw_proposal_median_rmsd_angstrom_including_rejections": statistics.median(raw),
                         "arm_specific_accepted_subset_median_rmsd_angstrom": statistics.median(passed) if passed else None,
                         "both_pass_same_subset_median_rmsd_angstrom": statistics.median(same_subset) if same_subset else None}
        result[mode] = {"arms": arms, "both_pass_count": len(both),
            "own_lower_raw_error_count_including_rejections": sum(r["own"]["scores"]["exact_name_rmsd_angstrom"] < r["template"]["scores"]["exact_name_rmsd_angstrom"] for r in rows),
            "own_lower_error_count_in_both_pass_subset": sum(r["own"]["scores"]["exact_name_rmsd_angstrom"] < r["template"]["scores"]["exact_name_rmsd_angstrom"] for r in both)}
    rows = list(record["results"].values())
    template = [r["template"] for r in rows]
    return {"modes": result, "template_helper_attempt_counts": dict(Counter(str(r["helper_attempts"]) for r in template)),
        "template_raw_non_L_or_nearplanar_ca_gate_count": sum("non_L_or_nearplanar_ca" in r["common_audit"]["raw"]["findings"] for r in template),
        "template_rounded_independent_CA_R_count": sum(r["common_audit"].get("rdkit_target_cip", {}).get("CA") == "R" for r in template),
        "template_raw_non_R_or_nearplanar_thr_cb_gate_count": sum("thr_cb_not_R_or_nearplanar" in r["common_audit"]["raw"]["findings"] for r in template),
        "template_rounded_independent_CB_S_count": sum(r["common_audit"].get("rdkit_target_cip", {}).get("CB") == "S" for r in template),
        "template_raw_val_name_orientation_count": sum("val_branch_name_orientation" in r["common_audit"]["raw"]["findings"] for r in template),
        "template_temporary_max_observed_motion_angstrom": max(r["observed_temporary_max_movement_angstrom"] for r in template),
        "not_independent_scientific_seeds": True, "not_full_legacy_pipeline": True,
        "retrospective_own_accepted_cases": True, "production_enabled": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = args.evidence.read_bytes()
    result = {"schema": "dvbfixer.retrospective-branched-comparison-summary.v1",
              "evidence": str(args.evidence), "evidence_sha256": hashlib.sha256(data).hexdigest(),
              "summary_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              **summarize(json.loads(data))}
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
