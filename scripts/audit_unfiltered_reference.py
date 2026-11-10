"""Reference-relative sanity control for a previously executed comparison gate.

No new reconstruction, candidate selection or threshold tuning. Held-out atoms
are used only now for post-execution reference audits and descriptive scoring.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

from scripts.compare_branched_completion import measure
from scripts.unfiltered_branched_comparison import reporting_scores, sha

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_source import plan_anchored_source
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.results.read_bytes()
    executed = json.loads(raw)
    root = Path(executed["retained_directory"]).resolve(strict=True)
    plan_root = Path(executed["plan_root"]).resolve(strict=True)
    raw_plan = read_private_file_bytes(plan_root, "plan.json", max_bytes=2_000_000)
    if sha(raw_plan) != executed["plan_sha256"]:
        raise ValueError("executed comparison plan changed")
    plan = json.loads(raw_plan)
    cases = {c["id"]: c for c in plan["cases"]}
    record = {"schema": "dvbfixer.unfiltered-native-reference-audit.v1", "comparison": str(args.results),
              "comparison_sha256": sha(raw), "plan_sha256": sha(raw_plan),
              "script_sha256": sha(Path(__file__).read_bytes()), "results": {}}
    for name, row in executed["results"].items():
        if row["status"] != "admitted":
            continue
        case = cases[name]
        for file, digest in row["retained_files"].items():
            if Path(file).name != file or sha(read_private_file_bytes(root/name, file, max_bytes=64*1024*1024)) != digest:
                raise ValueError("retained comparison artifact changed")
        source = read_private_file_bytes(root/name, "input.pdb", max_bytes=16*1024*1024)
        if sha(source) != plan["artifacts"][case["source"]]:
            raise ValueError("source differs from frozen plan")
        reference = read_private_file_bytes(plan_root, case["reference"], max_bytes=16*1024*1024)
        if sha(reference) != plan["artifacts"][case["reference"]]:
            raise ValueError("reference differs from frozen declaration")
        target = ComponentInstanceRef(**case["target"])
        anchored = case["mode"] == "anchored"
        source_plan = plan_anchored_source(source, target) if anchored else plan_branched_source(source, target)
        request = source_plan.request
        truth = {a.name: a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True))
                 if _component_ref(a) == target}
        native_additions = tuple(Addition(ExactAtomRef(target, n), n[0], truth[n]) for n in request.missing_atoms)
        audit, _ = measure(source, request, native_additions, anchored=anchored)
        attempts = json.loads((root/name/"own-attempts.json").read_text())
        own_reports = []
        for attempt in attempts["attempts"]:
            additions = tuple(Addition(ExactAtomRef(ComponentInstanceRef(**a["atom"]["component"]), a["atom"]["name"]),
                                       a["element"], tuple(a["xyz"])) for a in attempt["candidate"]["additions"])
            # Descriptive only: preserve attempt order, do not choose by truth.
            own_reports.append({"source_audit_accepted": attempt["audit"]["accepted"],
                                "scores_reporting_only": reporting_scores(request, additions, truth)})
        own_pass = row["own"].get("common_audit", {}).get("accepted", False)
        template_pass = row["template"].get("common_audit", {}).get("accepted", False)
        record["results"][name] = {"reference_common_audit": asdict(audit),
            "reference_self_rmsd_angstrom": 0., "own_attempts_reporting_only": own_reports,
            "reference_CG1_C_distance_angstrom": math.dist(truth["CG1"], truth["C"]) if request.residue_name == "VAL" else None,
            "native_gate_failure_not_a_reconstruction_error": not audit.accepted,
            "own_selected_common_pass": own_pass, "template_common_pass": template_pass}
    aggregates = {}
    for mode in ("paired", "anchored"):
        names = [n for n in executed["results"] if cases[n]["mode"] == mode]
        admitted = [n for n in names if executed["results"][n]["status"] == "admitted"]
        both = [n for n in admitted if record["results"][n]["own_selected_common_pass"]
                and record["results"][n]["template_common_pass"]]
        arms = {}
        for arm in ("own", "template"):
            same = [executed["results"][n][arm]["scores"] for n in both]
            arms[arm] = {"passes": sum(executed["results"][n].get(arm, {}).get("common_audit", {}).get("accepted", False) for n in admitted),
                "both_pass_exact_name_median_angstrom": statistics.median(s["exact_name_rmsd_angstrom"] for s in same) if same else None,
                "both_pass_generated_only_val_symmetry_median_angstrom": statistics.median(s["generated_only_val_symmetry_rmsd_angstrom"] for s in same) if same else None}
        aggregates[mode] = {"planned": len(names), "admitted": len(admitted), "both_pass": len(both),
            "reference_passes": sum(record["results"][n]["reference_common_audit"]["accepted"] for n in admitted),
            "template_protocol_blockers": sum(executed["results"][n].get("template", {}).get("status") == "execution_or_protocol_blocker" for n in admitted),
            "arms": arms}
    record["aggregate"] = aggregates
    with args.output.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(aggregates, indent=2))


if __name__ == "__main__":
    main()
