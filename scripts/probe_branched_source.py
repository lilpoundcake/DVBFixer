"""Predeclared explicit VAL/THR source projections and independent rounded audits.

No production registration, empirical statistics, observed-coordinate repair or
automatic occupancy normalization. Admission/raw/rounded refusals are retained.
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

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.anchored_branch import complete_anchored_branch
from dvbfixer.protein_completion.anchored_source import (
    ANCHORED_SOURCE_POLICY,
    audit_anchored_materialized_source,
    materialize_anchored_source,
    plan_anchored_source,
)
from dvbfixer.protein_completion.branched import BRANCHES
from dvbfixer.protein_completion.branched_source import (
    BRANCHED_SOURCE_POLICY,
    audit_branched_materialized_source,
    materialize_branched_source,
    plan_branched_source,
)
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import SourceCompletionError

FIXTURES = {"1UBQ": "tests/fixtures/protein_context/1UBQ.pdb",
            "1FR2": "tests/fixtures/regressions/1FR2.pdb",
            "1A3N": "tests/fixtures/protein_context/1A3N.pdb"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--mode", choices=("paired", "anchored"), default="paired")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    anchored = args.mode == "anchored"
    source_policy = ANCHORED_SOURCE_POLICY if anchored else BRANCHED_SOURCE_POLICY
    source_planner = plan_anchored_source if anchored else plan_branched_source
    complete = complete_anchored_branch if anchored else complete_sidechain
    materialize = materialize_anchored_source if anchored else materialize_branched_source
    audit_output = audit_anchored_materialized_source if anchored else audit_branched_materialized_source
    root = Path(__file__).resolve().parents[1]
    originals = {name: (root/path).read_bytes() for name, path in FIXTURES.items()}
    for name, data in originals.items():
        (args.output/f"original-{name}.pdb").write_bytes(data)
    cases = []
    for fixture, stop, residues in (("1UBQ", 72, (("VAL", 5), ("VAL", 26), ("THR", 7), ("THR", 55))),
                                     ("1FR2", None, (("VAL", 19), ("THR", 11))),
                                     ("1A3N", 59, (("VAL", 1), ("THR", 8)))):
        for residue, number in residues:
            branch = BRANCHES[residue]
            masks = ((branch.first,), (branch.second,)) if anchored else (
                (branch.first, branch.second), ("CB", branch.first, branch.second))
            for missing in masks:
                cases.append({"id": f"{fixture}-A1-{stop or 'end'}-{residue}{number}-{'-'.join(missing)}",
                              "fixture": fixture, "stop": stop, "residue_name": residue,
                              "number": str(number), "missing": missing, "admission_control": False})
    # Full chain A has partial occupancies, inventoried before proposal execution.
    # Retain explicit refusal controls, never rewrite occupancy to make them pass.
    for residue, number in (("VAL", 1), ("THR", 8)):
        branch = BRANCHES[residue]
        cases.append({"id": f"1A3N-full-A-{residue}{number}-occupancy-refusal", "fixture": "1A3N",
                      "stop": None, "residue_name": residue, "number": str(number),
                      "missing": (branch.first,) if anchored else ("CB", branch.first, branch.second),
                      "admission_control": True})
    plan = {"schema": "dvbfixer.private-anchored-source-controls.v1" if anchored else "dvbfixer.private-branched-source-controls.v1",
            "policy": source_policy, "source_mode": args.mode,
            "projection": "Unmodified ATOM chain A records within declared bounds, masks withheld; bare TER/END. Not intact deposited complex.",
            "fixture_sha256": {name: sha(data) for name, data in originals.items()}, "cases": cases,
            "default_cohort_evidence_capability_expanded": False, "production_enabled": False,
            "empirical_rotamer_data_loaded": False}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    sources = [Path(__file__).resolve(), root/"scripts/branched_stereo_oracle.py",
               *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan | {"retained_directory": str(args.output), "python": platform.python_version(),
                     "rdkit": rdBase.rdkitVersion,
                     "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources}, "results": {}}
    for case in cases:
        folder = args.output/case["id"]
        folder.mkdir()
        lines = [l for l in originals[case["fixture"]].splitlines(keepends=True)
                 if l.startswith(b"ATOM  ") and l[21:22] == b"A"
                 and (case["stop"] is None or int(l[22:26]) <= case["stop"])]
        source = b"".join(l for l in lines if not (l[22:26].strip().decode() == case["number"]
                          and l[12:16].strip().decode() in case["missing"]))+b"TER\nEND\n"
        (folder/"input.pdb").write_bytes(source)
        target = ComponentInstanceRef(1, "A", case["number"])
        finding = {}
        try:
            request = source_planner(source, target)
        except ValueError as exc:
            finding = {"status": "admission_refused", "reason": str(exc)}
        else:
            if case["admission_control"]:
                raise ValueError("partial-occupancy control was unexpectedly admitted")
            (folder/"source-plan.json").write_text(json.dumps(asdict(request), indent=2)+"\n")
            result = complete(request.request)
            (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
            if result.selected is None:
                finding = {"status": "raw_rejected", "reason": result.reason}
            else:
                try:
                    output = materialize(source, request, result.selected)
                except SourceCompletionError as exc:
                    finding = {"status": "materialization_refused", "reason": str(exc)}
                else:
                    checked = audit_output(source, request, result.selected, output.pdb)
                    atoms = _parse_atom_lines(output.pdb.decode().splitlines(keepends=True))
                    groups = {}
                    for atom in atoms:
                        if atom.residue_name == "THR" or _component_ref(atom) == target:
                            groups.setdefault(_component_ref(atom), {"name": atom.residue_name, "points": {}})["points"][atom.name] = atom.position
                    cip = []
                    for residue, group in groups.items():
                        labels = rdkit_labels(group["points"], group["name"])
                        expected = {"CA": "S", "CB": "R"} if group["name"] == "THR" else {"CA": "S"}
                        if labels != expected:
                            raise ValueError("independent rounded global THR/target CIP disagrees")
                        cip.append({"residue": asdict(residue), "labels": labels})
                    (folder/"completed.pdb").write_bytes(output.pdb)
                    record = output.dat_record()
                    record.save(folder/"completed.dat", verbose=False)
                    if DatRecord.load(folder/"completed.dat").to_dict() != record.to_dict():
                        raise ValueError("generated-only sidecar readback disagrees")
                    # Held-out reference coordinates enter reporting only now.
                    truth_atoms = _parse_atom_lines([l.decode() for l in lines])
                    truth = {a.name: a.position for a in truth_atoms if _component_ref(a) == target}
                    errors = [math.dist(a.xyz, truth[a.atom.name]) for a in checked.generated]
                    finding = {"status": "accepted", "plan_digest": request.digest,
                               "independent_readback": asdict(checked), "rdkit_rounded_checks": cip,
                               "selected_index": result.selected_index, "generated": len(checked.generated),
                               "heldout_exact_name_rmsd_not_used_for_selection": math.sqrt(sum(d*d for d in errors)/len(errors))}
        (folder/"audit.json").write_text(json.dumps(finding, indent=2, allow_nan=False)+"\n")
        finding["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][case["id"]] = finding
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({name: row["status"] for name, row in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
