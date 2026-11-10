"""Retain explicit mixed-method joint controls, not empirical/native qualification."""
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
from dvbfixer.pdb_component_reconstruction import (
    _component_ref,
    _format_generated_atom,
    _parse_atom_lines,
)
from dvbfixer.protein_completion.branched import BRANCHES, TETRAHEDRAL_ANGLE
from dvbfixer.protein_completion.geometry import extend, l_cb
from dvbfixer.protein_completion.mixed_cohort import (
    MIXED_POLICY,
    audit_materialized_mixed_cohort,
    complete_mixed_cohort,
    materialize_mixed_cohort,
    plan_mixed_cohort,
)
from dvbfixer.protein_completion.source import SourceCompletionError

FIXTURES = {"1UBQ": "tests/fixtures/protein_context/1UBQ.pdb",
            "1FR2": "tests/fixtures/regressions/1FR2.pdb", "1A3N": "tests/fixtures/protein_context/1A3N.pdb"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def synthetic(specs, *, water=False):
    lines, targets = [], []
    for j, (name, chain, mask, offset, rotate) in enumerate(specs):
        xyz = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
        cb = xyz["CB"] = l_cb(xyz["N"], xyz["CA"], xyz["C"])
        if name in BRANCHES:
            branch = BRANCHES[name]
            xyz[branch.first] = extend(xyz["N"], xyz["CA"], cb, branch.first_length, TETRAHEDRAL_ANGLE, 60.)
            xyz[branch.second] = extend(xyz["N"], xyz["CA"], cb, branch.second_length,
                                        TETRAHEDRAL_ANGLE, 60.+branch.second_chi_offset)
        targets.append(ComponentInstanceRef(1, chain, "1"))
        for i, (atom, p) in enumerate(xyz.items(), 1):
            if atom in mask:
                continue
            point = tuple((-x if rotate and k in {1, 2} else x)+offset[k] for k, x in enumerate(p))
            line = _format_generated_atom(serial=j*20+i, name=atom, element=atom[0], position=point,
                                          residue_name=name, chain_id=chain, sequence_number="1", insertion_code="")
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    if water:
        lines.append(_format_generated_atom(serial=80, name="O", element="O", position=(-2.169, -1.093, -2.192),
                                             residue_name="HOH", chain_id="W", sequence_number="1", insertion_code=""))
    return ("".join(lines)+"END\n").encode(), tuple(targets)


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
    declarations = (
        ("ubq-four-mixed", "1UBQ", 72, ((28, ("CB",)), (26, ("CG1", "CG2")), (7, ("OG1",)), (57, ("OG",))), "measure"),
        ("ubq-two-thr", "1UBQ", 72, ((7, ("CB", "OG1", "CG2")), (55, ("CG2",))), "measure"),
        ("fr2-paired-anchored", "1FR2", None, ((19, ("CG1", "CG2")), (11, ("OG1",))), "measure"),
        ("a3n-whole-val-anchored-thr", "1A3N", 59, ((1, ("CB", "CG1", "CG2")), (8, ("OG1",))), "measure"),
        ("a3n-occupancy-refusal", "1A3N", None, ((1, ("CG1",)), (8, ("OG1",))), "admission_refused"),
    )
    cases = []
    for label, fixture, stop, masks, expected in declarations:
        lines = [l for l in originals[fixture].splitlines(keepends=True) if l.startswith(b"ATOM  ")
                 and l[21:22] == b"A" and (stop is None or int(l[22:26]) <= stop)]
        selected = dict(masks)
        src = b"".join(l for l in lines if l[12:16].strip().decode() not in selected.get(int(l[22:26]), ()))+b"TER\nEND\n"
        truth = b"".join(lines)+b"TER\nEND\n"
        cases.append((label, src, tuple(ComponentInstanceRef(1, "A", str(n)) for n in selected), expected, truth))
    for label, first_mask, expected in (("synthetic-mixed-choice", ("CB", "OG1", "CG2"), "accepted"),
                                        ("synthetic-two-anchors-refusal", ("CG2",), "joint_rejected")):
        src, targets = synthetic((("THR", "A", first_mask, (0., 0., 0.), False),
                                  ("THR", "a", ("CG2",), (.8, -4.5, -4.95), True)))
        cases.append((label, src, targets, expected, None))
    src, targets = synthetic((("THR", "A", ("OG1",), (0., 0., 0.), False),
                              ("ALA", "B", ("CB",), (30., 0., 0.), False)), water=True)
    cases.append(("synthetic-rounding-refusal", src, targets, "rounded_refused", None))
    plan_record = {"schema": "dvbfixer.private-mixed-cohort-controls.v1", "policy": MIXED_POLICY,
        "fixture_sha256": {name: sha(data) for name, data in originals.items()},
        "projection": "Original ATOM chain A within declared bounds, declared masks withheld, bare TER/END; synthetic cases labelled. Not intact complexes.",
        "real_declarations": declarations,
        "cases": [{"id": label, "source_sha256": sha(src), "targets": [asdict(t) for t in ts], "expected": expected}
                  for label, src, ts, expected, _ in cases], "production_enabled": False, "old_evidence_format_expanded": False,
        "empirical_statistics_loaded": False}
    (args.output/"plan.json").write_text(json.dumps(plan_record, indent=2)+"\n")
    sources = [Path(__file__).resolve(), root/"scripts/branched_stereo_oracle.py",
               *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan_record | {"retained_directory": str(args.output), "python": platform.python_version(),
        "rdkit": rdBase.rdkitVersion, "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources}, "results": {}}
    for label, src, targets, expected, truth in cases:
        folder = args.output/label
        folder.mkdir()
        (folder/"input.pdb").write_bytes(src)
        row = {}
        try:
            plan = plan_mixed_cohort(src, targets)
        except ValueError as exc:
            row = {"status": "admission_refused", "reason": str(exc)}
        else:
            (folder/"source-plan.json").write_text(json.dumps(asdict(plan), indent=2)+"\n")
            result = complete_mixed_cohort(src, plan)
            (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
            row = {"status": "joint_rejected", "reason": result.reason, "plan_digest": plan.digest,
                   "modes": [t.mode for t in plan.targets], "joint_attempts": len(result.attempts),
                   "joint_passes": sum(a.accepted for a in result.attempts), "selected_index": result.selected_index}
            if result.selected is not None:
                try:
                    out = materialize_mixed_cohort(src, plan, result.selected)
                except SourceCompletionError as exc:
                    if not str(exc).startswith("rounded mixed PDB rejected:"):
                        raise
                    row.update(status="rounded_refused", reason=str(exc))
                else:
                    checked = audit_materialized_mixed_cohort(src, plan, result.selected, out.pdb)
                    atoms = _parse_atom_lines(out.pdb.decode().splitlines(keepends=True))
                    source_names = {_component_ref(a): a.residue_name for a in atoms}
                    cip = []
                    wanted = {r for r, name in source_names.items() if name == "THR"} | {
                        t.request.residue for t in plan.targets if t.request.residue_name in {"VAL", "THR"}}
                    for residue, name in source_names.items():
                        if residue not in wanted:
                            continue
                        xyz = {a.name: a.position for a in atoms if _component_ref(a) == residue}
                        labels = rdkit_labels(xyz, name)
                        if labels != ({"CA": "S", "CB": "R"} if name == "THR" else {"CA": "S"}):
                            raise ValueError("independent mixed target CIP disagrees")
                        cip.append({"residue": asdict(residue), "labels": labels})
                    (folder/"completed.pdb").write_bytes(out.pdb)
                    record = out.dat_record()
                    record.save(folder/"completed.dat", verbose=False)
                    if DatRecord.load(folder/"completed.dat").to_dict() != record.to_dict():
                        raise ValueError("mixed generated sidecar readback differs")
                    scores = []
                    if truth is not None:
                        # Held-out truth is used only after selection and readback.
                        reference = {(_component_ref(a), a.name): a.position
                                     for a in _parse_atom_lines(truth.decode().splitlines(keepends=True))}
                        for target in plan.targets:
                            added = [a for a in checked.generated if a.atom.component == target.request.residue]
                            errors = [math.dist(a.xyz, reference[(a.atom.component, a.atom.name)]) for a in added]
                            scores.append({"residue": asdict(target.request.residue),
                                "residue_name": source_names[target.request.residue],
                                "exact_name_rmsd_not_used_for_selection": math.sqrt(sum(d*d for d in errors)/len(errors))})
                    row.update(status="accepted", independent_readback=asdict(checked),
                               rdkit_target_cip=cip, heldout_report_only=scores)
        if expected != "measure" and row["status"] != expected:
            raise ValueError(f"declared control outcome mismatch: {label}: {row}")
        (folder/"audit.json").write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
        row["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][label] = row
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({name: row["status"] for name, row in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
