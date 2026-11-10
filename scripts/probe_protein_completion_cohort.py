"""Retain bounded joint-completion controls; not native packing qualification."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _format_generated_atom
from dvbfixer.protein_completion.cohort import (
    audit_materialized_cohort,
    complete_cohort,
    materialize_cohort,
    plan_cohort,
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def synthetic(name: str, offset: tuple[float, float, float], *, water: bool = False) -> bytes:
    points = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.), "C": (.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    lines = ["HEADER    SYNTHETIC JOINT COMPLETION CONTROL\n"]
    for j, chain in enumerate(("A", "B")):
        for i, (atom, p) in enumerate(points.items(), 1):
            xyz = p if j == 0 else (p[0]+offset[0], -p[1]+offset[1], -p[2]+offset[2])
            line = _format_generated_atom(serial=j*10+i, name=atom, element=atom[0], position=xyz,
                                          residue_name=name, chain_id=chain, sequence_number="1", insertion_code="")
            lines.append("ATOM  "+line[6:])
        lines.append("TER\n")
    if water:
        lines.append(_format_generated_atom(serial=25, name="O", element="O", position=(.589, -.6, 1.405),
                                             residue_name="HOH", chain_id="W", sequence_number="1", insertion_code=""))
    return ("".join(lines)+"END\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    fixture = root/"tests/fixtures/protein_context/1UBQ.pdb"
    original = fixture.read_bytes()
    (args.output/"original-1UBQ.pdb").write_bytes(original)
    cases = []
    for label, numbers in (("ubq-two-ala", (28, 46)), ("ubq-three-ser", (20, 57, 65)),
                           ("ubq-four-mixed", (20, 28, 46, 57))):
        src = b"".join(l for l in original.splitlines(keepends=True)
                       if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72
                       and not (int(l[22:26]) in numbers and l[12:16].strip() in {b"CB", b"OG"}))+b"TER\nEND\n"
        cases.append((label, src, tuple(ComponentInstanceRef(1, "A", str(n)) for n in numbers), "completed"))
    targets = (ComponentInstanceRef(1, "A", "1"), ComponentInstanceRef(1, "B", "1"))
    cases.extend((
        ("synthetic-ser-joint-options", synthetic("SER", (-3., 0., -4.)), targets, "completed"),
        ("synthetic-ala-no-joint-solution", synthetic("ALA", (0., -1.538, -5.)), targets, "joint-rejected"),
        ("synthetic-rounding-refusal", synthetic("ALA", (30., 0., 0.), water=True), targets, "rounded-rejected"),
    ))
    plan = {"schema": "dvbfixer.private-protein-cohort-controls.v1", "fixture_sha256": sha(original),
            "projection": "Unmodified 1UBQ ATOM A1--72, masks withheld, bare TER/END; synthetic cases labelled separately",
            "cases": [{"id": label, "input_sha256": sha(src), "targets": [asdict(t) for t in ts], "expected": expected}
                      for label, src, ts, expected in cases],
            "authentic_rotamer_data_loaded": False, "production_enabled": False}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    sources = [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan | {"retained_directory": str(args.output),
                     "source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources}, "results": {}}
    for label, src, ts, expected in cases:
        folder = args.output/label
        folder.mkdir()
        (folder/"input.pdb").write_bytes(src)
        request = plan_cohort(src, ts)
        (folder/"request.json").write_text(json.dumps(asdict(request), indent=2)+"\n")
        result = complete_cohort(src, request)
        (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
        finding = {"plan_digest": request.digest, "raw_status": result.status,
                   "joint_attempts": len(result.attempts),
                   "joint_passing": sum(a.accepted for a in result.attempts), "selected_index": result.selected_index}
        if result.selected is None:
            if expected != "joint-rejected":
                raise RuntimeError(f"unexpected raw refusal: {label}")
            finding["outcome"] = "joint-rejected"
        else:
            try:
                output = materialize_cohort(src, request, result.selected)
            except ValueError as exc:
                if expected != "rounded-rejected" or "rounded joint PDB rejected" not in str(exc):
                    raise
                finding.update(outcome="rounded-rejected", refusal=str(exc))
                (folder/"refusal.json").write_text(json.dumps(finding, indent=2)+"\n")
            else:
                if expected != "completed":
                    raise RuntimeError(f"unexpected successful materialization: {label}")
                (folder/"completed.pdb").write_bytes(output.pdb)
                record = output.dat_record()
                record.save(folder/"completed.dat")
                if record.load(folder/"completed.dat").added_atoms != record.added_atoms:
                    raise RuntimeError("generated sidecar readback disagrees")
                added = audit_materialized_cohort(src, request, result.selected, (folder/"completed.pdb").read_bytes())
                finding.update(outcome="completed", output_sha256=output.output_sha256,
                               generated=[asdict(a) for a in added], independent_readback_passed=True)
        finding["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][label] = finding
    report["all_expected_outcomes_matched"] = True
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: v["outcome"] for k, v in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
