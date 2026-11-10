"""Retain declared-link admission controls on an explicit real 1DQJ projection."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.pdb_component_reconstruction import _parse_atom_lines
from dvbfixer.protein_completion.links import parse_declared_links
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import (
    audit_materialized_source,
    materialize_source,
    plan_source,
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    fixture = root/"tests/fixtures/regressions/1DQJ.pdb"
    original = fixture.read_bytes()
    lines = original.decode("ascii").splitlines(keepends=True)
    ssbond = next(l for l in lines if l.startswith("SSBOND   1"))
    projected = []
    for number in (9, 23, 88):
        projected.extend(l for l in lines if l.startswith("ATOM  ") and l[21] == "A" and int(l[22:26]) == number
                         and not (number == 9 and l[12:16].strip() == "CB"))
        projected.append("TER\n")
    base = "".join(projected)+"END\n"
    atoms = _parse_atom_lines(base.splitlines(keepends=True))
    target = ComponentInstanceRef(1, "A", "9")
    n = next(a.serial for a in atoms if a.sequence_number == "9" and a.name == "N")
    ca = next(a.serial for a in atoms if a.sequence_number == "9" and a.name == "CA")
    cases = (
        ("unannotated-control", base, target, True),
        ("deposited-ssbond", ssbond+base, target, True),
        ("ssbond-and-target-conect", ssbond+base.replace("END\n", f"CONECT{n:5d}{ca:5d}\nEND\n"), target, True),
        ("nonidentity-symmetry", ssbond[:59]+"  2555"+ssbond[65:]+base, target, False),
        ("linked-target-refusal", ssbond+base, ComponentInstanceRef(1, "A", "23"), False),
    )
    plan = {"schema": "dvbfixer.private-protein-link-controls.v1", "fixture": str(fixture.relative_to(root)),
            "fixture_sha256": sha(original),
            "projection": "Unmodified A9 (CB withheld), A23 and A88 ATOM bytes; explicit bare TER boundaries; deposited first SSBOND where declared",
            "cases": [{"id": name, "target": asdict(t), "expected_admission": accepted}
                      for name, _, t, accepted in cases], "production_enabled": False,
            "authentic_rotamer_statistics_loaded": False}
    (args.output/"original.pdb").write_bytes(original)
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    sources = [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]
    report = plan | {"source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources},
                     "retained_directory": str(args.output), "results": {}}
    selected_coordinates = None
    for name, text, selected, expected_admission in cases:
        folder = args.output/name
        folder.mkdir()
        src = text.encode("ascii")
        (folder/"input.pdb").write_bytes(src)
        try:
            request = plan_source(src, selected)
        except ValueError as exc:
            finding = {"admitted": False, "refusal": str(exc)}
            (folder/"refusal.json").write_text(json.dumps(finding, indent=2)+"\n")
            if expected_admission:
                raise
        else:
            if not expected_admission:
                raise RuntimeError(f"unexpected admission: {name}")
            result = complete_sidechain(request)
            (folder/"request.json").write_text(json.dumps(asdict(request), indent=2)+"\n")
            (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
            if result.selected is None:
                raise RuntimeError(f"unexpected completion refusal: {name}")
            output = materialize_source(src, request, result.selected)
            (folder/"completed.pdb").write_bytes(output.pdb)
            output.dat_record().save(folder/"completed.dat")
            readback = audit_materialized_source(src, request, result.selected, (folder/"completed.pdb").read_bytes())
            if selected_coordinates is not None and readback != selected_coordinates:
                raise RuntimeError("metadata changed generated coordinates")
            selected_coordinates = readback
            parsed = parse_declared_links(text.splitlines(keepends=True), atoms)
            finding = {"admitted": True, "readback_passed": True, "request_digest": request.digest,
                       "links": [[asdict(a) for a in pair] for pair in request.links],
                       "declared_metadata": [asdict(a) for a in parsed],
                       "generated": [asdict(a) for a in readback], "output_sha256": output.output_sha256}
        finding["retained_files"] = {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}
        report["results"][name] = finding
    report["all_expected_outcomes_matched"] = True
    raw = (json.dumps(report, indent=2, allow_nan=False)+"\n").encode()
    (args.output/"evidence.json").write_bytes(raw)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(raw)}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({name: finding["admitted"] for name, finding in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
