"""Retain eight declared source-bound 1UBQ 1--72 reconstruction controls.

The explicit projection excludes partial-occupancy tail/water, without editing
occupancies. These are unchanged staggered controls, not authentic rotamer data
or production/whole-chemistry qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.source import (
    audit_materialized_source,
    materialize_source,
    plan_source,
)

CASES = ((28, ("CB",)), (46, ("CB",)), (20, ("OG",)), (20, ("CB", "OG")),
         (57, ("OG",)), (57, ("CB", "OG")), (65, ("OG",)), (65, ("CB", "OG")))


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    fixture = root/"tests/fixtures/protein_context/1UBQ.pdb"
    original = fixture.read_bytes()
    projection = b"".join(l for l in original.splitlines(keepends=True)
                          if l.startswith(b"ATOM  ") and int(l[22:26]) <= 72)+b"TER\nEND\n"
    (args.output/"original.pdb").write_bytes(original)
    (args.output/"projection.pdb").write_bytes(projection)
    plan = {"schema": "dvbfixer.private-protein-source-controls.v1", "fixture": str(fixture.relative_to(root)),
            "fixture_sha256": sha(original), "projection_sha256": sha(projection),
            "projection": "ATOM records for residues A1--72, unchanged bytes, followed by bare TER and END",
            "cases": [{"residue": n, "missing": names} for n, names in CASES],
            "sampling": "unchanged three-staggered-chi control; no authentic Dunbrack rows",
            "production_routing": False}
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    sources = [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py")),
               root/"src/dvbfixer/pdb_component_reconstruction.py"]
    report = plan | {"source_hashes": {str(p.relative_to(root)): sha(p.read_bytes()) for p in sources},
                     "retained_directory": str(args.output), "results": {}}
    for number, missing in CASES:
        name = f"A{number}-{'-'.join(missing)}"
        folder = args.output/name
        folder.mkdir()
        src = b"".join(l for l in projection.splitlines(keepends=True)
                       if not (l.startswith(b"ATOM  ") and int(l[22:26]) == number
                               and l[12:16].strip().decode() in missing))
        (folder/"input.pdb").write_bytes(src)
        req = plan_source(src, ComponentInstanceRef(1, "A", str(number)))
        result = complete_sidechain(req)
        (folder/"request.json").write_text(json.dumps(asdict(req), indent=2)+"\n")
        (folder/"attempts.json").write_text(json.dumps(asdict(result), indent=2)+"\n")
        if result.selected is None:
            raise RuntimeError(f"retained rejected control: {name}")
        materialized = materialize_source(src, req, result.selected)
        (folder/"completed.pdb").write_bytes(materialized.pdb)
        record = materialized.dat_record()
        record.save(folder/"completed.dat")
        if record.load(folder/"completed.dat").added_atoms != record.added_atoms:
            raise RuntimeError("canonical sidecar readback mismatch")
        generated = audit_materialized_source(src, req, result.selected, (folder/"completed.pdb").read_bytes())
        report["results"][name] = {"status": result.status, "selected_index": result.selected_index,
                                   "request_digest": req.digest, "input_sha256": sha(src),
                                   "output_sha256": materialized.output_sha256,
                                   "generated": [asdict(a) for a in generated],
                                   "readback_passed": True,
                                   "retained_files": {p.name: sha(p.read_bytes()) for p in sorted(folder.iterdir())}}
    report["all_eight_passed"] = len(report["results"]) == 8
    evidence = json.dumps(report, indent=2, allow_nan=False)+"\n"
    (args.output/"evidence.json").write_text(evidence)
    with args.summary.open("x") as handle:
        json.dump(report | {"raw_record_sha256": sha(evidence.encode())}, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: v["status"] for k, v in report["results"].items()}, indent=2))


if __name__ == "__main__":
    main()
