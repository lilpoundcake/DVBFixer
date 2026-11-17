#!/usr/bin/env python3
"""Digest-checked private review index/query; no selection, fitting or publication."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest


def scientific_task(request):
    raw = json.loads(request.to_json())
    selected = {key: raw[key] for key in ("target_sequences", "sequence_placements", "gaps", "fixed_atoms", "generated_atoms")}
    selected["normalized_source_sha256"] = request.normalized_pdb.sha256
    selected["retained_explicit_links"] = raw.get("retained_explicit_links", [])
    selected["heterogen_contexts"] = raw.get("heterogen_contexts", [])
    return hashlib.sha256(json.dumps(selected, sort_keys=True).encode()).hexdigest()


def filter_rows(rows, *, entry=None, backend=None, seed=None, stage=None):
    return [row for row in rows if (entry is None or row["entry_id"] == entry)
            and (backend is None or row["native_backend"] == backend)
            and (seed is None or row["native_seed"] == seed)
            and (stage is None or row["stage"] == stage)]


def candidate_id(row):
    identity = {key: row[key] for key in ("entry_id", "native_backend", "native_protocol", "native_seed",
        "optimizer_seed", "stage", "sha256", "scientific_task_sha256")}
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]


def build(reviews, output):
    rows, tasks, sources, bindings = [], {}, {}, {}
    for root in reviews:
        index_path = root / "review-index.json"
        index = json.loads(index_path.read_text())
        if index["kind"] != "private-unapproved-review-v1" or index["public_bundle"] or not index["all_manual_approvals_unset"]:
            raise ValueError("library only accepts explicitly unapproved private review exports")
        bindings[str(index_path)] = digest(index_path)
        for case in index["cases"]:
            review_path = root / case["review_path"]
            if not review_path.resolve().is_relative_to(root) or digest(review_path) != case["review_sha256"]:
                raise ValueError("review path/digest mismatch")
            recipe = json.loads(review_path.read_text())
            if recipe["selected_donor"] is not None or recipe["allowed_donors"] or recipe["manual_approval"] is not None:
                raise ValueError("never silently adopt or replace an existing manual selection")
            folder = review_path.parent
            for name, expected in recipe["files"].items():
                path = folder / name
                if Path(name).name != name or digest(path) != expected:
                    raise ValueError("private review payload changed")
                bindings[str(path)] = expected
            request = DiffusionRequest.from_json((folder / "request.json").read_text())
            task = scientific_task(request)
            entry = recipe["entry_id"]
            if (not re.fullmatch("[A-Za-z0-9]{4}", entry) or entry != case["entry_id"]
                    or recipe["native_backend"] not in {"Protenix", "MODELLER"}
                    or request.seeds != (recipe["native_seed"],)
                    or digest(folder / "01-observed-source.pdb") != request.normalized_pdb.sha256):
                raise ValueError("review entry/backend/seed/source provenance mismatch")
            if entry in tasks and tasks[entry] != task:
                raise ValueError("backend review cases do not have the same complete scientific mask/source")
            tasks[entry] = task
            sources.setdefault(entry, folder / "01-observed-source.pdb")
            generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
            for candidate in recipe["stages"]:
                if candidate["geometry_passed"] is None:
                    continue
                if candidate["manual_approval"] is not None:
                    raise ValueError("private review stage approval must remain explicit and unset")
                path = folder / candidate["file"]
                if (Path(candidate["file"]).name != candidate["file"] or path.suffix != ".pdb"
                        or recipe["files"].get(candidate["file"]) != candidate["sha256"]
                        or path.with_suffix(".dat").name not in recipe["files"]):
                    raise ValueError("review candidate/sidecar is not a bound local payload")
                dat = DatRecord.load(path.with_suffix(".dat"))
                actual = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in dat.added_atoms}
                if actual != generated or len(dat.added_atoms) != len(generated):
                    raise ValueError("review .dat reclassified source observations or lost generated atoms")
                row = {"entry_id": entry, "native_backend": recipe["native_backend"],
                    "native_protocol": recipe.get("native_protocol", "baseline"), "native_seed": recipe["native_seed"],
                    "optimizer_seed": recipe["optimizer_seed"], "stage": candidate["stage"],
                    "geometry_passed": candidate["geometry_passed"], "path": str(path), "sha256": digest(path),
                    "native_execution_status": recipe["native_full_inference_status"],
                    "scientific_task_sha256": task, "joint_group": recipe["joint_group"], "manual_approval": None}
                row["candidate_id"] = candidate_id(row)
                rows.append(row)
            bindings[str(review_path)] = digest(review_path)
    if len({row["candidate_id"] for row in rows}) != len(rows):
        raise ValueError("duplicate candidate lineage in library")
    output.mkdir(parents=True, exist_ok=False)
    with (output / "library-source.py").open("xb") as stream:
        stream.write(Path(__file__).read_bytes())
    columns = ("candidate_id", "entry_id", "native_backend", "native_protocol", "native_seed", "optimizer_seed", "stage", "geometry_passed", "native_execution_status", "path")
    with (output / "candidates.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    for entry, source in sources.items():
        commands = [f"load {json.dumps(str(source), ensure_ascii=False)}, dvbfixer_{entry}_observed_source"]
        for row in filter_rows(rows, entry=entry):
            name = "dvbfixer_"+re.sub("[^A-Za-z0-9_]", "_", f"{entry}_{row['native_backend']}_s{row['native_seed']}_{row['stage']}_{row['candidate_id']}")
            commands += [f"load {json.dumps(row['path'], ensure_ascii=False)}, {name}", f"disable {name}"]
        (output / f"{entry}-overlay.pml").write_text("\n".join(commands)+"\n", encoding="utf-8")
    (output / "README.md").write_text(
        "# Private unapproved proposal review\n\n"
        "Use candidates.tsv or the query command to filter by entry/backend/native seed/stage.\n"
        "Native and optimizer seeds, resource failures and geometry outcomes remain distinct.\n"
        "Per-entry overlays share the original source frame: no fitting is applied. All\n"
        "candidate objects are initially disabled; enable them manually for inspection.\n"
        "A joint all-gap proposal cannot be treated as independent per-gap donors.\n"
        "Geometry passes are not donor approvals, final evaluation or public publication.\n",
        encoding="utf-8")
    save(output / "selection-recipe.json", {"kind": "private-manual-selection-recipe-v1", "tasks": tasks,
        "selected_donors": {}, "allowed_donors": {}, "locked_joint_groups": {}, "manual_approval": None,
        "no_automatic_selection_from_geometry": True})
    save(output / "library.json", {"kind": "private-unapproved-review-library-v1", "candidates": rows,
        "tasks": tasks, "bindings": bindings, "builder_sha256": digest(Path(__file__)),
        "output_files": {p.name: digest(p) for p in output.iterdir() if p.is_file()},
        "source_snapshot": "library-source.py", "public_bundle": False, "manual_approvals": 0})
    print(json.dumps({"entries": len(tasks), "stage_records": len(rows),
                      "geometry_passed_stage_records": sum(row["geometry_passed"] for row in rows), "manual_approvals": 0}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "query"))
    parser.add_argument("--reviews", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--entry")
    parser.add_argument("--backend", choices=("Protenix", "MODELLER"))
    parser.add_argument("--seed", type=int)
    parser.add_argument("--stage")
    args = parser.parse_args()
    if args.mode == "build":
        build([path.resolve() for path in args.reviews], args.output.resolve())
    else:
        library = json.loads(args.library.read_text())
        for name, expected in library["output_files"].items():
            if Path(name).name != name or digest(args.library.parent / name) != expected:
                raise ValueError("library output changed; do not silently adopt a manual edit")
        for name, expected in library["bindings"].items():
            if digest(Path(name)) != expected:
                raise ValueError("library source artifact changed")
        print(json.dumps(filter_rows(library["candidates"], entry=args.entry, backend=args.backend, seed=args.seed, stage=args.stage),
                         indent=2, ensure_ascii=False))
