#!/usr/bin/env python3
"""Read actual full/projection context custody; distances never grant admission."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity
from dvbfixer.model.diffusion.gap_portfolio import _coordinates


def residue(atom):
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


def heavy(line):
    element = line[76:78].strip()
    return element not in {"H", "D"} if element else not line[12:16].strip().lstrip("0123456789").startswith(("H", "D"))


def context_facts(full_text, projection_text, request):
    full, records = _coordinates(full_text)
    projection, protein_records = _coordinates(projection_text)
    if not projection.keys() <= full.keys():
        raise ValueError("Original projection contains atoms absent from full context")
    for atom, point in projection.items():
        if (not np.array_equal(point, full[atom])
                or records[atom][17:20] != protein_records[atom][17:20]):
            raise ValueError("Original protein projection changed observations or residue identity")
    omitted = set(full)-set(projection)
    if any(not records[atom].startswith("HETATM") for atom in omitted):
        raise ValueError("Projection omitted an ATOM record; do not silently call it heterogen context")
    groups = defaultdict(list)
    counts = Counter()
    for atom in sorted(omitted):
        name = records[atom][17:20].strip()
        counts[name] += 1
        groups[(residue(atom), name)].append(atom)
    components = []
    for (identity, name), atoms in sorted(groups.items()):
        components.append({"residue": asdict(identity), "residue_name": name,
            "atoms": [asdict(atom) for atom in atoms], "atom_count": len(atoms),
            "heavy_atom_count": sum(heavy(records[atom]) for atom in atoms),
            "in_HOH_WAT_whitelist": name in {"HOH", "WAT"}})
    gap_rows = []
    omitted_heavy = tuple(sorted(atom for atom in omitted if heavy(records[atom])))
    trees = {}
    for label, atoms in (("all_omitted", omitted_heavy),
                         ("outside_HOH_WAT", tuple(atom for atom in omitted_heavy
                                                   if records[atom][17:20].strip() not in {"HOH", "WAT"}))):
        if atoms:
            trees[label] = (atoms, cKDTree(np.array([full[atom] for atom in atoms])))
    for index, gap in enumerate(request.gaps):
        anchors = {anchor for anchor in (gap.left_anchor, gap.right_anchor) if anchor is not None}
        observed = tuple(sorted(atom for atom in projection if residue(atom) in anchors and heavy(protein_records[atom])))
        if not observed:
            raise ValueError("Declared gap lacks original observed heavy anchor atoms")
        nearest = {}
        for label, (atoms, tree) in trees.items():
            distances, indices = tree.query(np.array([projection[atom] for atom in observed]))
            best = int(np.argmin(distances))
            partner = atoms[int(indices[best])]
            nearest[label] = {"distance_angstrom": float(distances[best]),
                "observed_anchor_atom": asdict(observed[best]), "omitted_atom": asdict(partner),
                "omitted_residue_name": records[partner][17:20].strip()}
        gap_rows.append({"gap_index": index, "chain": gap.chain,
            "anchor_observed_heavy_atoms": len(observed), "nearest": nearest})
    return {"original_full_atoms": len(full), "original_projection_atoms": len(projection),
        "omitted_atoms": len(omitted), "omitted_atom_counts_by_residue_name": dict(counts),
        "omitted_residue_instances": components, "gap_anchor_context": gap_rows,
        "projection_observations_exact_in_full_source_frame": True,
        "proximity_reference": "original pre-preparation observed anchor heavy atoms; not generated or repaired coordinates",
        "distance_does_not_establish_irrelevance_or_parameterization": True,
        "missing_region_positions_or_contacts_inferred": False, "scientific_admission": False}


def audit(inventory_path, output):
    inventory = json.loads(inventory_path.read_text())
    for name, expected in inventory["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Retained review input changed: {name}")
    output.mkdir(parents=True, exist_ok=False)
    (output / "entries").mkdir()
    rows, bindings = [], {str(inventory_path): digest(inventory_path)}
    for row in inventory["rows"]:
        entry = row["entry_id"]
        detail_path = inventory_path.parent / row["entry_report"]
        if not detail_path.resolve().is_relative_to(inventory_path.parent.resolve()):
            raise ValueError("Review report path escapes inventory")
        detail = json.loads(detail_path.read_text())
        request_path = Path(detail["request_path"])
        if digest(request_path) != detail["request_sha256"]:
            raise ValueError("Retained request changed")
        request = DiffusionRequest.from_json(request_path.read_text())
        work = request_path.parent
        metadata_path = work / "projection.json"
        metadata = json.loads(metadata_path.read_text())
        full_path, projection_path = work / "full-assembly-context.pdb", work / "protein-projection.pdb"
        for path, expected in ((full_path, metadata["full_context_sha256"]),
                               (projection_path, metadata["protein_projection_sha256"])):
            if digest(path) != expected:
                raise ValueError(f"Original full/projection custody mismatch: {path}")
        facts = context_facts(full_path.read_text(), projection_path.read_text(), request)
        expected_counts = metadata["omitted_heterogen_atom_counts"]
        if (facts["omitted_atom_counts_by_residue_name"] != expected_counts
                or expected_counts != detail["omitted_heterogen_atom_counts"]):
            raise ValueError("Actual context omission differs from the retained declaration")
        report = {"entry_id": entry, "scope": detail["current_scope_status"], **facts,
            "source_preparation_derivative": detail["prepared_source_derivative"],
            "current_prepared_source_path": detail["source_path"],
            "current_prepared_source_sha256": detail["source_sha256"],
            "prepared_source_not_substituted_for_original_frame": True,
            "original_full_context_path": str(full_path), "original_projection_path": str(projection_path)}
        save(output / "entries" / f"{entry}.json", report)
        outside = [item for item in facts["omitted_residue_instances"] if not item["in_HOH_WAT_whitelist"]]
        rows.append({"entry_id": entry, "scope": detail["current_scope_status"],
            "omitted_atoms": facts["omitted_atoms"], "omitted_residue_instances": len(facts["omitted_residue_instances"]),
            "outside_HOH_WAT_instances": len(outside), "outside_HOH_WAT_names": sorted({item["residue_name"] for item in outside}),
            "gap_count": len(request.gaps), "original_projection_custody_verified": True,
            "scientific_admission": False, "report": f"entries/{entry}.json"})
        for path in (detail_path, request_path, metadata_path, full_path, projection_path):
            bindings[str(path)] = digest(path)
    summary = {"entries": len(rows), "gap_anchor_reviews": sum(row["gap_count"] for row in rows),
        "entries_without_any_omitted_atoms": sum(row["omitted_atoms"] == 0 for row in rows),
        "entries_with_only_HOH_WAT_omissions": sum(row["omitted_atoms"] > 0 and row["outside_HOH_WAT_instances"] == 0 for row in rows),
        "entries_with_other_context": sum(row["outside_HOH_WAT_instances"] > 0 for row in rows),
        "actual_omission_counts_match_prior_declarations": True,
        "new_admissions_or_sources_or_restoration_runs": 0, "distance_cutoff_for_admission": None}
    with (output / "context-review.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    save(output / "context-audit.json", {"summary": summary, "rows": rows, "bindings": bindings,
        "auditor_sha256": digest(Path(__file__)),
        "reports_sha256": {row["report"]: digest(output / row["report"]) for row in rows}})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.inventory.resolve(), args.output.resolve())
