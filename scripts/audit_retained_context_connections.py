#!/usr/bin/env python3
"""Separate deposited assembly membership from actual projected connection loss."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import gemmi
import numpy as np
from scripts.audit_retained_construct_annotations import category_rows
from scripts.build_retained_admission_reviews import checked
from scripts.build_retained_gap_reviews import unique
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.scope import _parse_coordinate_records, _parse_explicit_links


def deposited_facts(block, assembly_id="1"):
    polymers = unique(category_rows(block, "_entity_poly."), "entity_id")
    chains = unique(category_rows(block, "_struct_asym."), "id")
    protein = {label for label, row in chains.items()
               if polymers.get(row["entity_id"], {}).get("type") == "polypeptide(L)"}
    generators = [row for row in category_rows(block, "_pdbx_struct_assembly_gen.")
                  if row["assembly_id"] == assembly_id]
    if not generators or not protein:
        raise ValueError("Missing declared assembly or protein membership")
    members = {label.strip() for row in generators for label in row["asym_id_list"].split(",")}
    if not members or not members <= chains.keys():
        raise ValueError("Assembly references unknown label chains")
    connections = category_rows(block, "_struct_conn.")
    unique(connections, "id")
    rows = []
    for connection in connections:
        labels = [connection[f"ptnr{i}_label_asym_id"] for i in (1, 2)]
        if not all(label in chains for label in labels):
            raise ValueError("Connection references unknown label chain")
        cross = (labels[0] in protein) != (labels[1] in protein)
        rows.append({"id": connection["id"], "depositor_connection_type": connection["conn_type_id"],
                     "raw_annotation": connection,
                     "protein_endpoint_flags": [label in protein for label in labels],
                     "assembly_member_label_flags": [label in members for label in labels],
                     "crosses_protein_nonprotein_boundary": cross,
                     "both_endpoint_labels_in_declared_assembly": all(label in members for label in labels),
                     "actual_operator_pair_or_coordinate_mapping_established": False})
    cross = [row for row in rows if row["crosses_protein_nonprotein_boundary"]]
    return {"assembly_id": assembly_id, "assembly_generators": generators,
            "protein_label_chains": sorted(protein), "assembly_member_label_chains": sorted(members),
            "connections": rows,
            "source_cross_boundary_types": dict(Counter(row["depositor_connection_type"] for row in cross)),
            "assembly_member_label_cross_boundary_types": dict(Counter(row["depositor_connection_type"] for row in cross
                if row["both_endpoint_labels_in_declared_assembly"])),
            "outside_assembly_member_label_cross_boundary_types": dict(Counter(row["depositor_connection_type"] for row in cross
                if not row["both_endpoint_labels_in_declared_assembly"])),
            "membership_does_not_establish_operator_pair_or_present_atoms": True,
            "absence_of_annotation_does_not_prove_no_interaction": True}


def projection_link_facts(full_text, projection_text):
    full_lines, projected_lines = full_text.splitlines(), projection_text.splitlines()
    full, _, _ = _parse_coordinate_records(full_lines)
    projection, _, _ = _parse_coordinate_records(projected_lines)
    full_map, projection_map = [{record.identity: record for record in records} for records in (full, projection)]
    if len(full_map) != len(full) or len(projection_map) != len(projection):
        raise ValueError("Duplicate coordinate identity")
    if not projection_map.keys() <= full_map.keys():
        raise ValueError("Projection contains an unknown atom")
    full_points, _ = _coordinates(full_text)
    projection_points, _ = _coordinates(projection_text)
    for atom, record in projection_map.items():
        original = full_map[atom]
        if not np.array_equal(projection_points[atom], full_points[atom]) or record.residue_name != original.residue_name:
            raise ValueError("Projection changed original coordinates or identity")
    omitted = full_map.keys() - projection_map.keys()
    links = _parse_explicit_links(full_lines, full)
    crossing = [(a, b) for a, b in sorted(links) if (a in projection_map) != (b in projection_map)]
    retained = {(a, b) for a, b in links if a in projection_map and b in projection_map}
    projected = _parse_explicit_links(projected_lines, projection)
    lost, added = retained - projected, projected - retained
    # Preserve every raw annotation, including those not resolved by the canonical
    # reader. A zero resolved-pair count is never claimed as a complete graph audit.
    raw = [line for line in full_lines if line.startswith(("LINK", "SSBOND", "CONECT"))]
    return {"original_full_atoms": len(full), "original_projection_atoms": len(projection),
            "omitted_atoms": len(omitted), "resolved_full_explicit_pairs": len(links),
            "resolved_full_pairs_between_retained_atoms": len(retained),
            "resolved_projection_pairs": len(projected),
            "lost_resolved_pairs_between_retained_atoms": [{"first": asdict(a), "second": asdict(b)} for a, b in sorted(lost)],
            "new_resolved_projection_pairs": [{"first": asdict(a), "second": asdict(b)} for a, b in sorted(added)],
            "resolved_retained_pair_custody_verified": not lost and not added,
            "resolved_projection_boundary_pairs": [{"first": asdict(a), "second": asdict(b),
                "first_retained": a in projection_map, "second_retained": b in projection_map,
                "first_residue_name": full_map[a].residue_name, "second_residue_name": full_map[b].residue_name}
                for a, b in crossing],
            "raw_full_connection_records": raw,
            "raw_full_connection_record_counts": dict(Counter(line[:6].strip() for line in raw)),
            "resolved_pair_chemical_classification": "not established from PDB LINK/SSBOND/CONECT alone",
            "unresolved_records_and_missing_annotations_not_qualified": True,
            "scientific_admission": False}


def audit(context_path, construct_path, output):
    context, construct = checked(context_path), checked(construct_path)
    contexts, constructs = [unique(data["rows"], "entry_id") for data in (context, construct)]
    if contexts.keys() != constructs.keys():
        raise ValueError("Entry denominators differ")
    bindings = {str(p): digest(p) for p in (context_path, construct_path)}
    reports, rows = [], []
    for entry, row in contexts.items():
        paths = [root.parent / source["report"] for root, source in
                 ((context_path, row), (construct_path, constructs[entry]))]
        for root, data, path in zip((context_path, construct_path), (context, construct), paths, strict=True):
            relative = str(path.relative_to(root.parent))
            if not path.resolve().is_relative_to(root.parent.resolve()) or digest(path) != data["reports_sha256"][relative]:
                raise ValueError("Per-entry report escaped or changed")
        ctx, cons = [json.loads(path.read_text()) for path in paths]
        if ctx["entry_id"] != entry or cons["entry_id"] != entry or ctx["scope"] != cons["scope"]:
            raise ValueError("Entry identity or scope changed")
        source, full, projection = [Path(name) for name in
            (cons["original_source_path"], ctx["original_full_context_path"], ctx["original_projection_path"])]
        if digest(source) != cons["original_source_sha256"]:
            raise ValueError("Original deposited source changed")
        for path in (full, projection):
            if str(path) not in context["bindings"] or digest(path) != context["bindings"][str(path)]:
                raise ValueError("Unbound original assembly/projection")
        text = full.read_text()
        if "REMARK 999 DVBfixer BIOLOGICAL ASSEMBLY 1" not in text.splitlines():
            raise ValueError("Original context does not declare the expected assembly")
        deposited = deposited_facts(gemmi.cif.read(str(source)).sole_block())
        actual = projection_link_facts(text, projection.read_text())
        if actual["omitted_atoms"] != ctx["omitted_atoms"]:
            raise ValueError("Actual omitted-atom denominator changed")
        report = {"entry_id": entry, "scope": ctx["scope"], "deposited": deposited, "actual_projection": actual,
                  "assembly_context_or_parameterization_qualified": False, "scientific_admission": False}
        reports.append(report)
        rows.append({"entry_id": entry, "scope": ctx["scope"],
                     "source_cross_boundary_annotations": sum(deposited["source_cross_boundary_types"].values()),
                     "assembly_member_label_cross_boundary_annotations": sum(deposited["assembly_member_label_cross_boundary_types"].values()),
                     "outside_assembly_member_label_cross_boundary_annotations": sum(deposited["outside_assembly_member_label_cross_boundary_types"].values()),
                     "actual_resolved_projection_boundary_pairs": len(actual["resolved_projection_boundary_pairs"]),
                     "resolved_pairs_between_retained_atoms": actual["resolved_full_pairs_between_retained_atoms"],
                     "lost_resolved_pairs_between_retained_atoms": len(actual["lost_resolved_pairs_between_retained_atoms"]),
                     "new_resolved_projection_pairs": len(actual["new_resolved_projection_pairs"]),
                     "raw_full_connection_records": len(actual["raw_full_connection_records"]),
                     "scientific_admission": False, "report": f"entries/{entry}.json"})
        for path in (*paths, source, full, projection):
            bindings[str(path)] = digest(path)
    summary = {"entries": len(rows),
               **{key: sum(row[key] for row in rows) for key in
                  ("source_cross_boundary_annotations", "assembly_member_label_cross_boundary_annotations",
                   "outside_assembly_member_label_cross_boundary_annotations", "actual_resolved_projection_boundary_pairs",
                   "resolved_pairs_between_retained_atoms", "lost_resolved_pairs_between_retained_atoms", "new_resolved_projection_pairs")},
               "source_cross_boundary_entries": sum(row["source_cross_boundary_annotations"] > 0 for row in rows),
               "scientific_admissions_or_new_sources_or_restoration_runs": 0,
               "missing_or_unresolved_annotations_qualified": False,
               "operator_pair_transport_or_forcefield_coverage_qualified": False}
    output.mkdir(parents=True, exist_ok=False)
    (output / "entries").mkdir()
    for row, report in zip(rows, reports, strict=True):
        save(output / row["report"], report)
    with (output / "connection-review.tsv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    save(output / "connection-audit.json", {"summary": summary, "rows": rows, "bindings": bindings,
         "auditor_sha256": digest(Path(__file__)),
         "reports_sha256": {row["report"]: digest(output / row["report"]) for row in rows},
         "TSV_sha256": digest(output / "connection-review.tsv")})
    save(output / "summary.json", summary)
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-audit", type=Path, required=True)
    parser.add_argument("--construct-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.context_audit.resolve(), args.construct_audit.resolve(), args.output.resolve())
