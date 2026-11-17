#!/usr/bin/env python3
"""Read deposited connection classes; never infer covalency from LINK distance."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import gemmi
from scripts.mac_reconstruction_comparison import digest, save


def declared_connections(block, protein_labels):
    columns = ["id", "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id",
               "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id"]
    rows = []
    for raw in block.find("_struct_conn.", columns):
        r = [gemmi.cif.as_string(v) for v in raw]
        first, second = r[2] in protein_labels, r[5] in protein_labels
        rows.append({"id": r[0], "depositor_connection_type": r[1],
            "first": {"label_chain": r[2], "component": r[3], "label_seq_id": r[4]},
            "second": {"label_chain": r[5], "component": r[6], "label_seq_id": r[7]},
            "crosses_protein_nonprotein_boundary": first != second,
            "interpretation": "deposited annotation only; chemistry/force-field routing not established"})
    return rows


def run(sources, output):
    source_data = sources / "source-results.json"
    completion = json.loads((sources / "completed.json").read_text())
    if digest(source_data) != completion["source_results_sha256"]:
        raise ValueError("source ledger changed")
    rows = []
    for row in json.loads(source_data.read_text())["rows"]:
        source = Path(row["source_path"])
        if digest(source) != row["source_sha256"]:
            raise ValueError("original mmCIF changed")
        document = gemmi.cif.read(str(source))
        block = document.sole_block()
        protein_labels = {c["label_chain"] for c in row["inventory"]["chains"]}
        connections = declared_connections(block, protein_labels)
        cross = Counter(c["depositor_connection_type"] for c in connections if c["crosses_protein_nonprotein_boundary"])
        rows.append({"entry_id": row["entry_id"], "source_sha256": row["source_sha256"],
            "connections": connections, "protein_nonprotein_connection_types": dict(cross),
            "protein_nonprotein_metal_coordination_annotated": cross.get("metalc", 0),
            "other_protein_nonprotein_connections_need_review": sum(n for kind, n in cross.items() if kind != "metalc"),
            "no_annotation_does_not_prove_no_bond": True, "scientific_admission_complete": False,
            "context_removed_or_modified": False})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "context-annotations.json", {"rows": rows, "source_results_sha256": digest(source_data),
        "classifier_sha256": digest(Path(__file__)), "not_a_parameterization_or_admission_decision": True})
    summary = {"entries": len(rows), "with_annotated_protein_metal_coordination": sum(bool(r["protein_nonprotein_metal_coordination_annotated"]) for r in rows),
        "with_other_protein_nonprotein_links": sum(bool(r["other_protein_nonprotein_connections_need_review"]) for r in rows),
        "deposited_cross_boundary_types": dict(sum((Counter(r["protein_nonprotein_connection_types"]) for r in rows), Counter())),
        "source_coordinates_modified": False, "forcefield_coverage_established": False}
    save(output / "completed.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.sources.resolve(), args.output.resolve())
