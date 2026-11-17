#!/usr/bin/env python3
"""Retain current RCSB sequence groups and whole-entry connected components."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from scripts.collect_real_gap_candidates import GRAPHQL, digest, post, save

QUERY = """
query NaturalGapGroups($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity_group_membership {
        aggregation_method similarity_cutoff group_id
      }
    }
  }
}
"""


def components(rows):
    """Conservative shared-protein30%-cluster OR exact-sequence connected groups."""
    owners = {}
    parent = {r["entry_id"]: r["entry_id"] for r in rows}

    def find(value):
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    for row in rows:
        for key in row["leakage_keys"]:
            if key in owners:
                first, second = find(row["entry_id"]), find(owners[key])
                parent[max(first, second)] = min(first, second)
            else:
                owners[key] = row["entry_id"]
    groups = defaultdict(list)
    for entry_id in parent:
        groups[find(entry_id)].append(entry_id)
    return [sorted(ids) for _, ids in sorted(groups.items())]


def run(census, output):
    complete = json.loads((census / "completed.json").read_text())
    data = (census / "census.json").read_bytes()
    if digest(data) != complete["census_sha256"]:
        raise ValueError("candidate census changed")
    rows = json.loads(data)["rows"]
    ids = [r["entry_id"] for r in rows]
    if len(ids) != len(set(ids)) or len(ids) != complete["entries"]:
        raise ValueError("candidate inventory incomplete/repeated")
    output.mkdir(parents=True, exist_ok=False)
    (output / "http").mkdir()
    save(output / "protocol.json", {"census_sha256": digest(data), "count": len(ids),
        "collector_sha256": digest(Path(__file__).read_bytes()),
        "HTTP_helper_sha256": digest(Path(post.__code__.co_filename).read_bytes()),
        "purpose": "modern RCSB group membership; deprecated cluster_membership returned null",
        "group_rule": "any shared protein30pct group or exact canonical protein sequence; transitive entry components",
        "final_splits_frozen": False, "sampling": False})
    entities = {}
    for batch, start in enumerate(range(0, len(ids), 20)):
        wanted = ids[start:start+20]
        response = post(output, f"group-{batch:03d}", GRAPHQL, {"query": QUERY, "variables": {"ids": wanted}})
        entries = response["data"]["entries"]
        if len(entries) != len(wanted) or {e["rcsb_id"] for e in entries if e is not None} != set(wanted):
            raise ValueError("RCSB grouping entries missing/repeated")
        for entry in entries:
            for entity in entry["polymer_entities"] or []:
                entities[entity["rcsb_id"]] = entity.get("rcsb_polymer_entity_group_membership") or []
        print(json.dumps({"grouping_entries_read": min(start+20, len(ids)), "total": len(ids)}), flush=True)
    enriched = []
    for row in rows:
        keys, unavailable = set(), set()
        for chain in row.get("chains", []):
            memberships = entities.get(chain["entity_id"], [])
            groups = [m["group_id"] for m in memberships
                      if m["aggregation_method"] == "sequence_identity" and m["similarity_cutoff"] == 30]
            keys.update(f"rcsb30:{g}" for g in groups)
            keys.add(f"exact-sequence:{chain['sequence_sha256']}")
            if not groups:
                unavailable.add(chain["entity_id"])
        enriched.append({**row, "leakage_keys": sorted(keys),
            "protein_entities_without_30pct_membership": sorted(unavailable),
            "modern_30pct_grouping_complete": not unavailable and bool(row.get("chains"))})
    groups = components(enriched)
    membership = {entry_id: index for index, ids_in_group in enumerate(groups) for entry_id in ids_in_group}
    for row in enriched:
        row["conservative_entry_component"] = membership[row["entry_id"]]
    save(output / "grouped-census.json", {"rows": enriched, "components": groups,
        "census_sha256": digest(data), "final_split_assignment": None,
        "limitation": "sequence clusters are grouping controls, not proof of independent folds or model-training absence"})
    summary = {"entries": len(enriched), "components": len(groups), "largest_component": max(map(len, groups)),
        "entries_with_complete_30pct_groups": sum(r["modern_30pct_grouping_complete"] for r in enriched),
        "grouped_census_sha256": digest((output / "grouped-census.json").read_bytes()),
        "scientific_admission_complete": False, "final_cohort_selected": False}
    save(output / "completed.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.census.resolve(), args.output.resolve())
