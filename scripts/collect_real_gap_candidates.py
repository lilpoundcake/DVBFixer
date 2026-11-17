#!/usr/bin/env python3
"""Immutable public-metadata census of natural gaps, NOT scientific admission.

No coordinates are hidden, edited or generated. RCSB label-sequence features are
retained as metadata evidence, pending deposited mmCIF/assembly/sequence audits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
SALT = "dvbfixer-real-natural-gap-census-2026-11-12-v1"
PINNED = ("1HE8", "1AHW", "1N8Z", "1G9M", "1GC1", "1JPS", "1I2M", "9CQD", "1KTZ", "1BRS")
BINS = ((1, 8), (9, 24), (25, 49), (50, 99), (100, None))
QUERY = """
query NaturalGapMetadata($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info {
      deposited_modeled_polymer_monomer_count
      deposited_unmodeled_polymer_monomer_count
    }
    assemblies { rcsb_id }
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can type }
      rcsb_cluster_membership { cluster_id identity }
      polymer_entity_instances {
        rcsb_id
        rcsb_polymer_entity_instance_container_identifiers { asym_id auth_asym_id }
        rcsb_polymer_instance_feature {
          type
          feature_positions { beg_seq_id end_seq_id }
        }
      }
    }
  }
}
"""


def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    with path.open("xb") as stream:
        stream.write(encoded(value))


def search_query(lower, upper):
    count = "rcsb_entry_info.deposited_unmodeled_polymer_monomer_count"
    nodes = [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": count, "operator": "greater_or_equal", "value": lower}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.deposited_modeled_polymer_monomer_count",
            "operator": "greater_or_equal", "value": 300}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "Protein"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "exptl.method", "operator": "in", "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]}},
    ]
    if upper is not None:
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": count, "operator": "less_or_equal", "value": upper}})
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry", "request_options": {"return_all_hits": True,
            "results_verbosity": "compact", "results_content_type": ["experimental"]}}


def post(root, name, url, payload):
    """Manifest-last receipts; never adopt a request/response without its receipt."""
    prefix = root / "http" / name
    request_path = prefix.with_suffix(".request.json")
    response_path = prefix.with_suffix(".response.json")
    receipt_path = prefix.with_suffix(".receipt.json")
    request_bytes = encoded(payload)
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_bytes())
        response = response_path.read_bytes()
        if (receipt["url"] != url or receipt["request_sha256"] != digest(request_bytes)
                or request_path.read_bytes() != request_bytes or receipt["response_sha256"] != digest(response)):
            raise ValueError("cached public-metadata receipt changed")
        return json.loads(response)
    if request_path.exists() or response_path.exists():
        raise ValueError("uncommitted metadata request/response retained; explicit recovery required")
    if shutil.disk_usage(root).free < 2 * 1024**3:
        raise RuntimeError("public-metadata acquisition lacks 2 GiB disk headroom")
    if sum(p.stat().st_size for p in root.rglob("*") if p.is_file()) > 512 * 1024**2:
        raise RuntimeError("512 MiB public-metadata evidence cap exceeded")
    with request_path.open("xb") as stream:
        stream.write(request_bytes)
    errors = []
    for attempt in range(1, 4):
        req = Request(url, data=request_bytes, headers={"Content-Type": "application/json",
                      "User-Agent": "dvbfixer-natural-gap-research/1"}, method="POST")
        try:
            with urlopen(req, timeout=60) as response:
                body = response.read(32 * 1024**2 + 1)
                status = response.status
            if len(body) > 32 * 1024**2:
                raise ValueError("public-metadata response exceeds 32 MiB cap")
            with response_path.open("xb") as stream:
                stream.write(body)
            parsed = json.loads(body)
            if not isinstance(parsed, dict) or parsed.get("errors"):
                raise ValueError(f"public service returned invalid/schema-error metadata: {name}")
            save(receipt_path, {"url": url, "request_sha256": digest(request_bytes),
                "response_sha256": digest(body), "status": status, "attempts": attempt,
                "received_utc": datetime.now(UTC).isoformat(), "previous_errors": errors})
            return parsed
        except (HTTPError, URLError, TimeoutError) as exc:
            failure = {"attempt": attempt, "error_type": type(exc).__name__, "message": str(exc)}
            errors.append(failure)
            save(root / "http" / f"{name}.attempt{attempt}.failure.json", failure)
    raise RuntimeError(f"public acquisition failed after three retained attempts: {name}")


def identifiers(response):
    hits = response.get("result_set", [])
    ids = [hit if isinstance(hit, str) else hit["identifier"] for hit in hits]
    if len(ids) != len(set(ids)) or len(ids) != response.get("total_count"):
        raise ValueError("public search returned duplicate/incomplete all-hit inventory")
    if any(not isinstance(i, str) or not i.isalnum() for i in ids):
        raise ValueError("public search returned malformed entry IDs")
    return ids


def spans(features, length):
    intervals = []
    for feature in features or []:
        if feature["type"] != "UNOBSERVED_RESIDUE_XYZ":
            continue
        for position in feature.get("feature_positions") or []:
            start = position["beg_seq_id"]
            stop = position.get("end_seq_id") or start
            if (type(start) is not int or type(stop) is not int
                    or not 1 <= start <= stop <= length):
                raise ValueError("unobserved-residue feature outside label sequence")
            intervals.append((start - 1, stop))
    merged = []
    for start, stop in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(stop, merged[-1][1])
        else:
            merged.append((start, stop))
    return [{"start": a, "stop": b, "length": b-a,
             "kind": "whole-chain" if a == 0 and b == length else
                     "N-terminal" if a == 0 else "C-terminal" if b == length else "internal"}
            for a, b in merged]


def summarize(entry):
    chains = []
    for entity in entry.get("polymer_entities") or []:
        poly = entity.get("entity_poly") or {}
        if poly.get("type") != "polypeptide(L)":
            continue
        sequence = "".join((poly.get("pdbx_seq_one_letter_code_can") or "").split())
        if not sequence:
            raise ValueError("protein metadata has no canonical sequence")
        memberships = entity.get("rcsb_cluster_membership") or []
        groups = [m["cluster_id"] for m in memberships if m["identity"] == 30]
        for instance in entity.get("polymer_entity_instances") or []:
            ids = instance["rcsb_polymer_entity_instance_container_identifiers"]
            gaps = spans(instance.get("rcsb_polymer_instance_feature"), len(sequence))
            chains.append({"instance_id": instance["rcsb_id"], "entity_id": entity["rcsb_id"],
                "label_chain": ids["asym_id"], "author_chain": ids["auth_asym_id"],
                "sequence": sequence, "sequence_sha256": digest(sequence.encode()),
                "noncanonical_sequence_symbols": sorted(set(sequence)-set("ACDEFGHIKLMNPQRSTVWY")),
                "target_residues": len(sequence), "gaps": gaps,
                "feature_inferred_observed_residues": len(sequence)-sum(g["length"] for g in gaps),
                "rcsb_30pct_clusters": groups, "cluster_assignment_available": bool(groups)})
    gaps = [g for c in chains for g in c["gaps"]]
    entry_id = entry["rcsb_id"]
    return {"entry_id": entry_id, "title": (entry.get("struct") or {}).get("title"),
        "methods": [m["method"] for m in entry.get("exptl") or []],
        "assembly_ids": [a["rcsb_id"] for a in entry.get("assemblies") or []],
        "chains": chains, "deposited_entry_info": entry.get("rcsb_entry_info"),
        "deposited_protein_target_residues": sum(c["target_residues"] for c in chains),
        "deposited_protein_observed_residues_feature_estimate": sum(c["feature_inferred_observed_residues"] for c in chains),
        "deposited_protein_chain_count": len(chains), "natural_gap_count_metadata": len(gaps),
        "longest_gap_metadata": max((g["length"] for g in gaps), default=0),
        "gap_kinds_metadata": dict(Counter(g["kind"] for g in gaps)),
        "previously_inspected": entry_id in PINNED, "small_control_not_primary": entry_id == "1BRS",
        "evidence_scope": "RCSB deposited-instance annotation, not source-coordinate/biological-assembly admission",
        "scientific_admission_complete": False, "split_assignment": None,
        "needs_source_assembly_sequence_occupancy_chemistry_audit": True,
        "needs_family_grouping_before_final_split": True, "checkpoint_sampling_started": False}


def collect(output, count=1000, resume=False):
    if not 10 <= count <= 1000:
        raise ValueError("declared public metadata count must be 10..1000")
    if resume:
        protocol = json.loads((output / "protocol.json").read_bytes())
        if protocol["collector_sha256"] != digest(Path(__file__).read_bytes()) or protocol["count"] != count:
            raise ValueError("immutable collector/target changed; use an explicit new attempt")
        if (output / "completed.json").exists():
            raise ValueError("completed census is immutable")
    else:
        output.mkdir(parents=True, exist_ok=False)
        (output / "http").mkdir()
        protocol = {"protocol": "natural-gap-public-metadata-census-v1", "count": count,
            "created_utc": datetime.now(UTC).isoformat(), "selection_salt": SALT,
            "collector_sha256": digest(Path(__file__).read_bytes()), "pinned_inspected_entries": PINNED,
            "search_strata": BINS, "strata_are_entry_total_missing_proxy_not_longest_gap": True,
            "primary_target": 200, "metadata_only": True, "checkpoint_sampling": False,
            "public_services": [SEARCH, GRAPHQL], "sequence_upload_to_MSA_service": False,
            "request_timeout_seconds": 60, "maximum_HTTP_attempts": 3,
            "maximum_response_bytes": 32*1024**2, "maximum_evidence_bytes": 512*1024**2,
            "min_free_disk_bytes": 2*1024**3, "download_coordinates": False}
        save(output / "protocol.json", protocol)
        (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    selected = list(PINNED)
    populations = []
    share, remainder = divmod(count-len(selected), len(BINS))
    for index, (lower, upper) in enumerate(BINS):
        response = post(output, f"search-{index}", SEARCH, search_query(lower, upper))
        ids = identifiers(response)
        ordered = sorted(ids, key=lambda i: digest(f"{SALT}:{index}:{i}".encode()))
        candidates = [i for i in ordered if i not in selected]
        needed = share + (index < remainder)
        if len(candidates) < needed:
            raise ValueError("public proxy stratum lacks the declared candidate count")
        selected.extend(candidates[:needed])
        populations.append({"entry_total_missing_proxy": [lower, upper], "population": len(ids), "selected_new": needed})
    selection = {"entry_ids": selected, "count": len(selected), "populations": populations,
                 "not_a_final_cohort_or_split": True}
    if (output / "selection.json").exists():
        if (output / "selection.json").read_bytes() != encoded(selection):
            raise ValueError("frozen metadata selection differs")
    else:
        save(output / "selection.json", selection)
    rows = []
    for batch, start in enumerate(range(0, len(selected), 12)):
        ids = selected[start:start+12]
        response = post(output, f"metadata-{batch:03d}", GRAPHQL, {"query": QUERY, "variables": {"ids": ids}})
        entries = response["data"]["entries"]
        found = {e["rcsb_id"]: e for e in entries if e is not None}
        if len(entries) != len(ids) or set(found) != set(ids):
            raise ValueError("public entry metadata missing, repeated or undeclared")
        for entry_id in ids:
            try:
                row = summarize(found[entry_id])
            except (ValueError, KeyError, TypeError) as exc:
                row = {"entry_id": entry_id, "status": "metadata_inventory_refused", "error": str(exc),
                       "scientific_admission_complete": False, "checkpoint_sampling_started": False}
            rows.append(row)
        print(json.dumps({"metadata_entries_read": len(rows), "declared": len(selected)}), flush=True)
    save(output / "census.json", {"rows": rows, "selection_sha256": digest(encoded(selection)),
        "qualification": False, "source_admission_not_performed": True})
    summary = {"entries": len(rows), "metadata_refusals": sum(r.get("status") == "metadata_inventory_refused" for r in rows),
        "with_reported_natural_gaps": sum(r.get("natural_gap_count_metadata", 0) > 0 for r in rows),
        "longest_gap_bins": {f"{lo}-{hi or 'plus'}": sum(lo <= r.get("longest_gap_metadata", 0) and
            (hi is None or r.get("longest_gap_metadata", 0) <= hi) for r in rows) for lo, hi in BINS},
        "scientific_admission_complete": False, "final_200_selected": False, "checkpoint_sampling_started": False}
    save(output / "completed.json", {**summary, "census_sha256": digest((output / "census.json").read_bytes())})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    collect(args.output.resolve(), args.count, args.resume)
