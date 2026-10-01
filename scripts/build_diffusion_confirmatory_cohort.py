#!/usr/bin/env python3
"""Lock a nonredundant post-release diffusion benchmark candidate cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

from dvbfixer.model.diffusion.benchmark import paired_noninferiority_power

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
FINAL_MINIMUM_GROUPS = 180
FINAL_MAXIMUM_GROUPS = 300
SCREENING_POOL_GROUPS = 500
DEFAULT_GROUPS = SCREENING_POOL_GROUPS
DEFAULT_CANDIDATES = 800
MINIMUM_COMPLETENESS = 0.90
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
SELECTION_SEED = "dvbfixer-diffusion-confirmatory-v1"

ENTITY_QUERY = """
query EntityMetadata($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    entity_poly {
      rcsb_sample_sequence_length
      pdbx_seq_one_letter_code_can
    }
    entry {
      rcsb_id
      rcsb_accession_info { initial_release_date }
      exptl { method }
      rcsb_entry_info { resolution_combined }
    }
    polymer_entity_instances {
      rcsb_id
      rcsb_polymer_entity_instance_container_identifiers {
        asym_id
        auth_asym_id
        auth_to_entity_poly_seq_mapping
      }
    }
  }
}
"""


def search_request(start_date: str, end_date: str, rows: int) -> dict[str, Any]:
    """Return the frozen RCSB grouped-search request."""
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "range",
                        "value": {
                            "from": start_date,
                            "to": end_date,
                            "include_lower": True,
                            "include_upper": True,
                        },
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "range",
                        "value": {
                            "from": 80,
                            "to": 500,
                            "include_lower": True,
                            "include_upper": True,
                        },
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "in",
                        "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": 3.5,
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "results_verbosity": "compact",
            "group_by": {
                "aggregation_method": "sequence_identity",
                "similarity_cutoff": 30,
                "ranking_criteria_type": {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                },
            },
            "group_by_return_type": "representatives",
        },
    }


def post_json(url: str, payload: dict[str, Any]) -> dict[str, Any]:
    request = Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "dvbfixer-diffusion-cohort/1",
        },
        method="POST",
    )
    with urlopen(request, timeout=120) as response:
        result = json.load(response)
    if not isinstance(result, dict):
        raise RuntimeError(f"{url} returned a non-object JSON response")
    return result


def _canonical_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _entity_ids(search_response: dict[str, Any]) -> list[str]:
    identifiers: list[str] = []
    for result in search_response.get("result_set", []):
        identifier = result if isinstance(result, str) else result.get("identifier", "")
        if not isinstance(identifier, str) or "_" not in identifier:
            raise RuntimeError("RCSB grouped search returned an invalid polymer-entity ID")
        identifiers.append(identifier.upper())
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError("RCSB grouped search repeated a representative")
    return identifiers


def _sequence(value: Any) -> str:
    if not isinstance(value, str):
        return ""
    return "".join(value.split()).upper()


def _modeled_count(mapping: Any) -> int:
    if not isinstance(mapping, list):
        return 0
    return sum(value not in {None, "?", "."} for value in mapping)


def _record(entity: dict[str, Any]) -> dict[str, Any] | None:
    entity_id = entity.get("rcsb_id")
    entity_poly = entity.get("entity_poly") or {}
    entry = entity.get("entry") or {}
    sequence = _sequence(entity_poly.get("pdbx_seq_one_letter_code_can"))
    sequence_length = entity_poly.get("rcsb_sample_sequence_length")
    if not isinstance(entity_id, str) or not isinstance(sequence_length, int):
        return None
    if len(sequence) != sequence_length or not set(sequence) <= CANONICAL_AMINO_ACIDS:
        return None

    chains: list[dict[str, Any]] = []
    for instance in entity.get("polymer_entity_instances") or []:
        identifiers = instance.get("rcsb_polymer_entity_instance_container_identifiers") or {}
        label_chain = identifiers.get("asym_id")
        author_chain = identifiers.get("auth_asym_id")
        modeled = _modeled_count(identifiers.get("auth_to_entity_poly_seq_mapping"))
        if not isinstance(label_chain, str) or not label_chain:
            continue
        chains.append(
            {
                "label_asym_id": label_chain,
                "auth_asym_id": author_chain if isinstance(author_chain, str) else "",
                "modeled_residues": modeled,
            }
        )
    if not chains:
        return None
    chain = min(chains, key=lambda item: (-item["modeled_residues"], item["label_asym_id"]))
    completeness = chain["modeled_residues"] / sequence_length
    if completeness < MINIMUM_COMPLETENESS:
        return None

    entry_id = entry.get("rcsb_id")
    accession = entry.get("rcsb_accession_info") or {}
    release_date = accession.get("initial_release_date")
    resolutions = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    methods = sorted(
        method["method"]
        for method in entry.get("exptl") or []
        if isinstance(method, dict) and isinstance(method.get("method"), str)
    )
    if (
        not isinstance(entry_id, str)
        or not isinstance(release_date, str)
        or not resolutions
        or not methods
    ):
        return None
    resolution = min(float(value) for value in resolutions)
    return {
        "independence_group": f"rcsb-30pct-representative:{entity_id.upper()}",
        "polymer_entity_id": entity_id.upper(),
        "pdb_id": entry_id.upper(),
        "entity_id": entity_id.rsplit("_", 1)[-1],
        "label_asym_id": chain["label_asym_id"],
        "auth_asym_id": chain["auth_asym_id"],
        "method": methods,
        "resolution_angstrom": resolution,
        "initial_release_date": release_date[:10],
        "sequence": sequence,
        "sequence_length": sequence_length,
        "modeled_residues": chain["modeled_residues"],
        "completeness": round(completeness, 6),
        "mmcif_url": f"https://files.rcsb.org/download/{entry_id.upper()}.cif",
    }


def build_manifest(
    *,
    start_date: str,
    end_date: str,
    group_count: int = DEFAULT_GROUPS,
    candidate_count: int = DEFAULT_CANDIDATES,
    request_json: Callable[[str, dict[str, Any]], dict[str, Any]] = post_json,
) -> dict[str, Any]:
    """Query RCSB and return a locked candidate-cohort manifest."""
    if group_count != SCREENING_POOL_GROUPS:
        raise ValueError(
            f"blind coordinate screening must lock exactly {SCREENING_POOL_GROUPS} metadata "
            f"groups; the final inferential sample is limited to "
            f"[{FINAL_MINIMUM_GROUPS}, {FINAL_MAXIMUM_GROUPS}] accepted groups"
        )
    if candidate_count < group_count:
        raise ValueError("candidate count must not be below the requested group count")

    query = search_request(start_date, end_date, candidate_count)
    search_response = request_json(SEARCH_URL, query)
    entity_ids = _entity_ids(search_response)
    if len(entity_ids) < group_count:
        raise RuntimeError("RCSB returned fewer sequence-cluster representatives than requested")

    entities: dict[str, dict[str, Any]] = {}
    for start in range(0, len(entity_ids), 50):
        batch = entity_ids[start : start + 50]
        response = request_json(
            GRAPHQL_URL,
            {"query": ENTITY_QUERY, "variables": {"ids": batch}},
        )
        if response.get("errors"):
            raise RuntimeError(f"RCSB GraphQL error: {response['errors']}")
        for entity in (response.get("data") or {}).get("polymer_entities") or []:
            if isinstance(entity, dict) and isinstance(entity.get("rcsb_id"), str):
                entities[entity["rcsb_id"].upper()] = entity

    candidates = [
        record
        for entity_id in entity_ids
        if (record := _record(entities.get(entity_id, {}))) is not None
    ]
    candidates.sort(
        key=lambda item: hashlib.sha256(
            f"{SELECTION_SEED}:{item['polymer_entity_id']}".encode()
        ).hexdigest()
    )
    selected: list[dict[str, Any]] = []
    selected_entries: set[str] = set()
    for candidate in candidates:
        if candidate["pdb_id"] in selected_entries:
            continue
        candidate["planned_gap_length"] = 5 if len(selected) % 2 == 0 else 10
        candidate["request_sequence_basis"] = "observed-coordinate-sequence"
        candidate["screening_index"] = len(selected)
        selected.append(candidate)
        selected_entries.add(candidate["pdb_id"])
        if len(selected) == group_count:
            break
    if len(selected) < group_count:
        raise RuntimeError(
            f"only {len(selected)} unique-entry canonical representatives passed metadata filters"
        )

    return {
        "schema_version": 1,
        "status": "candidate-metadata-locked-coordinate-screening-pending",
        "materialized_on": end_date,
        "statistical_unit": "one-withheld-gap-per-sequence-cluster",
        "screening_pool_count": len(selected),
        "final_minimum_independence_groups": FINAL_MINIMUM_GROUPS,
        "final_maximum_independence_groups": FINAL_MAXIMUM_GROUPS,
        "request_sequence_basis": "observed-coordinate-sequence",
        "entity_sequence_role": "rcsb-clustering-and-provenance-only",
        "screening_order_locked": True,
        "screening_order": "cases-array-order-after-selection-seed-ranking",
        "screen_all_metadata_groups_before_inference": True,
        "sequence_identity_cluster_percent": 30,
        "minimum_chain_completeness": MINIMUM_COMPLETENESS,
        "selection_seed": SELECTION_SEED,
        "release_window": {"from": start_date, "to": end_date},
        "search_api": SEARCH_URL,
        "metadata_api": GRAPHQL_URL,
        "search_request": query,
        "search_response_sha256": _canonical_digest(search_response),
        "matching_polymer_entity_count": search_response.get("total_count"),
        "matching_sequence_cluster_count": search_response.get("group_by_count"),
        "returned_representative_count": len(entity_ids),
        "eligible_backends": ["protenix-v1"],
        "proxy_only_backends": ["boltz-2"],
        "excluded_backends": ["rfdiffusion-v1"],
        "boltz_leakage_note": (
            "Post-release PDB timing is only a proxy. Boltz-2 remains ineligible because "
            "its non-PDB training membership and sequence homology cannot be audited."
        ),
        "coordinate_screening_required": True,
        "power_model": {
            "method": "paired-normal-approximation-for-two-sided-ci-noninferiority",
            "confidence_level": 0.95,
            "noninferiority_margin": 0.05,
            "assumed_true_pass_rate_difference": 0.0,
            "assumed_discordance_rate": 0.05,
            "estimated_power": paired_noninferiority_power(
                FINAL_MINIMUM_GROUPS,
                discordance_rate=0.05,
            ),
            "power_sample_size": FINAL_MINIMUM_GROUPS,
            "final_inference_maximum_groups": FINAL_MAXIMUM_GROUPS,
        },
        "cases": selected,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--start-date", default="2025-06-07")
    parser.add_argument("--end-date", default="2026-09-25")
    parser.add_argument("--groups", type=int, default=DEFAULT_GROUPS)
    parser.add_argument("--candidates", type=int, default=DEFAULT_CANDIDATES)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = build_manifest(
        start_date=args.start_date,
        end_date=args.end_date,
        group_count=args.groups,
        candidate_count=args.candidates,
    )
    args.output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
