#!/usr/bin/env python3
"""Weightless shared-context planning on retained loops, not native preflight."""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import BackendOption, DiffusionRequest
from dvbfixer.model.diffusion.protein_context import (
    ProteinContextError,
    prepare_cc94_request,
    prepare_protenix_request,
    validate_context_capability,
    validate_context_source,
)

CASES = ("1NON", "2BE1")
ROUTES = ("protenix-full", "protenix-auto", "cc94-auto")


def profile_view(request, profile):
    return replace(request, backend_options=tuple(o for o in request.backend_options if o.name != "profile")
                   + (BackendOption("profile", profile),))


def run(scope, output):
    scope_path = scope / "loop-scope-inventory.json"
    data = json.loads(scope_path.read_text())
    inventory = {r["entry_id"]: r for r in data["rows"]}
    bindings = {str(scope_path): digest(scope_path)}
    cases = []
    for entry in CASES:
        row = inventory[entry]
        if (row["current_scope_status"] != "loop_scope_and_context_review_pending"
                or row["terminal_construct_review_required"] or row["small_assembly_control"]):
            raise ValueError("planning controls must be retained non-small, internal-gap-only requests")
        request_path = Path(row["request_path"])
        if digest(request_path) != row["request_sha256"]:
            raise ValueError("retained source request changed")
        # Find the source via the readiness result, not by guessing a nested path.
        result = json.loads((request_path.parent / "result.json").read_text())
        source = Path(result["source_input"])
        if digest(source) != result["source_input_sha256"]:
            raise ValueError("retained source bytes changed")
        projection_path = request_path.parent / "projection.json"
        projection = json.loads(projection_path.read_text())
        if set(projection["omitted_heterogen_atom_counts"]) - {"HOH", "WAT"}:
            raise ValueError("planning controls must not silently omit non-water chemical context")
        request = DiffusionRequest.from_json(request_path.read_text())
        if request.normalized_pdb.sha256 != digest(source):
            raise ValueError("request bound to different source")
        for path in (request_path, source, projection_path):
            bindings[str(path)] = digest(path)
        cases.append((entry, request, source.read_bytes(), request_path))
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"cases": CASES, "routes": ROUTES, "bindings": bindings,
        "implementation_sha256": digest(Path(__file__)), "checkpoint_loads": 0, "new_sampling_tasks": 0,
        "purpose": "shared Python context-contract feasibility only; not scientific cohort admission or pinned native consumer execution",
        "profiles_are_derived_views_not_edits_to_retained_requests": True,
        "full_vs_local_is_not_a_matched_method_comparison": True})
    rows = []
    for entry, original, source, parent in cases:
        for route in ROUTES:
            directory = output / entry / route
            directory.mkdir(parents=True, exist_ok=False)
            row = {"entry_id": entry, "route": route, "parent_request_sha256": digest(parent),
                   "target_residues": sum(len(t.sequence) for t in original.target_sequences),
                   "gaps": len(original.gaps), "native_consumer_preflight_executed": False,
                   "checkpoint_loaded": False, "scientific_admission": False}
            try:
                if route == "cc94-auto":
                    view = profile_view(original, "protpardelle-1c-cc94-mps-v1")
                    planned = prepare_cc94_request(view, source, "auto")
                    model = "protpardelle"
                else:
                    view = original
                    planned = prepare_protenix_request(view, source, "full" if route.endswith("full") else "auto")
                    model = "protenix"
                if (planned.fixed_atoms != original.fixed_atoms or planned.generated_atoms != original.generated_atoms
                        or planned.target_sequences != original.target_sequences or planned.gaps != original.gaps):
                    raise RuntimeError("context builder changed parent atom/gap ownership")
                validate_context_source(planned, source)
                validate_context_capability(planned, model)
                plans = planned.gap_contexts or ((planned.protein_context,) if planned.protein_context else ())
                if not plans or any(p is None for p in plans):
                    raise RuntimeError("successful context plan lacks explicit axes")
                (directory / "derived-context-request.json").write_text(planned.to_json())
                row.update(status="shared_context_contract_passed", context_axes=len(plans),
                    context_residues_per_axis=[sum(s.stop-s.start for s in p.spans) for p in plans],
                    coupled_gap_sampling_established=False,
                    derived_request_sha256=digest(directory / "derived-context-request.json"))
            except ProteinContextError as exc:
                row.update(status="shared_context_contract_refused", error=str(exc))
            save(directory / "result.json", row)
            rows.append(row)
            print(json.dumps(row), flush=True)
    save(output / "results.json", {"rows": rows, "no_source_or_frozen_request_mutations": True,
        "no_checkpoint_or_native_consumer_executed": True, "final_cohort_selected": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.scope.resolve(), args.output.resolve())
