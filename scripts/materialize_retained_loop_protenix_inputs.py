#!/usr/bin/env python3
"""Exercise the shipped multichain template encoder, without model/data loading."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import gemmi
import numpy as np
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.protein_context import sampling_entities, source_coordinates

PROJECT = Path(__file__).resolve().parents[1]
ENCODER = PROJECT / "deploy/protenix-v1/build-template-input.py"


def audit_template(payload, entity, coordinates):
    if len(payload) != 1:
        raise ValueError("template encoder produced ambiguous template entries")
    template = payload[0]
    if (template["queryIndices"] != list(entity.placement.observed_target_indices)
            or template["templateIndices"] != list(range(len(entity.placement.observed_residues)))):
        raise ValueError("template alignment differs from source-only context placement")
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(template["mmcif"]).sole_block())
    if len(structure) != 1 or len(structure[0]) != 1 or structure[0][0].name != entity.engine_chain:
        raise ValueError("encoder changed engine entity identity")
    observed = entity.placement.observed_residues
    residues = list(structure[0][0])
    if [(str(r.seqid.num), r.seqid.icode.strip()) for r in residues] != [(r.residue_number, r.insertion_code) for r in observed]:
        raise ValueError("encoder changed observed residue identity/order")
    seen = set()
    for residue in residues:
        for atom in residue:
            identity = AtomIdentity(entity.target.chain, str(residue.seqid.num), residue.seqid.icode.strip(), atom.name)
            if identity in seen or identity not in coordinates:
                raise ValueError("template contains duplicate or non-observed coordinates")
            seen.add(identity)
            if not np.array_equal(np.array([atom.pos.x, atom.pos.y, atom.pos.z]), coordinates[identity]):
                raise ValueError("template changed an observed coordinate")
    residue_set = set(observed)
    expected = {a for a in coordinates if any(a.chain == r.chain and a.residue_number == r.residue_number
                and a.insertion_code == r.insertion_code for r in residue_set)}
    if seen != expected:
        raise ValueError("template lost observed atoms including terminal OXT")
    return len(seen)


def materialize(planning, scope, output):
    planned = json.loads((planning / "results.json").read_text())
    source_rows = {r["entry_id"]: r for r in json.loads((scope / "loop-scope-inventory.json").read_text())["rows"]}
    selected = [r for r in planned["rows"] if r["route"] == "protenix-auto"
                and r["status"] == "shared_context_contract_passed"]
    spec = importlib.util.spec_from_file_location("retained_loop_template_encoder", ENCODER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load existing template encoder")
    encoder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(encoder)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"case_ids": [r["entry_id"] for r in selected],
        "planning_sha256": digest(planning / "results.json"), "scope_sha256": digest(scope / "loop-scope-inventory.json"),
        "encoder_sha256": digest(ENCODER), "auditor_sha256": digest(Path(__file__)),
        "model_weights_loaded": False, "native_dataloader_executed": False,
        "sampling_tasks_authorized": 0, "scientific_cohort_admission": False})
    rows = []
    for row in selected:
        entry = row["entry_id"]
        source_row = source_rows[entry]
        if source_row["current_scope_status"] != "loop_scope_and_context_review_pending":
            raise ValueError("scope excludes the retained planning control")
        parent_request = Path(source_row["request_path"])
        if digest(parent_request) != row["parent_request_sha256"]:
            raise ValueError("parent source request changed")
        ready = json.loads((parent_request.parent / "result.json").read_text())
        source = Path(ready["source_input"])
        derived_path = planning / entry / row["route"] / "derived-context-request.json"
        if digest(derived_path) != row["derived_request_sha256"] or digest(source) != ready["source_input_sha256"]:
            raise ValueError("derived request/source changed")
        request = DiffusionRequest.from_json(derived_path.read_text())
        work = output / entry
        work.mkdir()
        staged_source = (work / request.normalized_pdb.path).resolve()
        if not staged_source.is_relative_to(work.resolve()):
            raise ValueError("source artifact path escapes owned staging directory")
        staged_source.parent.mkdir(parents=True, exist_ok=True)
        with staged_source.open("xb") as stream:
            stream.write(source.read_bytes())
        staged_request = work / "context-request.json"
        with staged_request.open("xb") as stream:
            stream.write(derived_path.read_bytes())
        input_path = encoder.build_input(staged_request, work / "encoded")
        data = json.loads(input_path.read_text())
        entities = sampling_entities(request)
        sequences = data[0]["sequences"]
        if len(data) != 1 or len(sequences) != len(entities):
            raise ValueError("template encoder changed entity count")
        coordinates = source_coordinates(request, source.read_bytes())
        counts = []
        for record, entity in zip(sequences, entities, strict=True):
            protein = record["proteinChain"]
            if protein["sequence"] != entity.target.sequence or protein["id"] != [entity.engine_chain]:
                raise ValueError("template encoder changed native entity sequence/chain")
            template_path = Path(protein["templatesPath"])
            if not template_path.resolve().is_relative_to((work / "encoded").resolve()):
                raise ValueError("template path escapes owned encoder output")
            counts.append(audit_template(json.loads(template_path.read_text()), entity, coordinates))
        result = {"entry_id": entry, "engine_entities": len(entities),
            "context_residues": sum(len(e.target.sequence) for e in entities),
            "observed_template_atoms_per_entity": counts, "source_and_request_staged_byte_exact": True,
            "source_only_templates_independently_checked": True, "input_json_sha256": digest(input_path),
            "native_dataloader_executed": False, "model_weights_loaded": False,
            "sampling": False, "scientific_admission": False}
        save(work / "encoder-readback.json", result)
        rows.append(result)
        print(json.dumps(result), flush=True)
    save(output / "results.json", {"rows": rows, "not_native_feature_preflight_or_sampling": True})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--planning", type=Path, required=True)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    materialize(args.planning.resolve(), args.scope.resolve(), args.output.resolve())
