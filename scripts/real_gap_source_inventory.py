"""Deposited mmCIF label-axis inventory, not pipeline-specific CIF normalization.

This research reader never emits a PDB, reassigns IDs or repairs observations.
Scientific pipelines still use the shared CLI boundary and admission services.
"""
from __future__ import annotations

import math
from collections import defaultdict

import gemmi


def missing_intervals(length, observed):
    if any(type(i) is not int or not 1 <= i <= length for i in observed):
        raise ValueError("observed label_seq_id is outside the deposited sequence")
    missing = set(range(1, length+1))-set(observed)
    intervals = []
    for index in sorted(missing):
        if intervals and index == intervals[-1][1]+1:
            intervals[-1][1] = index
        else:
            intervals.append([index, index])
    return [{"start": a-1, "stop": b, "length": b-a+1,
             "kind": "whole-chain" if a == 1 and b == length else
             "N-terminal" if a == 1 else "C-terminal" if b == length else "internal"}
            for a, b in intervals]


def inventory_mmcif(text, metadata):
    document = gemmi.cif.read_string(text)
    if len(document) != 1:
        raise ValueError("source mmCIF must contain exactly one entry block")
    block = document.sole_block()

    def value(raw):
        return gemmi.cif.as_string(raw)

    entry_id = value(block.find_value("_entry.id") or "?")
    if entry_id != metadata["entry_id"]:
        raise ValueError("deposited source entry identity differs from metadata")
    entities = {}
    for row in block.find("_entity_poly.", ["entity_id", "type", "pdbx_seq_one_letter_code_can"]):
        entity, kind, sequence = map(value, row)
        if kind == "polypeptide(L)":
            sequence = "".join(sequence.split())
            if not sequence or entity in entities:
                raise ValueError("source protein sequence missing/duplicated")
            entities[entity] = sequence
    chains = {}
    for row in block.find("_struct_asym.", ["id", "entity_id"]):
        label, entity = map(value, row)
        if entity in entities:
            if label in chains:
                raise ValueError("duplicate deposited label chain")
            chains[label] = entity
    columns = ["group_PDB", "type_symbol", "label_asym_id", "label_entity_id", "label_seq_id",
               "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "occupancy", "pdbx_PDB_model_num",
               "Cartn_x", "Cartn_y", "Cartn_z"]
    table = block.find("_atom_site.", columns)
    if not table or not chains:
        raise ValueError("source lacks required protein/coordinate identity categories")
    models = sorted({value(row[9]) for row in table}, key=lambda s: int(s))
    selected_model = "1" if "1" in models else models[0]
    observed, positive = defaultdict(set), defaultdict(set)
    identities = defaultdict(lambda: defaultdict(set))
    partial, zero = defaultdict(set), defaultdict(set)
    atom_records, heterogen_records = 0, 0
    for raw in table:
        row = [value(v) for v in raw]
        if row[9] != selected_model:
            continue
        if row[0] == "HETATM":
            heterogen_records += 1
        if row[1] in {"H", "D"} or row[2] not in chains:
            continue
        label, entity = row[2], row[3]
        if entity != chains[label]:
            raise ValueError("atom-site polymer entity does not match its label chain")
        index = int(row[4])
        if not 1 <= index <= len(entities[entity]):
            raise ValueError("atom-site label position outside deposited protein sequence")
        occupancy = float(row[8])
        if not math.isfinite(occupancy) or not 0 <= occupancy <= 1 or not all(math.isfinite(float(v)) for v in row[10:13]):
            raise ValueError("source coordinates/occupancy invalid")
        observed[label].add(index)
        if occupancy > 0:
            positive[label].add(index)
        if occupancy == 0:
            zero[label].add(index)
        elif occupancy < 1:
            partial[label].add(index)
        identities[label][index].add((row[5], row[6], row[7]))
        atom_records += 1
    by_label = {c["label_chain"]: c for c in metadata["chains"]}
    if set(by_label) != set(chains):
        raise ValueError("metadata and deposited protein label-chain inventory differ")
    result = []
    for label, entity in chains.items():
        expected = by_label[label]
        sequence = entities[entity]
        if sequence != expected["sequence"]:
            raise ValueError("deposited canonical protein sequence changed since metadata collection")
        if any(len(ids) != 1 for ids in identities[label].values()):
            raise ValueError("one label position maps to ambiguous author residue identity")
        gaps = missing_intervals(len(sequence), observed[label])
        meta_gaps = expected["gaps"]
        axis = [{"label_seq_id": i, "author_chain": next(iter(ids))[0],
                 "author_residue_number": next(iter(ids))[1], "insertion_code": next(iter(ids))[2]}
                for i, ids in sorted(identities[label].items())]
        result.append({"label_chain": label, "entity_id": entity, "sequence": sequence,
            "author_chain_ids": sorted({a["author_chain"] for a in axis}), "observed_author_axis": axis,
            "target_residues": len(sequence), "observed_residues_any_heavy_record": len(observed[label]),
            "observed_residues_positive_occupancy": len(positive[label]), "natural_absent_coordinate_spans": gaps,
            "agrees_with_RCSB_unobserved_features": gaps == meta_gaps,
            "residues_with_partial_occupancy": sorted(partial[label]), "residues_with_zero_occupancy_atoms": sorted(zero[label]),
            "terminal_construct_review_required": any(g["kind"] != "internal" for g in gaps),
            "all_observations_retained": True})
    assembly_generators = [{"assembly_id": value(r[0]), "label_chain_list": value(r[1]), "operator_expression": value(r[2])}
                           for r in block.find("_pdbx_struct_assembly_gen.", ["assembly_id", "asym_id_list", "oper_expression"])]
    all_gaps = [g for c in result for g in c["natural_absent_coordinate_spans"]]
    return {"entry_id": entry_id, "selected_coordinate_model_for_inventory": selected_model, "source_model_ids": models,
        "multi_model_source_needs_explicit_policy": len(models) > 1, "chains": result,
        "deposited_target_residues": sum(c["target_residues"] for c in result),
        "deposited_observed_protein_residues": sum(c["observed_residues_any_heavy_record"] for c in result),
        "protein_heavy_coordinate_records": atom_records, "HETATM_records_in_selected_model": heterogen_records,
        "natural_absent_gap_count": len(all_gaps), "longest_absent_gap": max((g["length"] for g in all_gaps), default=0),
        "feature_disagreement_chains": [c["label_chain"] for c in result if not c["agrees_with_RCSB_unobserved_features"]],
        "assembly_generators": assembly_generators, "biological_assembly_materialized": False,
        "scientific_admission_complete": False, "observed_heavy_atom_completeness_checked": False,
        "source_coordinate_changes": False, "checkpoint_sampling_started": False}
