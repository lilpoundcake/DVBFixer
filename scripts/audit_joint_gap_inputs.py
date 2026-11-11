#!/usr/bin/env python3
"""Offline deposited-source / actual pinned motif-consumer audit, no sampling.

Reconstructs target rows independently of the joint builder and executes the
pinned motif loading/placement/crop helpers without loading model weights or
calling model.sample. Keeps input validity distinct from conditioning coverage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    import numpy as np
    import torch
    from Bio import SeqIO
    from Bio.Data.PDBData import protein_letters_3to1
    from protpardelle.common import residue_constants
    from protpardelle.core.models import apply_crop_cond_strategy, parse_fixed_pos_str
    from protpardelle.data.pdb_io import load_feats_from_pdb

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--deposited", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root, out = args.sample, args.output
    parent = json.loads((root / "parent-request.json").read_text())
    view = json.loads((root / "engine-view.json").read_text())
    trace = json.loads((root / "trace.json").read_text())
    group, engine = view["request"], view["engine_view"]["request"]

    def residue(r):
        return r["chain"], r["residue_number"], r.get("insertion_code", "")

    def atom(a):
        return (*residue(a), a["atom_name"])

    def records(path):
        result = {}
        names = {}
        for line in path.read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")) and line[16] in {" ", "A"}:
                identity = (line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
                if identity not in result:
                    result[identity] = np.asarray([float(line[i:i + 8]) for i in (30, 38, 46)])
                    names[identity[:3]] = line[17:20].strip()
        return result, names

    source, source_names = records(root / "source.pdb")
    deposited, _ = records(args.deposited)
    sequences = {t["chain"]: t["sequence"] for t in parent["target_sequences"]}
    seqres = {r.annotations["chain"]: str(r.seq) for r in SeqIO.parse(args.deposited, "pdb-seqres")}
    seqres_checks = {chain: seqres.get(chain) == sequence for chain, sequence in sequences.items()}
    assert all(seqres_checks.values()), "deposited SEQRES differs from parent target"
    target_residues = {}
    mismatches = []
    for placement in parent["sequence_placements"]:
        for index, r in zip(placement["observed_target_indices"], placement["observed_residues"], strict=True):
            identity = residue(r)
            target_residues[placement["chain"], index] = identity
            if protein_letters_3to1.get(source_names[identity]) != sequences[placement["chain"]][index]:
                mismatches.append({"residue": identity, "target_index": index})
    assert not mismatches, "observed amino acids differ from requested target placement"
    gap_for_residue = {}
    for i, gap in enumerate(parent["gaps"]):
        for index, r in zip(range(gap["target_interval"]["start"], gap["target_interval"]["stop"]),
                            gap["generated_residues"], strict=True):
            assert (gap["chain"], index) not in target_residues
            target_residues[gap["chain"], index] = residue(r)
            gap_for_residue[residue(r)] = i
    assert len(target_residues) == sum(len(s) for s in sequences.values())
    original_errors = [float(np.linalg.norm(point - deposited[a])) for a, point in source.items() if a in deposited]
    missing_original = [a for a in source if a not in deposited]
    assert not missing_original and max(original_errors, default=0) == 0, "source coordinates differ from deposit"
    expected_rows = []
    expected_letters = []
    source_indices = []
    spans = group["protein_context"]["spans"]
    mapping = {atom(e): atom(s) for e, s in view["engine_view"]["atom_to_source"]}
    assert len(mapping) == len(set(mapping.values()))
    for span, target, placement in zip(spans, engine["target_sequences"], engine["sequence_placements"], strict=True):
        expected_sequence = sequences[span["chain"]][span["start"]:span["stop"]]
        assert target["sequence"] == expected_sequence
        for index in range(span["start"], span["stop"]):
            expected_rows.append(target_residues[span["chain"], index])
            expected_letters.append(sequences[span["chain"]][index])
            source_indices.append((span["chain"], index))
    generated = {atom(a) for a in parent["generated_atoms"]}
    expected_generated_rows = {i for i, r in enumerate(expected_rows) if r in gap_for_residue}
    expected_fixed_rows = [i for i in range(len(expected_rows)) if i not in expected_generated_rows]
    assert trace["fixed_positions"] == expected_fixed_rows
    assert set(trace["fixed_positions"]).isdisjoint(expected_generated_rows)
    assert len(expected_generated_rows) == sum(len(g["generated_residues"]) for g in group["gaps"])
    assert len(expected_rows) == trace["tokens"]
    # Independent physical feature check, retaining source sequence separations.
    for i in range(1, len(source_indices)):
        if source_indices[i][0] == source_indices[i - 1][0]:
            assert trace["physical_chains"][i] == trace["physical_chains"][i - 1]
            assert trace["physical_positions"][i] - trace["physical_positions"][i - 1] == source_indices[i][1] - source_indices[i - 1][1]
    features, _ = load_feats_from_pdb(root / "motif.pdb", include_pos_feats=True)
    selected = []
    for placement in trace["motif_placement"].split(";/;"):
        segments = [s for s in placement.split("/") if s[0].isalpha()]
        selected.extend(parse_fixed_pos_str(",".join(segments), features["chain_id_mapping"],
            features["residue_index_orig"], features["chain_index"]))
    assert len(selected) == len(set(selected)) == len(expected_fixed_rows)
    assert set(selected) == set(range(len(features["aatype"])))
    actual_positions = features["atom_positions"][selected]
    actual_mask = features["atom_mask"][selected]
    actual_types = features["aatype"][selected]
    expected_aatype = [residue_constants.restype_order[expected_letters[i]] for i in expected_fixed_rows]
    assert actual_types.tolist() == expected_aatype, "pinned consumer reordered motif amino acids"
    expected_atoms = {atom(c["atom"]): np.asarray(c["coordinate_angstrom"])
                      for c in group["protein_context"]["fixed_conditioning"]}
    consumer_atoms = set()
    max_error = 0.0
    for row, index in enumerate(expected_fixed_rows):
        r = expected_rows[index]
        for name, col in residue_constants.atom_order.items():
            identity = (*r, name)
            present = bool(actual_mask[row, col])
            assert present == (identity in expected_atoms), "pinned consumer changed fixed atom mask"
            if present:
                consumer_atoms.add(identity)
                max_error = max(max_error, float(np.linalg.norm(actual_positions[row, col].numpy() - expected_atoms[identity])))
    assert consumer_atoms == set(expected_atoms) and max_error < 1e-5
    # Replay before the sampler's common random rigid rotation. A shared rotation
    # preserves these identity/mask/pairwise geometry checks; it is NOT recovered
    # from the old trajectory or claimed as an archived denoiser input tensor.
    centered = (actual_positions - actual_positions[:, 1:2, :].mean(dim=0, keepdim=True)) * actual_mask.unsqueeze(-1)
    conditional = torch.zeros((1, len(expected_rows), 37, 3), dtype=centered.dtype)
    conditional[0, expected_fixed_rows] = centered
    aatypes = [[residue_constants.restype_1to3[residue_constants.restypes[i]] for i in actual_types.tolist()]]
    cropped = apply_crop_cond_strategy(conditional, [expected_fixed_rows], aatypes, "backbone-sidechain")
    assert torch.equal(cropped, conditional)
    assert torch.count_nonzero(cropped[0, sorted(expected_generated_rows)]) == 0
    np.savez_compressed(out / "replayed-consumer-conditioning.npz", atom_positions=cropped.numpy(),
                        fixed_rows=np.asarray(expected_fixed_rows), generated_rows=np.asarray(sorted(expected_generated_rows)))
    # Per-gap source coverage is reported separately from semantic correctness.
    selected_residues = {a[:3] for a in consumer_atoms}
    omissions = [r for r in source_names if r not in selected_residues]
    coverage = {}
    for number in (150, 241, 242, 284, 474):
        key = ("A", str(number), "")
        coverage[str(number)] = {"joint_conditioned": key in selected_residues,
            "independent_contexts": [i for i, plan in enumerate(parent["gap_contexts"])
                if any(residue(c["atom"]) == key for c in plan["fixed_conditioning"])]}
    token_table = [{"token": i, "source_residue": r, "amino_acid": expected_letters[i],
                    "source_target_index": source_indices[i][1], "generated": i in expected_generated_rows,
                    "gap_index": gap_for_residue.get(r)} for i, r in enumerate(expected_rows)]
    (out / "tokens.json").write_text(json.dumps(token_table, indent=2) + "\n")
    (out / "omitted-observed-residues.json").write_text(json.dumps(omissions, indent=2) + "\n")
    report = {"scope": __doc__, "new_sampling": False, "model_weights_loaded": False,
        "deposited_seqres_matches": seqres_checks, "observed_sequence_mismatches": mismatches,
        "source_coordinate_max_error_from_deposit_A": max(original_errors, default=0),
        "tokens": len(expected_rows), "generated_residues": len(expected_generated_rows),
        "fixed_residue_rows": len(expected_fixed_rows), "fixed_conditioning_atoms": len(consumer_atoms),
        "pinned_consumer_mask_order_sequence_and_coordinates_match": True,
        "consumer_replay_frame": "before common random rigid rotation; not an archived denoiser tensor",
        "consumer_coordinate_max_error_A": max_error, "fixed_generated_overlap": False,
        "physical_positions_and_chain_separations_match": True,
        "conditioning_by_chain": dict(Counter(r[0] for r in selected_residues)),
        "omitted_observed_by_chain": dict(Counter(r[0] for r in omissions)),
        "problem_partner_coverage": coverage,
        "source_sha256": hashlib.sha256((root / "source.pdb").read_bytes()).hexdigest(),
        "deposited_sha256": hashlib.sha256(args.deposited.read_bytes()).hexdigest(),
        "input_integrity_proves_structural_accuracy": False}
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
