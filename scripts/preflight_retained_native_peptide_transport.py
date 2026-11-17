#!/usr/bin/env python3
"""Data-only native covalent graph audit; valid masks do not prove peptide transport."""
from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

from scripts import preflight_retained_loop_protenix_native as native
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import AtomIdentity


def junction_transport(request, axis, native_bonds, native_chains):
    if len(set(axis)) != len(axis) or len(axis) != len(native_chains):
        raise ValueError("native graph identity axis invalid")
    index = {atom: i for i, atom in enumerate(axis)}
    edges = {frozenset((int(row[0]), int(row[1]))) for row in native_bonds}
    rows = []
    for gap_index, gap in enumerate(request.gaps):
        for side, left, right in (("left", gap.left_anchor, gap.generated_residues[0]),
                                  ("right", gap.generated_residues[-1], gap.right_anchor)):
            if left is None or right is None:
                continue
            first = AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")
            second = AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")
            if first not in index or second not in index:
                raise ValueError("native context omitted intended junction endpoint")
            a, b = index[first], index[second]
            rows.append({"gap_index": gap_index, "side": side,
                "source_atoms": [asdict(first), asdict(second)], "native_indices": [a, b],
                "native_chains": [str(native_chains[a]), str(native_chains[b])],
                "same_native_polymer_chain": bool(native_chains[a] == native_chains[b]),
                "native_covalent_bond_present": frozenset((a, b)) in edges})
    return {"junctions": rows, "intended_internal_junctions": len(rows),
        "native_bonds_present": sum(row["native_covalent_bond_present"] for row in rows),
        "cross_native_chain_junctions": sum(not row["same_native_polymer_chain"] for row in rows),
        "structural_template_mask_validity_is_not_polymer_continuity": True,
        "causality_of_raw_geometry_failure_established": False, "public_admission_or_repair_changed": False}


def worker(case, output):
    sys.path.insert(0, str(native.CLONE))
    sys.path.insert(0, str(native.PROJECT / "deploy/protenix-v1"))
    import socket

    import torch

    def forbidden(*args, **kwargs):
        raise RuntimeError("native graph preflight forbids model/checkpoint/network")

    torch.load = forbidden
    socket.create_connection = socket.socket.connect = socket.socket.connect_ex = forbidden
    handle = torch.nn.modules.module.register_module_forward_pre_hook(forbidden)
    import checkpoint_gap_smoke as helper

    original_axis, original_save = helper._model_atom_axis, native.save
    evidence = {}

    def capture(array, request):
        result = original_axis(array, request)
        if array.bonds is None:
            raise ValueError("native atom array lacks a covalent graph")
        bonds = array.bonds.as_array()
        evidence.update(junction_transport(request, result[0], bonds, array.chain_id))
        save(output / "native-covalent-transport.json", {**evidence,
            "native_atom_axis": [asdict(atom) for atom in result[0]],
            "native_bonds": bonds.tolist(), "native_chain_ids": array.chain_id.tolist()})
        return result

    def enriched_save(path, value):
        if path.name == "result.json":
            value = {**value, "native_peptide_transport": evidence}
        original_save(path, value)

    helper._model_atom_axis, native.save = capture, enriched_save
    try:
        native.worker(case, output)
    finally:
        helper._model_atom_axis, native.save = original_axis, original_save
        handle.remove()


def declare(output):
    encoded = native.PROJECT / "tmp/mac-real-gap-continuation-2026-11-12/retained-loop-protenix-encoder-attempt01"
    files = [Path(__file__), Path(native.__file__), native.PROJECT / "deploy/protenix-v1/checkpoint_gap_smoke.py"]
    files += list((native.CLONE / "configs").rglob("*.py")) + list((native.CLONE / "protenix/data").rglob("*.py"))
    files += [p for p in encoded.rglob("*") if p.is_file()]
    output.mkdir(parents=True, exist_ok=False)
    jobs = [{"case_id": case+"-native-covalent", "case": str(encoded / case), "seed": 0,
        "command": [str(native.NATIVE_PYTHON), "-m", "scripts.preflight_retained_native_peptide_transport", "worker",
            "--case", str(encoded / case), "--output", str(output / "execution" / f"{case}-native-covalent-seed-0")]}
        for case in ("1NON", "2BE1")]
    save(output / "plan.json", {"python": str(native.NATIVE_PYTHON), "python_sha256": digest(native.NATIVE_PYTHON),
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
        "scope": "data-only actual atom-array covalent graph for previously encoded context; no corrected input or new proposals",
        "checkpoint_or_model_calls": 0, "limits_changed": False,
        "resource_guard": "existing serial CPU1800s/job,RSS8GiB,pressure/headroom guard"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "worker"))
    parser.add_argument("--case", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    else:
        worker(args.case.resolve(), args.output.resolve())
