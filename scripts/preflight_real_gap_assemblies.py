#!/usr/bin/env python3
"""Guarded assembly/altloc/observed-chemistry inventory; no source repair or ML."""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.altloc_selection import select_altlocs
from dvbfixer.biological_assembly import parse_biological_assemblies, render_biological_assembly
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.scope import _parse_coordinate_records, _parse_explicit_links
from dvbfixer.structure_input import PDB_CHAIN_IDS, normalize_structure

PROJECT = Path(__file__).resolve().parents[1]


def inspect(row, directory):
    directory.mkdir(parents=True, exist_ok=False)
    source = Path(row["source_path"])
    result = {"entry_id": row["entry_id"], "source_sha256": row["source_sha256"],
        "scientific_admission_complete": False, "model_sampling": False, "observed_atom_repair": False}
    stage = "normalization"
    try:
        if digest(source) != row["source_sha256"]:
            raise RuntimeError("original coordinate bytes changed")
        normalized = normalize_structure(source, directory / "normalized.pdb")
        text = normalized.pdb_path.read_text()
        stage = "assembly1"
        lines = text.splitlines(keepends=True)
        assemblies = parse_biological_assemblies(lines)
        if "1" not in assemblies:
            raise ValueError("no declared assembly1; no asymmetric-unit substitution")
        assembly = assemblies["1"]
        source_counts = Counter(line[21] for line in lines if line.startswith(("ATOM  ", "HETATM")))
        copies = sum(len(group.chains)*len(group.operators()) for group in assembly.groups)
        atom_count = sum(sum(source_counts[c] for c in group.chains)*len(group.operators()) for group in assembly.groups)
        if atom_count+copies > 99999:
            raise ValueError("complete biological assembly exceeds PDB serial capacity; no cropping")
        rendered = "".join(render_biological_assembly(lines, assembly, list(PDB_CHAIN_IDS), renumber=False))
        (directory / "assembly1.pdb").write_text(rendered)
        stage = "altloc_selection"
        selected, altloc = select_altlocs(rendered, source_sha256=digest(directory / "assembly1.pdb"), policy="occupancy")
        save(directory / "altloc-report.json", altloc)
        if selected is None:
            raise ValueError("residue-coherent occupancy selection refused")
        (directory / "selected.pdb").write_text(selected)
        stage = "source_route_inventory"
        selected_lines = selected.splitlines()
        atoms, _, _ = _parse_coordinate_records(selected_lines)
        links = _parse_explicit_links(selected_lines, atoms)
        by_residue = defaultdict(list)
        for line in selected_lines:
            if line.startswith(("ATOM  ", "HETATM")):
                key = line[21], line[22:26].strip(), line[26].strip()
                by_residue[key].append(line)
        protein_keys = {key for key, group in by_residue.items()
                        if any(line.startswith("ATOM  ") and line[17:20].strip() in CANONICAL_HEAVY_ATOMS for line in group)}
        unknown_atom_residues = sorted({line[17:20].strip() for line in selected_lines
            if line.startswith("ATOM  ") and line[17:20].strip() not in CANONICAL_HEAVY_ATOMS})
        partial = [line for key, group in by_residue.items() if key in protein_keys for line in group
                   if line[54:60].strip() and float(line[54:60]) < 1]
        incomplete = []
        for key in sorted(protein_keys):
            group = by_residue[key]
            names = {line[17:20].strip() for line in group}
            if len(names) != 1:
                raise ValueError("one residue key contains mixed residue identities")
            name = next(iter(names))
            present = {line[12:16].strip() for line in group}
            missing = sorted(set(CANONICAL_HEAVY_ATOMS[name])-present)
            if missing:
                incomplete.append({"chain": key[0], "residue_number": key[1], "insertion_code": key[2],
                                   "residue_name": name, "missing_heavy_atoms": missing})
        def protein(atom):
            return (atom.chain, atom.residue_number, atom.insertion_code) in protein_keys
        linked_heterogens = [(first, second) for first, second in links if protein(first) != protein(second)]
        atom_records = {atom.identity: atom for atom in atoms}
        normalized_polymer_keys = {(normalized.chain_map.get(a["author_chain"], a["author_chain"]),
                                   a["author_residue_number"], a["insertion_code"])
                                  for chain in row["inventory"]["chains"] for a in chain["observed_author_axis"]}
        noncanonical_polymer = sorted({line[17:20].strip() for line in lines if line.startswith(("ATOM  ", "HETATM"))
            and (line[21], line[22:26].strip(), line[26].strip()) in normalized_polymer_keys
            and line[17:20].strip() not in CANONICAL_HEAVY_ATOMS})
        obstacles = []
        if partial:
            obstacles.append("selected protein observations retain partial/zero occupancy")
        if unknown_atom_residues or noncanonical_polymer:
            obstacles.append("noncanonical polymer/ATOM chemistry requires explicit routing")
        if linked_heterogens:
            obstacles.append("protein/nonprotein explicit links need covalent-versus-coordination review")
        result.update(status="assembly_inventory_complete", assembly_id="1", chain_map=normalized.chain_map,
            assembly_atom_records=atom_count, canonical_observed_protein_residues=len(protein_keys),
            partial_or_zero_occupancy_protein_atoms=len(partial), missing_observed_heavy=incomplete,
            needs_observed_heavy_completion=bool(incomplete), unknown_ATOM_residue_names=unknown_atom_residues,
            noncanonical_polymer_residue_names=noncanonical_polymer,
            protein_nonprotein_explicit_links=len(linked_heterogens), strict_protein_route_obstacles=obstacles,
            protein_nonprotein_link_details=[{"first": asdict(first), "second": asdict(second),
                "first_residue_name": atom_records[first].residue_name,
                "second_residue_name": atom_records[second].residue_name,
                "chemical_bond_classification": "not established by LINK/CONECT existence alone"}
                for first, second in sorted(linked_heterogens)],
            heterogen_atom_records=sum(a.record_name == "HETATM" for a in atoms),
            source_geometry_checked=False, all_original_context_still_retained=True,
            selected_sha256=digest(directory / "selected.pdb"))
    except Exception as exc:
        result.update(status="assembly_inventory_refused", stage=stage, error_type=type(exc).__name__, error=str(exc))
        (directory / "failure.txt").write_text(traceback.format_exc())
    save(directory / "inventory.json", result)
    return result


def worker(plan_root, batch, output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("assembly worker requires fresh executor-owned reservation")
    rows = json.loads((plan_root / "batches.json").read_text())[batch]
    results = [inspect(row, output / "cases" / row["entry_id"]) for row in rows]
    save(output / "result.json", {"rows": results, "case_count": len(results),
        "no_checkpoint_sampling_or_repair": True})


def declare(sources, output):
    completion = json.loads((sources / "completed.json").read_text())
    raw = sources / "source-results.json"
    if digest(raw) != completion["source_results_sha256"]:
        raise ValueError("source inventory changed")
    rows = json.loads(raw.read_text())["rows"]
    if len(rows) != completion["downloaded"]:
        raise ValueError("original source inventory incomplete")
    output.mkdir(parents=True, exist_ok=False)
    batches = [rows[start:start+64] for start in range(0, len(rows), 64)]
    save(output / "batches.json", batches)
    jobs = [{"case_id": f"assembly-batch-{index:02d}", "case": str(raw), "seed": 0,
        "command": [sys.executable, "-m", "scripts.preflight_real_gap_assemblies", "worker", "--batch", str(index),
            "--plan", str(output), "--output", str(output / "execution" / f"assembly-batch-{index:02d}-seed-0")]}
        for index in range(len(batches))]
    files = [Path(__file__), raw, output / "batches.json", PROJECT / "src/dvbfixer/structure_input.py",
        PROJECT / "src/dvbfixer/altloc_selection.py", PROJECT / "src/dvbfixer/biological_assembly.py",
        PROJECT / "src/dvbfixer/model/diffusion/scope.py", PROJECT / "src/dvbfixer/model/diffusion/protein_atoms.py"]
    files += [Path(r["source_path"]) for r in rows]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)), "jobs": jobs,
        "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "count": len(jobs),
        "original_sources": len(rows), "seed_zero_is_executor_reservation_not_sampler_seed": True,
        "purpose": "assembly/altloc/observed inventory only; no inference, chemistry rewrite or observed repair",
        "resource_guard": "serial CPU executor, 1800s/batch, RSS8GiB, sustained pressure/headroom, 10GiB evidence cap"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sources", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--batch", type=int)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.sources.resolve(), args.output.resolve())
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        worker(args.plan.resolve(), args.batch, args.output.resolve())
