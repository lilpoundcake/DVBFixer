#!/usr/bin/env python3
"""Export one visualization-only multi-model PDB per gap from portfolio donors."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.selection
    request_bytes = (root / "request.json").read_bytes()
    request = json.loads(request_bytes)
    source = (root / "source.pdb").read_bytes()
    parents = json.loads((root / "parents.json").read_text())
    selected = json.loads((root / "result.json").read_text())["selected"]["choices"]
    selected_seeds = {c["gap_index"]: parents[c["parent_index"]]["seed"] for c in selected}

    def identity(atom):
        return atom["chain"], atom["residue_number"], atom.get("insertion_code", ""), atom["atom_name"]

    def line_identity(line):
        return line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()

    source_records = {line_identity(line): line for line in source.decode().splitlines(keepends=True)
                      if line.startswith(("ATOM  ", "HETATM"))}
    positions = {}
    for placement in request["sequence_placements"]:
        for index, residue in zip(placement["observed_target_indices"], placement["observed_residues"], strict=True):
            positions[residue["chain"], residue["residue_number"], residue.get("insertion_code", "")] = index
    for gap in request["gaps"]:
        for index, residue in zip(range(gap["target_interval"]["start"], gap["target_interval"]["stop"]),
                                  gap["generated_residues"], strict=True):
            positions[residue["chain"], residue["residue_number"], residue.get("insertion_code", "")] = index
    donor_records = []
    for parent in parents:
        data = (root / f'parent-{parent["seed"]}' / "stage1.pdb").read_bytes()
        assert hashlib.sha256(data).hexdigest() == parent["coordinate_sha256"]
        assert hashlib.sha256(request_bytes).hexdigest() == parent["request_sha256"]
        assert hashlib.sha256(source).hexdigest() == parent["source_sha256"]
        assert parent["stage"] == "generated-only"
        lines = [line for line in data.decode().splitlines(keepends=True) if line.startswith(("ATOM  ", "HETATM"))]
        ids = [line_identity(line) for line in lines]
        assert len(ids) == len(set(ids))
        for line in lines:
            atom = line_identity(line)
            if atom in source_records:
                assert line[30:54] == source_records[atom][30:54]
        donor_records.append(lines)
    outputs = []
    for i, (gap, context) in enumerate(zip(request["gaps"], request["gap_contexts"], strict=True)):
        generated_residues = {(r["chain"], r["residue_number"], r.get("insertion_code", "")) for r in gap["generated_residues"]}
        generated = {identity(a) for a in request["generated_atoms"] if identity(a)[:3] in generated_residues}
        fixed = {identity(c["atom"]) for c in context["fixed_conditioning"]}
        wanted = fixed | generated
        first, last = gap["generated_residues"][0], gap["generated_residues"][-1]
        filename = f'gap-{i:02d}_{gap["chain"]}{first["residue_number"]}{first.get("insertion_code", "")}-{last["residue_number"]}{last.get("insertion_code", "")}_all-seeds.pdb'
        contents = ["REMARK 999 VISUALIZATION ONLY: GENERATED-ONLY PORTFOLIO DONORS\n",
            "REMARK 999 LOCAL OBSERVED MOTIF + ONE COMPLETE GENERATED GAP\n",
            "REMARK 999 ORIGINAL FRAME; NO FITTING, MINIMIZATION OR RESAMPLING\n",
            f'REMARK 999 GAP INDEX {i}; PATCHWORK SELECTED SEED {selected_seeds[i]}\n']
        for model, (parent, lines) in enumerate(zip(parents, donor_records, strict=True), 1):
            contents.append(f'REMARK 999 MODEL {model} = SEED {parent["seed"]}\n')
        for model, lines in enumerate(donor_records, 1):
            contents.append(f"MODEL     {model:4d}\n")
            found = set()
            previous_residue = None
            for line in lines:
                atom = line_identity(line)
                if atom not in wanted:
                    continue
                residue = atom[:3]
                if previous_residue is not None and residue != previous_residue:
                    if residue[0] != previous_residue[0] or positions[residue] != positions[previous_residue] + 1:
                        contents.append("TER\n")
                contents.append(line if line.endswith("\n") else line + "\n")
                found.add(atom)
                previous_residue = residue
            assert found == wanted, (i, model, wanted - found)
            contents.extend(("TER\n", "ENDMDL\n"))
        contents.append("END\n")
        data = "".join(contents).encode()
        with (args.output / filename).open("xb") as stream:
            stream.write(data)
        outputs.append({"file": filename, "gap_index": i, "chain": gap["chain"],
            "first_generated_residue": first, "last_generated_residue": last,
            "models": len(parents), "generated_atoms_per_model": len(generated),
            "observed_motif_atoms_per_model": len(fixed), "selected_patchwork_seed": selected_seeds[i],
            "sha256": hashlib.sha256(data).hexdigest()})
    manifest = {"scope": __doc__, "stage": "generated-only, before multi-parent selection/repair",
        "model_to_seed": {str(i): p["seed"] for i, p in enumerate(parents, 1)},
        "source_sha256": hashlib.sha256(source).hexdigest(), "donors": parents,
        "original_frame": True, "coordinates_changed": False, "new_sampling": False,
        "publication_qualified": False, "files": outputs}
    (args.output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    (args.output / "README.txt").write_text(
        "1HE8: one multi-model PDB per original gap/local motif.\n"
        "MODEL 1=seed 7; 2=101; 3=102; 4=103; 5=104.\n"
        "Each model: the same observed local conditioning motif + this gap's generated atoms.\n"
        "Stage: generated-only refinement, the actual donors used by portfolio selection.\n"
        "All coordinates retained exactly; no per-fragment fitting or further minimization.\n"
        "Original chain/residue IDs retained; TER separates omitted sequence spans.\n"
        "Load all models as an ensemble to compare, or switch models individually.\n"
        "Geometry failures remain included. This is visualization, not five accepted structures.\n")
    print(json.dumps({"output": str(args.output), "files": len(outputs), "models_per_file": len(parents)}))


if __name__ == "__main__":
    main()
