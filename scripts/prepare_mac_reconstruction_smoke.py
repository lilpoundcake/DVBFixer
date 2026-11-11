#!/usr/bin/env python3
"""Prepare a labeled 1UBQ A1-72 / gap10-13 Mac smoke without reference leakage.

This is a new projected-input research case, never a modified frozen replay.
Missing-coordinate atoms in the RFD3 input carry sequence-derived identities and
zero dummy coordinates, not held-out native coordinates. Only fixed observations
are supplied in the other backend templates. The reference is stored separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def prepare(reference: Path, output: Path):
    data = reference.read_bytes()
    records = {}
    names = {}
    for line in data.decode().splitlines(keepends=True):
        if not line.startswith("ATOM  ") or line[21] != "A" or not 1 <= int(line[22:26]) <= 72:
            continue
        if line[16] != " " or line[26] != " " or float(line[54:60]) != 1.0:
            raise ValueError("selected reference has alternate identity/partial occupancy")
        residue = int(line[22:26])
        name = line[12:16].strip()
        key = residue, name
        if key in records:
            raise ValueError("duplicate atom")
        records[key] = line
        names[residue] = line[17:20].strip()
    if set(names) != set(range(1, 73)):
        raise ValueError("explicit A1-72 projection is incomplete")
    from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE

    one_letter = {v: k for k, v in _CANONICAL_ONE_TO_THREE.items()}
    sequence = "".join(one_letter[names[i]] for i in range(1, 73))
    for i in range(1, 73):
        expected = set(CANONICAL_HEAVY_ATOMS[names[i]])
        if {name for residue, name in records if residue == i} != expected:
            raise ValueError(f"noncanonical/incomplete source atom inventory at A{i}")
    observed = {a: line for a, line in records.items() if not 10 <= a[0] <= 13}
    seqres = [f"SEQRES {row:3d} A   72  {' '.join(names[i] for i in range(start, min(start + 13, 73)))}\n"
              for row, start in enumerate(range(1, 73, 13), 1)]
    header = ["HEADER    EXPLICIT 1UBQ A1-72 RESEARCH PROJECTION\n",
              "REMARK 999 A10-13 COORDINATES WITHHELD; NOT ORIGINAL FULL ENTRY\n", *seqres]
    source_text = "".join([*header, *observed.values(), "TER\n", "END\n"])
    dummy = []
    for i in range(1, 73):
        if not 10 <= i <= 13:
            dummy.extend(line for (residue, _), line in records.items() if residue == i)
        else:
            for k, name in enumerate(CANONICAL_HEAVY_ATOMS[names[i]]):
                atom_field = f" {name:<3}" if len(name) < 4 else name
                dummy.append(f"ATOM  {9000 + i * 20 + k:5d} {atom_field} {names[i]:3s} A{i:4d}    "
                             f"{0:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}          {name[0]:>2s}\n")
    # AtomWorks/Biotite requires monotonically increasing serials. This adapter
    # input has no CONECT records; reserialize only this explicit dummy projection,
    # never source/reference records or scientific identity/coordinate fields.
    dummy = [line[:6] + f"{i:5d}" + line[11:] for i, line in enumerate(dummy, 1)]
    output.mkdir(parents=True, exist_ok=False)

    def write(name, content):
        value = content if isinstance(content, bytes) else content.encode()
        with (output / name).open("xb") as stream:
            stream.write(value)
        return hashlib.sha256(value).hexdigest()

    files = {}
    files["reference-projection-NOT-CONDITIONING.pdb"] = write("reference-projection-NOT-CONDITIONING.pdb",
        "".join([*header, *records.values(), "TER\n", "END\n"]))
    files["source.pdb"] = write("source.pdb", source_text)
    files["rfd3-sequence-dummy.pdb"] = write("rfd3-sequence-dummy.pdb", "".join([*header, *dummy, "TER\n", "END\n"]))
    files["target.fasta"] = write("target.fasta", ">A explicit_1UBQ_A1-72\n" + sequence + "\n")
    spec = {"ubq72_gap10_13": {"input": str((output / "rfd3-sequence-dummy.pdb").absolute()),
        "contig": "A1-72", "select_fixed_atoms": {"A1-9": "ALL", "A10-13": "", "A14-72": "ALL"},
        "select_unfixed_sequence": False,
        "plddt_enhanced": True}}
    files["rfd3-input.json"] = write("rfd3-input.json", json.dumps(spec, indent=2) + "\n")
    boltz = {"version": 1, "sequences": [{"protein": {"id": "A", "sequence": sequence, "msa": "empty"}}],
        "templates": [{"pdb": str((output / "source.pdb").absolute()), "chain_id": "A", "template_id": "A1"}]}
    files["boltz-input.yaml"] = write("boltz-input.yaml", json.dumps(boltz, indent=2) + "\n")
    manifest = {"scope": __doc__, "reference_entry_sha256": hashlib.sha256(data).hexdigest(),
        "projection": "A1-72", "withheld": "A10-13", "sequence": sequence, "seed": 7,
        "observed_atoms": len(observed), "generated_atoms": len(records) - len(observed),
        "rfd3_hidden_coordinates": "all zero dummy; sequence-derived canonical names, no native coordinates",
        "boltz_mode": "template-conditioned proposal, NOT established strict coordinate-preserving inpainting",
        "files": files, "new_sampling": False}
    write("manifest.json", json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.reference, args.output)))
