#!/usr/bin/env python3
"""Explicit contiguous local projection from an existing immutable whole-gap request.

No reference coordinates, generated donors or guessed sequence placement enter
conditioning. Sequence identities use the request's established placement. The
span includes every observed residue between selected gaps and fixed flanks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS


def prepare(request_path, source_path, indices, output, *, flanks=8, include_declared_context=False):
    request = DiffusionRequest.from_json(request_path.read_text())
    source_bytes = source_path.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("source digest differs from established request")
    if not indices or len(set(indices)) != len(indices) or any(not 0 <= i < len(request.gaps) for i in indices):
        raise ValueError("select valid unique whole gaps")
    gaps = tuple(request.gaps[i] for i in indices)
    chain = gaps[0].chain
    if any(g.chain != chain for g in gaps) or flanks < 1:
        raise ValueError("this contiguous projection requires one chain and positive flanks")
    target = next(t for t in request.target_sequences if t.chain == chain)
    placement = next(p for p in request.sequence_placements if p.chain == chain)
    mapping = dict(zip(placement.observed_target_indices, placement.observed_residues))
    for gap in request.gaps:
        if gap.chain != chain:
            continue
        for index, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues):
            if index in mapping:
                raise ValueError("overlapping established target placement")
            mapping[index] = residue
    start = max(0, min(g.target_interval.start for g in gaps)-flanks)
    stop = min(len(target.sequence), max(g.target_interval.stop for g in gaps)+flanks)
    selected_indices = set(range(start, stop))
    if include_declared_context:
        if request.gap_contexts is None:
            raise ValueError("request has no established local context plans")
        for gap_index in indices:
            if request.gap_contexts[gap_index].source_sha256 != request.normalized_pdb.sha256:
                raise ValueError("established context plan is bound to another source")
            for span in request.gap_contexts[gap_index].spans:
                if span.chain != chain:
                    raise ValueError("this research adapter does not yet transport multiple native chains")
                selected_indices.update(range(span.start, span.stop))
    if len(selected_indices) > 256:
        raise ValueError("local projection exceeds prospectively bounded 256 residues")
    for i, gap in enumerate(request.gaps):
        if gap.chain == chain and selected_indices & set(range(gap.target_interval.start, gap.target_interval.stop)) and i not in indices:
            raise ValueError("projection contains an unselected gap; declare it or use another span")
    ordered_indices = tuple(sorted(selected_indices))
    ranges = []
    for index in ordered_indices:
        if not ranges or ranges[-1][1] != index:
            ranges.append([index, index+1])
        else:
            ranges[-1][1] += 1
    projected = tuple(mapping[i] for i in ordered_indices)
    if any(r.insertion_code for r in projected):
        raise ValueError("insertion-code RFD3 contig selection is not yet audited")
    for left, right in ranges:
        if [int(mapping[i].residue_number) for i in range(left, right)] != list(range(int(mapping[left].residue_number), int(mapping[right-1].residue_number)+1)):
            raise ValueError("RFD3 contiguous numeric selection would not preserve this placement")
    residue_keys = {(r.chain, r.residue_number, r.insertion_code) for r in projected}
    if any((a.chain, a.residue_number, a.insertion_code) in residue_keys
           for link in request.retained_explicit_links for a in (link.atom1, link.atom2)):
        raise ValueError("local RFD3 adapter does not yet transport retained explicit crosslinks")
    observed = {}
    for line in source_bytes.decode().splitlines(True):
        if not line.startswith(("ATOM  ", "HETATM")) or (line[21], line[22:26].strip(), line[26].strip()) not in residue_keys:
            continue
        identity = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if identity in observed or line[16].strip() or float(line[54:60]) != 1:
            raise ValueError("ambiguous/partial-occupancy source atom")
        observed[identity] = line
    generated = {a for a in request.generated_atoms if (a.chain, a.residue_number, a.insertion_code) in residue_keys}
    rows, axis, mask = [], [], {}
    for index, residue in zip(ordered_indices, projected):
        name = _CANONICAL_ONE_TO_THREE[target.sequence[index]]
        expected = {AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, a) for a in CANONICAL_HEAVY_ATOMS[name]}
        actual = {a for a in observed if (a.chain, a.residue_number, a.insertion_code) ==
                  (residue.chain, residue.residue_number, residue.insertion_code)}
        if actual | (expected & generated) != expected or actual & generated:
            raise ValueError(f"source lacks canonical complete observations at {residue}")
        mask[f"{chain}{residue.residue_number}"] = "" if expected <= generated else "ALL"
        for atom_name in CANONICAL_HEAVY_ATOMS[name]:
            atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
            if atom in observed:
                line = observed[atom]
                if line[17:20].strip() != name:
                    raise ValueError("observed sequence differs from established placement")
            else:
                field = f" {atom_name:<3}" if len(atom_name) < 4 else atom_name
                line = (f"ATOM  {len(rows)+1:5d} {field} {name:3s} {chain}{int(residue.residue_number):4d}    "
                        f"{0:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}          {atom_name[0]:>2s}\n")
            rows.append(line[:6]+f"{len(rows)+1:5d}"+line[11:])
            axis.append({"atom": {"chain": atom.chain, "residue_number": atom.residue_number,
                         "insertion_code": atom.insertion_code, "atom_name": atom.atom_name},
                         "residue_name": name, "generated": atom in generated})
    output.mkdir(parents=True, exist_ok=False)
    header = "REMARK 999 EXPLICIT CONTIGUOUS LOCAL PROJECTION; NOT FULL SOURCE\n"
    (output / "conditioning.pdb").write_text(header+"".join(rows)+"TER\nEND\n")
    (output / "source-projection.pdb").write_text(header+"".join(observed.values())+"TER\nEND\n")
    spec = {"local_projection": {"input": str((output / "conditioning.pdb").absolute()),
        "contig": ",/0,".join(f"{chain}{mapping[left].residue_number}-{mapping[right-1].residue_number}" for left, right in ranges),
        "select_fixed_atoms": mask, "select_unfixed_sequence": False, "plddt_enhanced": True}}
    (output / "rfd3-input.json").write_text(json.dumps(spec, indent=2)+"\n")
    manifest = {"source_sha256": request.normalized_pdb.sha256, "request_sha256": hashlib.sha256(request_path.read_bytes()).hexdigest(),
        "selected_gap_indices": indices, "flanks": flanks, "target_interval": {"start": start, "stop": stop},
        "sequence": "".join(target.sequence[i] for i in ordered_indices), "generated_atoms": len(generated), "observed_atoms": len(observed),
        "source_target_spans": ranges, "declared_context_included": include_declared_context,
        "model_entity_cuts": "explicit /0 between noncontiguous source segments; source identity/frame retained, not full-chain topology" if len(ranges)>1 else "none within contiguous source span",
        "axis": axis, "heldout_reference_supplied": False, "other_generated_donors_supplied": False,
        "omitted_source_not_guaranteed_irrelevant": True, "raw_generated_stage_only": True}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "source", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--gaps", type=int, nargs="+", required=True)
    parser.add_argument("--declared-context", action="store_true")
    args = parser.parse_args()
    result = prepare(args.request, args.source, tuple(args.gaps), args.output, include_declared_context=args.declared_context)
    print(json.dumps({k:v for k,v in result.items() if k != "axis"}))
