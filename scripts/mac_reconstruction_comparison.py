"""Private reference-blind contracts for the next Mac comparison slice.

Generation and identity transport read source/request only. Held-out coordinates
belong exclusively to post-generation scoring. No public capability is enabled.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE

DEVELOPMENT_SEEDS = (7, 42, 101, 102, 103, 104)
VALIDATION_SEEDS = (211, 223, 227, 229, 233)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    serialized = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    with path.open("x") as stream:
        stream.write(serialized)


def target_axis(request: DiffusionRequest):
    """Exact residue ordinals from the frozen mask, never from reference fitting."""
    result = {}
    for target in request.target_sequences:
        placement = next(p for p in request.sequence_placements if p.chain == target.chain)
        mapping = dict(zip(placement.observed_target_indices, placement.observed_residues))
        for gap in request.gaps:
            if gap.chain != target.chain:
                continue
            for index, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues):
                if index in mapping:
                    raise ValueError("overlapping target mask")
                mapping[index] = residue
        if set(mapping) != set(range(len(target.sequence))):
            raise ValueError("request does not map every target ordinal")
        result[target.chain] = tuple(mapping[i] for i in range(len(target.sequence)))
    return result


def modeller_alignment(request: DiffusionRequest) -> str:
    """Use the already established sequence placement, not another align2d fit."""
    target_axis(request)
    template, target = [], []
    for sequence in request.target_sequences:
        placement = next(p for p in request.sequence_placements if p.chain == sequence.chain)
        observed = set(placement.observed_target_indices)
        template.append("".join(symbol if i in observed else "-" for i, symbol in enumerate(sequence.sequence)))
        target.append(sequence.sequence)
    return (">P1;source\nstructureX:source:FIRST:@:LAST:@::::\n" + "/".join(template) + "*\n"
            + ">P1;target\nsequence:target::::::::\n" + "/".join(target) + "*\n")


def normalize_modeller(request: DiffusionRequest, raw: str):
    """Verified ordinal transport; report and exclude only terminal extra OXT.

    The unmodified native PDB must be retained by the caller. No coordinates are
    fitted, reflected, restored or read from the held-out reference here.
    """
    axis = target_axis(request)
    chains = {}
    for line in raw.splitlines(keepends=True):
        if line.startswith("HETATM"):
            raise ValueError("MODELLER comparator is canonical protein-only")
        if not line.startswith("ATOM  "):
            continue
        if line[16].strip() or float(line[54:60]) != 1:
            raise ValueError("ambiguous native atom")
        chain = chains.setdefault(line[21], {})
        chain.setdefault((line[22:26], line[26]), []).append(line)
    if len(chains) != len(request.target_sequences):
        raise ValueError("MODELLER changed chain inventory")
    expected = set(request.fixed_atoms) | set(request.generated_atoms)
    seen, output, mapping, excluded = set(), [], [], []
    for (native_chain, residues), target in zip(chains.items(), request.target_sequences):
        if len(residues) != len(target.sequence):
            raise ValueError("MODELLER changed residue inventory")
        for index, ((number, code), lines) in enumerate(residues.items()):
            residue = axis[target.chain][index]
            mapping.append({"native": [native_chain, number.strip(), code.strip()],
                            "target": [residue.chain, residue.residue_number, residue.insertion_code]})
            for line in lines:
                if line[17:20].strip() != _CANONICAL_ONE_TO_THREE[target.sequence[index]]:
                    raise ValueError("MODELLER changed target sequence")
                atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, line[12:16].strip())
                if atom in seen:
                    raise ValueError("MODELLER duplicate atom")
                seen.add(atom)
                if atom not in expected:
                    if atom.atom_name != "OXT" or index != len(residues)-1:
                        raise ValueError("MODELLER unexpected atom inventory")
                    excluded.append(line.rstrip())
                    continue
                if not np.isfinite([float(line[i:i+8]) for i in (30, 38, 46)]).all():
                    raise ValueError("MODELLER nonfinite coordinates")
                output.append(line[:21] + residue.chain + f"{int(residue.residue_number):4d}"
                              + (residue.insertion_code or " ") + line[27:])
        output.append("TER\n")
    if not expected <= seen:
        raise ValueError("MODELLER missing requested heavy atoms")
    return "".join(output) + "END\n", {"residue_mapping": mapping, "excluded_native_terminal_OXT": excluded,
        "source_restored": False, "chirality_repaired": False, "reference_coordinates_read": False}
