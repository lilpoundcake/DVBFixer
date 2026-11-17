"""Private complete-context export, not a public diffusion bundle or new request.

Original full-context lines remain an ordered byte-identical subsequence. Only
new-atom serials are reassigned to avoid collisions with original waters/TERs;
new target S-S annotations are explicitly distinguished from original records.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    Pair,
    build_ledger,
    digest,
    read_records,
)
from dvbfixer.model.diffusion.rfd3_refinement import _candidate
from dvbfixer.model.diffusion.rfd3_water_context import _water_records


@dataclass(frozen=True)
class FullContextExport:
    text: str
    dat: DatRecord
    new_atom_serial_mapping: tuple[tuple[AtomIdentity, int, int], ...]
    added_target_pairs: tuple[Pair, ...]
    original_protein_atoms: int
    original_water_atoms: int
    scientific_qualification: bool = False


def target_ssbond_record(number: int, pair: Pair, distance_angstrom: float) -> str:
    """A derived target annotation, never overwrite an original deposited line."""
    if not 1 <= number <= 999 or not math.isfinite(distance_angstrom) or not 0 < distance_angstrom < 100:
        raise CustodyError("target S-S annotation exceeds PDB capacity")
    line = list(" " * 80)
    line[:6] = "SSBOND"
    line[7:10] = f"{number:3d}"
    for atom, chain_column, residue_column, code_column, name_column in (
        (pair[0], 15, 17, 21, 11), (pair[1], 29, 31, 35, 25)
    ):
        if (len(atom.chain) != 1 or len(atom.insertion_code) > 1 or atom.atom_name != "SG"
                or len(f"{int(atom.residue_number):4d}") != 4):
            raise CustodyError("target S-S identity exceeds fixed-column PDB limits")
        line[name_column:name_column + 3] = "CYS"
        line[chain_column] = atom.chain
        line[residue_column:residue_column + 4] = f"{int(atom.residue_number):4d}"
        line[code_column] = atom.insertion_code or " "
    line[59:65], line[66:72], line[73:78] = f"{'1555':>6}", f"{'1555':>6}", f"{distance_angstrom:5.2f}"
    return "".join(line).rstrip() + "\n"


def materialize_full_context(
    request: DiffusionRequest, source: bytes, context: bytes, projection: bytes,
    *, context_sha256: str, projection_sha256: str, target_pairs: tuple[Pair, ...],
) -> FullContextExport:
    """Merge already-audited heavy coordinates into unchanged original water context."""
    if digest(context) != context_sha256 or digest(projection) != projection_sha256:
        raise CustodyError("full-context export input hashes changed")
    ledger = build_ledger(request, source, context)
    _candidate(request, ledger, projection)
    protein, _, _, _ = read_records(projection.decode())
    waters, _ = _water_records(context)
    if not context.endswith(b"END\n"):
        raise CustodyError("full-context export requires a terminal END newline")
    lines = context.decode().splitlines(keepends=True)
    if any(line[:6].strip() in {"MODEL", "ENDMDL", "CONECT", "MASTER", "ANISOU"} for line in lines):
        raise CustodyError("this private export does not support model/conect/master/anisou context")
    generated = set(request.generated_atoms)
    current = {a.identity: a for a in protein}
    normalized = tuple(tuple(sorted(p)) for p in target_pairs)
    endpoints = [a for pair in normalized for a in pair]
    sg = {a.identity for a in protein if a.identity.atom_name == "SG"}
    if (any(len(pair) != 2 or pair[0] == pair[1] for pair in normalized)
            or len(set(normalized)) != len(normalized) or len(set(endpoints)) != len(endpoints)
            or set(endpoints) != sg or not set(ledger.links) <= set(normalized)
            or any(current[a].residue_name != "CYS" for a in endpoints)):
        raise CustodyError("full-context target graph is not exact singly partnered canonical CYS")
    extra_pairs = tuple(sorted(set(normalized) - set(ledger.links)))
    reserved = {int(line[6:11]) for line in lines
                if line[:6].strip() in {"ATOM", "HETATM", "TER"} and line[6:11].strip()}
    start = max(reserved, default=0) + 1
    if start + len(generated) - 1 > 99999:
        raise CustodyError("full-context new atom serial capacity exceeded")
    serials, mapping = {}, []
    for serial, atom in enumerate((a for a in protein if a.identity in generated), start):
        serials[atom.identity] = serial
        mapping.append((atom.identity, atom.serial, serial))
    if serials.keys() != generated:
        raise CustodyError("full-context new atom inventory changed")
    last_observed = {}
    for atom in ledger.atoms:
        last_observed[atom.identity.chain] = atom.identity
    original_by_record = {atom.record: atom.identity for atom in ledger.atoms}
    water_by_record = {water.record for water in waters}
    queue, index = protein, 0
    output = []
    remarks = (
        "REMARK 999 PRIVATE DVBFIXER HEAVY-ATOM RESTORATION CANDIDATE\n"
        f"REMARK 999 RESTORED MODELLED HEAVY ATOMS: {len(generated)}; NOT EXPERIMENTAL DATA\n"
        "REMARK 999 ORIGINAL PROTEIN/WATER RECORDS UNCHANGED; NO H/PROTONATION ADDED\n"
        f"REMARK 999 ADDED NATIVE-TARGET SSBOND RECORDS: {len(extra_pairs)}; NOT SOURCE ANNOTATIONS\n"
        "REMARK 999 PRIVATE CASE GEOMETRY ONLY; NOT A QUALIFIED DIFFUSION BUNDLE\n"
    )
    ss_lines = [i for i, line in enumerate(lines) if line.startswith("SSBOND")]
    annotation_index = ss_lines[-1] if ss_lines else next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM"))) - 1
    next_ss = max((int(lines[i][7:10]) for i in ss_lines), default=0) + 1
    additions = [target_ssbond_record(number, pair, math.dist(current[pair[0]].xyz, current[pair[1]].xyz))
                 for number, pair in enumerate(extra_pairs, next_ss)]
    written, seen_observed, seen_waters = set(), set(), set()

    def add_new(atom):
        if atom.identity not in generated or atom.identity in written:
            raise CustodyError("invalid new-atom insertion during context merge")
        output.append(atom.record[:6] + f"{serials[atom.identity]:5d}" + atom.record[11:])
        written.add(atom.identity)

    for position, line in enumerate(lines):
        if position == 0 and not line.startswith("HEADER"):
            output.append(remarks)
        identity = original_by_record.get(line)
        if identity is not None:
            while index < len(queue) and queue[index].identity != identity:
                if queue[index].identity.chain != identity.chain:
                    raise CustodyError("context/projection protein order differs")
                add_new(queue[index])
                index += 1
            if index == len(queue) or queue[index].record != line:
                raise CustodyError("original observation not preserved by full-context projection")
            output.append(line)
            index += 1
            seen_observed.add(identity)
            if identity == last_observed[identity.chain]:
                while index < len(queue) and queue[index].identity.chain == identity.chain:
                    add_new(queue[index])
                    index += 1
        else:
            if line.startswith(("ATOM  ", "HETATM")):
                if line not in water_by_record:
                    raise CustodyError("full context contains an unrecognized protein/water record")
                seen_waters.add(line)
            output.append(line)
        if position == 0 and line.startswith("HEADER"):
            output.append(remarks)
        if position == annotation_index:
            output.extend(additions)
    if index != len(queue) or written != generated or seen_observed != set(request.fixed_atoms) or seen_waters != water_by_record:
        raise CustodyError("full-context merge lost observations/waters/new atoms")
    text = "".join(output)
    final_atoms, _, links, _ = read_records(text)
    if len(final_atoms) != len(protein) + len(waters) or set(links) != set(normalized):
        raise CustodyError("full-context merge lost inventory or target links")
    dat = DatRecord(description="Private restored heavy atoms only; not measured observations or MD-ready preparation",
                    target_chains={target.chain: len(target.sequence) for target in request.target_sequences})
    for atom in protein:
        if atom.identity not in generated:
            continue
        a = atom.identity
        dat.added_atoms.append({"chain": a.chain, "resid": a.residue_number, "icode": a.insertion_code,
                                "resname": atom.residue_name, "atom": a.atom_name, "element": atom.element})
        summary = dat.residue_summary.setdefault(f"{a.chain}/{atom.residue_name}{a.residue_number}{a.insertion_code}", {"heavy": 0, "hydrogen": 0})
        summary["heavy"] += 1
    return FullContextExport(text, dat, tuple(mapping), extra_pairs, len(ledger.atoms), len(waters))
