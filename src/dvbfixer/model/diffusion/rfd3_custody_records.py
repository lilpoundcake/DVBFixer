"""Private observed-protein source ledger; no sampler or public capability.

Strict admission wraps the canonical PDB/link readers instead of silently
accepting their deliberately permissive handling of malformed source links.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionContractError, DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import _conect_serials
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.scope import (
    DiffusionScopeError,
    _parse_coordinate_records,
    _parse_explicit_links,
)

REVISION = "private-rfd3-source-link-custody-v1"
Pair = tuple[AtomIdentity, AtomIdentity]


class CustodyError(DiffusionContractError):
    """Failed private custody check, never permission to repair or retry."""


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@dataclass(frozen=True)
class SourceAtom:
    identity: AtomIdentity
    residue_name: str
    element: str
    serial: int
    record: str
    xyz: tuple[float, float, float]


@dataclass(frozen=True)
class SourceLedger:
    revision: str
    source_sha256: str
    canonical_request_sha256: str
    context_sha256: str
    atoms: tuple[SourceAtom, ...]
    link_records: tuple[str, ...]
    links: tuple[Pair, ...]
    auxiliary_records: tuple[str, ...]


def expected_names(request: DiffusionRequest) -> dict[AtomIdentity, str]:
    """Resolve the already-frozen placement, without a new sequence alignment."""
    targets = {target.chain: target.sequence for target in request.target_sequences}
    residues = {}
    positions = {}
    for placement in request.sequence_placements:
        if placement.target_length != len(targets[placement.chain]):
            raise CustodyError("frozen placement target length changed")
        for index, residue in zip(placement.observed_target_indices, placement.observed_residues, strict=True):
            key = residue.chain, residue.residue_number, residue.insertion_code
            if key in residues or (residue.chain, index) in positions:
                raise CustodyError("overlapping frozen source placement")
            residues[key] = _CANONICAL_ONE_TO_THREE[targets[residue.chain][index]]
            positions[residue.chain, index] = key
    for gap in request.gaps:
        for index, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True):
            key = residue.chain, residue.residue_number, residue.insertion_code
            if key in residues or (residue.chain, index) in positions:
                raise CustodyError("overlapping frozen gap placement")
            residues[key] = _CANONICAL_ONE_TO_THREE[targets[residue.chain][index]]
            positions[residue.chain, index] = key
    if set(positions) != {(chain, i) for chain, sequence in targets.items() for i in range(len(sequence))}:
        raise CustodyError("incomplete frozen target placement")
    expected = {AtomIdentity(*key, name): residue_name for key, residue_name in residues.items()
                for name in CANONICAL_HEAVY_ATOMS[residue_name]}
    for atom in request.fixed_atoms:
        if atom.atom_name == "OXT":
            terminal = positions[atom.chain, len(targets[atom.chain]) - 1]
            if (atom.chain, atom.residue_number, atom.insertion_code) != terminal:
                raise CustodyError("preserved-only OXT is not a true target terminus")
            expected[atom] = residues[terminal]
    if set(expected) != set(request.fixed_atoms) | set(request.generated_atoms):
        raise CustodyError("request is not complete canonical protein plus observed terminal OXT")
    observed_residues = {(a.chain, a.residue_number, a.insertion_code) for a in request.fixed_atoms}
    if any((a.chain, a.residue_number, a.insertion_code) in observed_residues for a in request.generated_atoms):
        raise CustodyError("partial-residue completion is outside this private contract")
    return expected


def read_records(text: str) -> tuple[tuple[SourceAtom, ...], tuple[str, ...], tuple[Pair, ...], tuple[str, ...]]:
    """Strict source/final syntax; canonical link semantics, no proximity inference."""
    lines = text.splitlines(keepends=True)
    if any(line[:6].strip() in {"MODEL", "ENDMDL"} for line in lines):
        raise CustodyError("explicit multi-model records are outside this private contract")
    try:
        records, _, altlocs = _parse_coordinate_records(lines)
    except DiffusionScopeError as error:
        raise CustodyError(str(error)) from error
    if altlocs:
        raise CustodyError("alternate locations are outside this private contract")
    by_identity = {record.identity: record for record in records}
    atoms, serials, links, auxiliary = [], set(), [], []
    annotation_pairs: set[tuple[str, Pair]] = set()
    for line in lines:
        kind = line[:6].strip()
        if kind in {"ATOM", "HETATM"}:
            serial = int(line[6:11])
            if serial in serials:
                raise CustodyError("duplicate source/final atom serial")
            serials.add(serial)
            atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
            record = by_identity[atom]
            if float(line[54:60]) != 1.0:
                raise CustodyError("partial occupancy is outside this private contract")
            xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            atoms.append(SourceAtom(atom, record.residue_name, record.element, serial, line, xyz))
        elif kind in {"ANISOU", "SIGATM", "SIGUIJ"}:
            auxiliary.append(line)
    atom_map = {a.identity: a for a in atoms}
    coordinate_lines = [line for line in lines if line[:6].strip() in {"ATOM", "HETATM"}]
    if len(auxiliary) != len(set(auxiliary)):
        raise CustodyError("duplicate atom-associated auxiliary record")
    for line in auxiliary:
        atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if (atom not in atom_map or int(line[6:11]) != atom_map[atom].serial
                or line[17:20].strip() != atom_map[atom].residue_name):
            raise CustodyError("dangling or inconsistent atom-associated auxiliary record")
    for line in lines:
        kind = line[:6].strip()
        if kind not in {"SSBOND", "LINK", "CONECT"}:
            continue
        if line in links:
            raise CustodyError("duplicate explicit link record")
        if kind == "CONECT":
            endpoints = _conect_serials(line)
            if (endpoints is None or len(endpoints) < 2 or len(set(endpoints)) != len(endpoints)
                    or not set(endpoints) <= serials):
                raise CustodyError("dangling, repeated or malformed CONECT endpoint")
        else:
            if any(line[start:stop].strip() not in {"", "1555"} for start, stop in ((59, 65), (66, 72))):
                raise CustodyError("unsupported crystallographic image/operator")
            if kind == "LINK" and (line[16:17].strip() or line[46:47].strip()):
                raise CustodyError("alternate LINK endpoint")
        # The canonical CONECT reader resolves serials from coordinate LINES,
        # not from its parsed-atom argument. SS/LINK use identities directly.
        parsed = _parse_explicit_links([*coordinate_lines, line] if kind == "CONECT" else [line], records)
        expected_count = len(endpoints) - 1 if kind == "CONECT" else 1
        if len(parsed) != expected_count:
            raise CustodyError("unresolved or malformed explicit link")
        for pair in parsed:
            if kind != "CONECT":
                if (kind, pair) in annotation_pairs:
                    raise CustodyError("duplicate annotated endpoint pair")
                annotation_pairs.add((kind, pair))
                fields = ((15, 17, 21, 21, 11, 14, "SG"), (29, 31, 35, 35, 25, 28, "SG")) if kind == "SSBOND" else (
                    (21, 22, 26, 26, 17, 20, line[12:16].strip()), (51, 52, 56, 56, 47, 50, line[42:46].strip()))
                for chain, start, stop, code, ns, ne, name in fields:
                    atom = AtomIdentity(line[chain], line[start:stop].strip(), line[code].strip(), name)
                    if atom_map[atom].residue_name != line[ns:ne].strip():
                        raise CustodyError("explicit link residue name differs from endpoint")
        links.append(line)
    pairs = tuple(sorted(_parse_explicit_links(lines, records)))
    if any(not math.isfinite(value) for atom in atoms for value in atom.xyz):
        raise CustodyError("non-finite source/final coordinate")
    return tuple(atoms), tuple(links), pairs, tuple(auxiliary)


def build_ledger(request: DiffusionRequest, source: bytes, context: bytes) -> SourceLedger:
    if digest(source) != request.normalized_pdb.sha256 or not context:
        raise CustodyError("source binding changed or full context absent")
    names = expected_names(request)
    atoms, records, pairs, auxiliary = read_records(source.decode("utf-8"))
    if {atom.identity for atom in atoms} != set(request.fixed_atoms):
        raise CustodyError("source observations differ from frozen fixed inventory")
    if any(atom.residue_name != names[atom.identity] or atom.element != atom.identity.atom_name[0] for atom in atoms):
        raise CustodyError("source sequence/element differs from frozen canonical target")
    declared = tuple(sorted({tuple(sorted((link.atom1, link.atom2))) for link in request.retained_explicit_links}))
    if len(declared) != len(request.retained_explicit_links) or pairs != declared:
        raise CustodyError("source explicit links differ from frozen request")
    if any(endpoint not in request.fixed_atoms for pair in pairs for endpoint in pair):
        raise CustodyError("link endpoint is not an immutable observation")
    return SourceLedger(REVISION, digest(source), digest(request.to_json().encode()), digest(context), atoms, records, pairs, auxiliary)
