"""Insert generated-only atoms while retaining source observations verbatim.

This private materializer is not scientific validation or publication admission.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from dvbfixer.model.diffusion.atom_completion_contract import (
    AtomCompletionRequest,
    audit_completion_request_source,
)
from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    DiffusionRequest,
    ResidueIdentity,
)
from dvbfixer.model.diffusion.pdb_materialize import (
    _conect_serials,
    _coordinate_identity,
    _generated_lines,
    _serial,
    materialize_candidate_pdb,
)


def materialize_atom_completion(
    source: bytes,
    request: AtomCompletionRequest,
    atom_axis: tuple[AtomIdentity, ...],
    coordinates: Sequence[Sequence[float]],
) -> bytes:
    """Accept exactly the generated mask, never even an unchanged fixed atom.

    Residue names/elements derive from the validated request, not runner labels.
    Source atom/ANISOU records are byte-preserved, including line endings. Other
    metadata may be reconciled by the existing whole-gap materializer.
    """
    audit_completion_request_source(request, source)
    text = source.decode("utf-8")
    _source_layout(text)
    if len(atom_axis) != len(coordinates) or len(set(atom_axis)) != len(atom_axis):
        raise DiffusionContractError("completion atom-axis lengths or duplicate identities are invalid")
    if set(atom_axis) != set(request.generated_atoms):
        raise DiffusionContractError("completion axis must contain exactly generated atoms; fixed replacements forbidden")
    represented: dict[AtomIdentity, tuple[str, str, tuple[float, float, float]]] = {}
    names = {partial.residue: partial.residue_name for partial in request.masks.partial_residues}
    from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE

    targets = {item.chain: item.sequence for item in request.target_sequences}
    for gap in request.gaps:
        for index, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True):
            names[residue] = _CANONICAL_ONE_TO_THREE[targets[gap.chain][index]]
    for atom, values in zip(atom_axis, coordinates, strict=True):
        if len(values) != 3:
            raise DiffusionContractError("completion coordinate must contain three values")
        xyz = (float(values[0]), float(values[1]), float(values[2]))
        if any(not math.isfinite(value) or len(f"{value:8.3f}") != 8 for value in xyz):
            raise DiffusionContractError("completion coordinate cannot be represented in PDB format")
        number = int(atom.residue_number)
        if (len(f"{number:4d}") != 4 or str(number) != atom.residue_number
                or len(atom.chain) != 1 or len(atom.insertion_code) > 1):
            raise DiffusionContractError("completion residue identity cannot be represented in PDB format")
        represented[atom] = (names[_residue(atom)], atom.atom_name[0], xyz)
    if request.gaps:
        # A legitimate whole-gap-only adapter, never synthetic gaps for partial
        # residues. Only generated coordinates are passed, so fixed updates are
        # impossible despite the legacy materializer supporting that operation.
        gaps_only = DiffusionRequest(
            4, request.normalized_pdb, request.target_sequences, request.sequence_placements,
            request.gaps, request.masks.fixed_atoms, request.masks.gap_generated_atoms,
            (), len(request.seeds), request.seeds,
        )
        axis = request.masks.gap_generated_atoms
        text = materialize_candidate_pdb(text, axis, [represented[a][2] for a in axis],
                                         tuple(represented[a][0] for a in axis), gaps_only)
    lines, used_serials, last_record = _source_layout(text)
    available = (serial for serial in range(1, 100_000) if serial not in used_serials)
    insertions: dict[int, list[str]] = {}
    for partial in request.masks.partial_residues:
        insertions[last_record[partial.residue]] = _generated_lines(partial.missing_atoms, represented, available)
    output: list[str] = []
    for index, line in enumerate(lines):
        if line[:6].strip() == "MASTER":
            # Its atom/record counts describe the incomplete source, not this
            # candidate. Source bytes remain retained by the request digest.
            continue
        output.append(line)
        if index in insertions:
            if not line.endswith("\n"):
                raise DiffusionContractError("cannot append completion atoms to unterminated source record")
            newline = "\r\n" if line.endswith("\r\n") else "\n"
            output.extend(item.rstrip("\n") + newline for item in insertions[index])
    result = "".join(output).encode("utf-8")
    original_observations = [line for line in source.splitlines(keepends=True)
                             if line[:6].strip() in {b"ATOM", b"HETATM", b"ANISOU"}]
    generated = set(request.generated_atoms)
    retained = [line for line in result.splitlines(keepends=True)
                if line[:6].strip() == b"ANISOU" or (
                    line[:6].strip() in {b"ATOM", b"HETATM"}
                    and _coordinate_identity(line.decode("utf-8")) not in generated)]
    if retained != original_observations:
        raise DiffusionContractError("completion changed original source observations")
    return result


def _residue(atom: AtomIdentity) -> ResidueIdentity:
    return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)


def _source_layout(text: str) -> tuple[list[str], set[int], dict[ResidueIdentity, int]]:
    """Refuse ambiguous record layout before assigning any new PDB serials."""
    lines = text.splitlines(keepends=True)
    serials: set[int] = set()
    atom_serials: set[int] = set()
    last: dict[ResidueIdentity, int] = {}
    previous: ResidueIdentity | None = None
    previous_atom: str | None = None
    ended = False
    for index, line in enumerate(lines):
        record = line[:6].strip()
        if ended and record:
            raise DiffusionContractError("completion source contains records after END")
        if record == "END":
            ended = True
        if record in {"MODEL", "ENDMDL"}:
            raise DiffusionContractError("completion requires a selected single-model PDB without MODEL records")
        if record in {"ATOM", "HETATM", "TER"}:
            serial = _serial(line)
            if record != "TER" or line[6:11].strip():
                if serial is None or not 1 <= serial <= 99_999 or serial in serials:
                    raise DiffusionContractError("completion source serials must be valid and unique")
                serials.add(serial)
                if record != "TER":
                    atom_serials.add(serial)
        atom = _coordinate_identity(line)
        if atom is not None:
            residue = _residue(atom)
            if residue != previous and residue in last:
                raise DiffusionContractError("completion source contains a split or repeated residue block")
            previous = residue
            previous_atom = line
            last[residue] = index
        elif record == "ANISOU":
            if previous_atom is None or previous_atom[6:27] != line[6:27]:
                raise DiffusionContractError("completion source ANISOU must immediately follow its matching atom")
            assert previous is not None
            last[previous] = index
            previous_atom = None
        else:
            previous_atom = None
            if record in {"TER", "END"}:
                previous = None
    for line in lines:
        if line[:6].strip() == "CONECT":
            endpoints = _conect_serials(line)
            if endpoints is None or not set(endpoints) <= atom_serials:
                raise DiffusionContractError("completion source has malformed or dangling CONECT endpoints")
    return lines, serials, last
