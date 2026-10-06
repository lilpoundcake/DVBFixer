"""Materialize generated coordinates without discarding source PDB records."""

from __future__ import annotations

import math
from collections.abc import Iterator, Sequence

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionContractError,
    DiffusionRequest,
    GapKind,
    ResidueIdentity,
)


def materialize_candidate_pdb(
    source_text: str,
    atom_axis: tuple[AtomIdentity, ...],
    coordinates: Sequence[Sequence[float]],
    residue_names: tuple[str, ...],
    request: DiffusionRequest,
    *,
    elements: tuple[str, ...] | None = None,
) -> str:
    """Preserve source records while updating fixed and inserting generated atoms."""
    if len(atom_axis) != len(coordinates) or len(atom_axis) != len(residue_names):
        raise DiffusionContractError("candidate atom-axis metadata lengths do not match")
    if elements is not None and len(elements) != len(atom_axis):
        raise DiffusionContractError("candidate element metadata length does not match")
    if len(set(atom_axis)) != len(atom_axis):
        raise DiffusionContractError("candidate atom axis contains duplicate identities")

    represented: dict[AtomIdentity, tuple[str, str, tuple[float, float, float]]] = {}
    for index, (identity, residue_name, xyz_values) in enumerate(
        zip(atom_axis, residue_names, coordinates)
    ):
        if len(xyz_values) != 3:
            raise DiffusionContractError("candidate coordinate must contain three values")
        xyz = (
            float(xyz_values[0]),
            float(xyz_values[1]),
            float(xyz_values[2]),
        )
        if any(not math.isfinite(value) or len(f"{value:8.3f}") != 8 for value in xyz):
            raise DiffusionContractError("candidate coordinate cannot be represented in PDB format")
        element = elements[index].strip().upper() if elements is not None else _element(identity)
        represented[identity] = (residue_name, element, xyz)

    allowed = set(request.fixed_atoms) | set(request.generated_atoms)
    unexpected = set(represented) - allowed
    if unexpected:
        raise DiffusionContractError(
            f"candidate atom axis contains {len(unexpected)} unrequested atoms"
        )
    required = set(request.generated_atoms)
    if not required <= set(represented):
        raise DiffusionContractError(
            f"candidate atom axis omits {len(required - set(represented))} generated atoms"
        )
    source_lines = source_text.splitlines(keepends=True)
    source_atom_serials = {
        serial
        for line in source_lines
        if line[:6].strip() in {"ATOM", "HETATM"}
        for serial in [_serial(line)]
        if serial is not None
    }
    used_serials = {
        serial
        for line in source_lines
        if line[:6].strip() in {"ATOM", "HETATM", "TER"}
        for serial in [_serial(line)]
        if serial is not None
    }
    available_serials = (
        serial for serial in range(1, 100_000) if serial not in used_serials
    )
    gap_by_generated_residue = {
        residue: index
        for index, gap in enumerate(request.gaps)
        for residue in gap.generated_residues
    }
    if len(gap_by_generated_residue) != sum(
        len(gap.generated_residues) for gap in request.gaps
    ):
        raise DiffusionContractError("generated residues overlap between gaps")
    generated_atoms_by_gap: list[list[AtomIdentity]] = [
        [] for _gap in request.gaps
    ]
    for identity in request.generated_atoms:
        residue = ResidueIdentity(
            identity.chain,
            identity.residue_number,
            identity.insertion_code,
        )
        try:
            gap_index = gap_by_generated_residue[residue]
        except KeyError as exc:
            raise DiffusionContractError(
                "generated atom does not belong to a requested gap"
            ) from exc
        generated_atoms_by_gap[gap_index].append(identity)
    generated_lines_by_gap = [
        _generated_lines(tuple(atoms), represented, available_serials)
        for atoms in generated_atoms_by_gap
    ]

    fixed = set(request.fixed_atoms)
    generated = set(request.generated_atoms)
    fixed_seen: set[AtomIdentity] = set()
    inserted: set[int] = set()
    pending_internal: set[int] = set()
    left_anchors: dict[ResidueIdentity, list[int]] = {}
    right_anchors: dict[ResidueIdentity, list[int]] = {}
    n_terminal_anchors: dict[ResidueIdentity, list[int]] = {}
    c_terminal_anchors: dict[ResidueIdentity, list[int]] = {}
    for index, gap in enumerate(request.gaps):
        if gap.gap_kind is GapKind.N_TERMINAL:
            assert gap.right_anchor is not None
            n_terminal_anchors.setdefault(gap.right_anchor, []).append(index)
        elif gap.gap_kind is GapKind.C_TERMINAL:
            assert gap.left_anchor is not None
            c_terminal_anchors.setdefault(gap.left_anchor, []).append(index)
        else:
            assert gap.left_anchor is not None and gap.right_anchor is not None
            left_anchors.setdefault(gap.left_anchor, []).append(index)
            right_anchors.setdefault(gap.right_anchor, []).append(index)
    output: list[str] = []
    previous_residue: ResidueIdentity | None = None

    def insert_c_terminal_after(residue: ResidueIdentity | None) -> None:
        if residue is None:
            return
        for gap_index in c_terminal_anchors.get(residue, ()):
            if gap_index not in inserted:
                output.extend(generated_lines_by_gap[gap_index])
                inserted.add(gap_index)

    for line in source_lines:
        record = line[:6].strip()
        if record == "MASTER":
            continue
        if record == "CONECT":
            endpoints = _conect_serials(line)
            if endpoints is None or not set(endpoints) <= source_atom_serials:
                continue
        if record in {"END", "ENDMDL"}:
            insert_c_terminal_after(previous_residue)
            previous_residue = None
        if record == "TER":
            insert_c_terminal_after(previous_residue)
            previous_residue = None
        if record == "TER" and pending_internal:
            raise DiffusionContractError(
                "source PDB contains TER between the diffusion gap anchors"
            )
        source_identity = _coordinate_identity(line)
        if source_identity is not None:
            residue = ResidueIdentity(
                source_identity.chain,
                source_identity.residue_number,
                source_identity.insertion_code,
            )
            if previous_residue is not None and residue != previous_residue:
                insert_c_terminal_after(previous_residue)
            for gap_index in n_terminal_anchors.get(residue, ()):
                if gap_index not in inserted:
                    output.extend(generated_lines_by_gap[gap_index])
                    inserted.add(gap_index)
            pending_internal.update(left_anchors.get(residue, ()))
            for gap_index in right_anchors.get(residue, ()):
                if gap_index not in inserted:
                    output.extend(generated_lines_by_gap[gap_index])
                    inserted.add(gap_index)
                    pending_internal.discard(gap_index)
            if source_identity in generated:
                raise DiffusionContractError(
                    "source PDB already contains a requested generated atom"
                )
            if source_identity in fixed:
                if source_identity in fixed_seen:
                    raise DiffusionContractError(
                        "source PDB contains a duplicate requested fixed atom"
                    )
                fixed_seen.add(source_identity)
                if source_identity in represented:
                    xyz = represented[source_identity][2]
                    line = (
                        line[:30]
                        + "".join(f"{value:8.3f}" for value in xyz)
                        + line[54:]
                    )
            previous_residue = residue
        output.append(line)

    insert_c_terminal_after(previous_residue)

    if fixed_seen != fixed:
        raise DiffusionContractError(
            f"source PDB omits {len(fixed - fixed_seen)} requested fixed atoms"
        )
    if len(inserted) != len(request.gaps):
        raise DiffusionContractError(
            f"source PDB omits {len(request.gaps) - len(inserted)} gap insertion point(s)"
        )
    return "".join(output)


def _generated_lines(
    generated_atoms: tuple[AtomIdentity, ...],
    represented: dict[AtomIdentity, tuple[str, str, tuple[float, float, float]]],
    available_serials: Iterator[int],
) -> list[str]:
    lines: list[str] = []
    serials = iter(available_serials)
    for identity in generated_atoms:
        try:
            residue_number = int(identity.residue_number)
        except ValueError as exc:
            raise DiffusionContractError(
                "generated PDB residue numbers must be integers"
            ) from exc
        try:
            serial = next(serials)
        except StopIteration as exc:
            raise DiffusionContractError("PDB atom serial space is exhausted") from exc
        residue_name, element, xyz = represented[identity]
        lines.append(
            f"ATOM  {serial:5d} {_atom_name_field(identity.atom_name, element)} "
            f"{residue_name:>3} {identity.chain}{residue_number:4d}"
            f"{identity.insertion_code or ' ':1}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {element:>2}\n"
        )
    return lines


def _coordinate_identity(line: str) -> AtomIdentity | None:
    if line[:6].strip() not in {"ATOM", "HETATM"} or len(line) < 54:
        return None
    altloc = line[16:17]
    if altloc not in {"", " ", "A"}:
        return None
    return AtomIdentity(
        line[21:22],
        line[22:26].strip(),
        line[26:27].strip(),
        line[12:16].strip(),
    )


def _serial(line: str) -> int | None:
    try:
        return int(line[6:11])
    except ValueError:
        return None


def _conect_serials(line: str) -> tuple[int, ...] | None:
    fields = [line[index:index + 5] for index in range(6, len(line.rstrip("\r\n")), 5)]
    try:
        serials = tuple(int(field) for field in fields if field.strip())
    except ValueError:
        return None
    return serials if len(serials) >= 2 else None


def _element(identity: AtomIdentity) -> str:
    name = identity.atom_name.lstrip("0123456789")
    if not name:
        raise DiffusionContractError("cannot infer element from atom name")
    return name[0].upper()


def _atom_name_field(atom_name: str, element: str) -> str:
    if len(atom_name) == 4 or (atom_name and atom_name[0].isdigit()) or len(element) == 2:
        return f"{atom_name:<4}"
    return f" {atom_name:<3}"
