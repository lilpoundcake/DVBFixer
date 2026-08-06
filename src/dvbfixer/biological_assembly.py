"""Parse PDB REMARK 350 records and render biological assemblies."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


class AssemblyError(ValueError):
    """The deposited biological-assembly description cannot be used safely."""


@dataclass
class AssemblyGroup:
    chains: list[str]
    rows: dict[int, dict[int, tuple[float, float, float, float]]] = field(
        default_factory=dict
    )

    def operators(self) -> list[tuple[int, np.ndarray, np.ndarray]]:
        result = []
        for operator_id in sorted(self.rows):
            rows = self.rows[operator_id]
            if set(rows) != {1, 2, 3}:
                raise AssemblyError(f"BIOMT operator {operator_id} is incomplete")
            matrix = np.array([rows[i][:3] for i in (1, 2, 3)], dtype=float)
            offset = np.array([rows[i][3] for i in (1, 2, 3)], dtype=float)
            result.append((operator_id, matrix, offset))
        if not result:
            raise AssemblyError("assembly chain group has no BIOMT operators")
        return result


@dataclass
class BiologicalAssembly:
    identifier: str
    groups: list[AssemblyGroup] = field(default_factory=list)


_BIOMOLECULE = re.compile(r"BIOMOLECULE:\s*(.+)$")
_CHAINS = re.compile(r"(?:APPLY THE FOLLOWING TO CHAINS:|AND CHAINS:)\s*(.+)$")
_BIOMT = re.compile(
    r"BIOMT([123])\s+(\d+)\s+([+-]?\d+(?:\.\d*)?)\s+"
    r"([+-]?\d+(?:\.\d*)?)\s+([+-]?\d+(?:\.\d*)?)\s+"
    r"([+-]?\d+(?:\.\d*)?)"
)


def parse_biological_assemblies(lines: list[str]) -> dict[str, BiologicalAssembly]:
    """Parse REMARK 350 biomolecules, chain groups, and BIOMT matrices."""
    assemblies: dict[str, BiologicalAssembly] = {}
    current_ids: list[str] = []
    current_groups: list[AssemblyGroup] = []
    for line in lines:
        if not line.startswith("REMARK 350"):
            continue
        payload = line[10:].strip()
        match = _BIOMOLECULE.search(payload)
        if match:
            current_ids = [part.strip() for part in match.group(1).split(",") if part.strip()]
            if not current_ids:
                raise AssemblyError("empty REMARK 350 BIOMOLECULE identifier")
            current_groups = []
            for identifier in current_ids:
                assembly = assemblies.setdefault(identifier, BiologicalAssembly(identifier))
                current_groups.append(AssemblyGroup([]))
                assembly.groups.append(current_groups[-1])
            continue
        match = _CHAINS.search(payload)
        if match:
            if not current_ids:
                raise AssemblyError("REMARK 350 chain list appears before BIOMOLECULE")
            chains = [part.strip().rstrip(".") for part in match.group(1).split(",")]
            chains = [chain for chain in chains if chain]
            if payload.startswith("AND CHAINS:") and current_groups:
                for group in current_groups:
                    group.chains.extend(chains)
            else:
                current_groups = []
                for identifier in current_ids:
                    group = AssemblyGroup(list(chains))
                    assemblies[identifier].groups.append(group)
                    current_groups.append(group)
            continue
        match = _BIOMT.search(payload)
        if match:
            if not current_groups or not any(group.chains for group in current_groups):
                raise AssemblyError("BIOMT operator appears before an APPLY TO CHAINS record")
            row = int(match.group(1))
            operator_id = int(match.group(2))
            values = tuple(float(match.group(i)) for i in range(3, 7))
            for group in current_groups:
                group.rows.setdefault(operator_id, {})[row] = values

    # A placeholder group is created at BIOMOLECULE so continuation parsing is
    # straightforward; remove it when the real APPLY group follows.
    for assembly in assemblies.values():
        assembly.groups = [group for group in assembly.groups if group.chains]
        if not assembly.groups:
            raise AssemblyError(f"biomolecule {assembly.identifier} has no chain groups")
        for group in assembly.groups:
            group.operators()
    if not assemblies:
        raise AssemblyError("no biological assemblies found in REMARK 350")
    return assemblies


def _transform_xyz(line: str, matrix: np.ndarray, offset: np.ndarray) -> str:
    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    x, y, z = matrix @ xyz + offset
    return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"


def _transform_anisou(line: str, matrix: np.ndarray) -> str:
    values = [int(line[start:start + 7]) for start in (28, 35, 42, 49, 56, 63)]
    u11, u22, u33, u12, u13, u23 = values
    tensor = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]], dtype=float)
    rotated = matrix @ tensor @ matrix.T
    result = [
        rotated[0, 0], rotated[1, 1], rotated[2, 2],
        rotated[0, 1], rotated[0, 2], rotated[1, 2],
    ]
    return line[:28] + "".join(f"{round(value):7d}" for value in result) + line[70:]


def _set_serial_chain_resid(
    line: str, serial: int, chain: str, residue_number: int | None,
) -> str:
    line = f"{line[:6]}{serial:5d}{line[11:21]}{chain}{line[22:]}"
    if residue_number is not None:
        line = f"{line[:22]}{residue_number:4d} {line[27:]}"
    return line


def _atom_models(lines: list[str]) -> list[tuple[str | None, list[str]]]:
    models: list[tuple[str | None, list[str]]] = []
    current: list[str] = []
    label: str | None = None
    saw_model = False
    for line in lines:
        if line.startswith("MODEL"):
            saw_model = True
            current = []
            label = line[10:14].strip() or str(len(models) + 1)
        elif line.startswith("ENDMDL"):
            models.append((label, current))
            current = []
            label = None
        elif line.startswith(("ATOM  ", "HETATM", "ANISOU", "CONECT")):
            current.append(line)
    if not saw_model:
        return [(None, current)]
    return models


def render_biological_assembly(
    lines: list[str], assembly: BiologicalAssembly, chain_ids: list[str], *,
    renumber: bool = False,
) -> list[str]:
    """Return a standalone PDB for one assembly."""
    source_chains = {
        line[21] for line in lines if line.startswith(("ATOM  ", "HETATM"))
    }
    copies: list[tuple[str, int, bytes, bytes, np.ndarray, np.ndarray]] = []
    seen = set()
    for group_index, group in enumerate(assembly.groups):
        missing = sorted(set(group.chains) - source_chains)
        if missing:
            raise AssemblyError(
                f"biomolecule {assembly.identifier} references missing chain(s): "
                + ", ".join(missing)
            )
        for operator_id, matrix, offset in group.operators():
            signature = (operator_id, matrix.tobytes(), offset.tobytes())
            for chain in group.chains:
                key = (chain, signature)
                if key not in seen:
                    seen.add(key)
                    copies.append((
                        chain, group_index, matrix.tobytes(), offset.tobytes(), matrix, offset,
                    ))
    if len(copies) > len(chain_ids):
        raise AssemblyError(
            f"assembly needs {len(copies)} chain IDs; PDB supports only {len(chain_ids)}"
        )

    used: set[str] = set()
    copy_ids: list[str] = []
    for source_chain, _group, _matrix_key, _offset_key, _matrix, _offset in copies:
        if source_chain in chain_ids and source_chain not in used:
            output_chain = source_chain
        else:
            output_chain = next((candidate for candidate in chain_ids if candidate not in used), "")
        if not output_chain:
            raise AssemblyError("no unused PDB chain ID remains for generated copy")
        used.add(output_chain)
        copy_ids.append(output_chain)

    header = []
    for line in lines:
        if line.startswith((
            "ATOM  ", "HETATM", "ANISOU", "TER", "CONECT", "MODEL", "ENDMDL",
            "MASTER", "END",
        )):
            continue
        if line.startswith(("SEQRES", "SSBOND", "LINK ")):
            continue
        header.append(line)
    header.extend([
        f"REMARK 999 DVBfixer BIOLOGICAL ASSEMBLY {assembly.identifier}\n",
        "REMARK 999 GENERATED FROM REMARK 350 BIOMT OPERATORS\n",
    ])

    copy_by_context = {
        ((matrix_key, offset_key), source_chain): output_chain
        for (
            source_chain, _group, matrix_key, offset_key, _matrix, _offset
        ), output_chain in zip(copies, copy_ids)
    }
    for (
        source_chain, _group, _matrix_key, _offset_key, _matrix, _offset
    ), output_chain in zip(copies, copy_ids):
        for line in lines:
            if line.startswith("SEQRES") and line[11] == source_chain:
                header.append(line[:11] + output_chain + line[12:])
    for line in lines:
        if line.startswith("SSBOND"):
            chain_a, chain_b = line[15], line[29]
            for context in sorted({key[0] for key in copy_by_context}):
                mapped_a = copy_by_context.get((context, chain_a))
                mapped_b = copy_by_context.get((context, chain_b))
                if mapped_a and mapped_b:
                    record = line[:15] + mapped_a + line[16:29] + mapped_b + line[30:]
                    header.append(record)
        elif line.startswith("LINK  "):
            chain_a, chain_b = line[21], line[51]
            for context in sorted({key[0] for key in copy_by_context}):
                mapped_a = copy_by_context.get((context, chain_a))
                mapped_b = copy_by_context.get((context, chain_b))
                if mapped_a and mapped_b:
                    record = line[:21] + mapped_a + line[22:51] + mapped_b + line[52:]
                    header.append(record)

    output = header
    for label, records in _atom_models(lines):
        if label is not None:
            output.append(f"MODEL     {int(label):4d}\n" if label.isdigit() else f"MODEL     {label}\n")
        atoms = [line for line in records if line.startswith(("ATOM  ", "HETATM"))]
        anisou = {
            int(line[6:11]): line for line in records if line.startswith("ANISOU")
        }
        conect = [line for line in lines if line.startswith("CONECT")]
        serial_map: dict[tuple[int, int], int] = {}
        next_serial = 1
        for copy_index, (
            (source_chain, _group, _matrix_key, _offset_key, matrix, offset), output_chain,
        ) in enumerate(
            zip(copies, copy_ids)
        ):
            residue_map: dict[tuple[str, str], int] = {}
            last_atom: str | None = None
            for line in atoms:
                if line[21] != source_chain:
                    continue
                source_serial = int(line[6:11])
                residue_number = None
                if renumber:
                    residue_key = (line[22:26], line[26])
                    residue_number = residue_map.setdefault(residue_key, len(residue_map) + 1)
                transformed = _transform_xyz(line, matrix, offset)
                transformed = _set_serial_chain_resid(
                    transformed, next_serial, output_chain, residue_number
                )
                output.append(transformed)
                serial_map[(copy_index, source_serial)] = next_serial
                if source_serial in anisou:
                    transformed_anisou = _transform_anisou(anisou[source_serial], matrix)
                    output.append(_set_serial_chain_resid(
                        transformed_anisou, next_serial, output_chain, residue_number
                    ))
                last_atom = transformed
                next_serial += 1
            if last_atom is not None:
                output.append(
                    f"TER   {next_serial:5d}      {last_atom[17:20]} {output_chain}"
                    f"{last_atom[22:26]}{last_atom[26]}\n"
                )
                next_serial += 1
        for copy_index in range(len(copies)):
            for line in conect:
                values = [int(value) for value in line[6:].split()]
                mapped = [serial_map.get((copy_index, value)) for value in values]
                if mapped and mapped[0] is not None:
                    bonded = [value for value in mapped[1:] if value is not None]
                    if bonded:
                        output.append("CONECT" + "".join(f"{value:5d}" for value in [mapped[0], *bonded]) + "\n")
        if label is not None:
            output.append("ENDMDL\n")
    output.append("END\n")
    return output


def assembly_output_path(input_path: Path, output: str | None, identifier: str, all_mode: bool) -> Path:
    if output and not all_mode:
        return Path(output)
    base = Path(output) if output else input_path.with_suffix("")
    if base.suffix.lower() == ".pdb":
        base = base.with_suffix("")
    return base.with_name(f"{base.name}_assembly_{identifier}.pdb")
