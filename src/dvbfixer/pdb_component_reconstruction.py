"""Whole-PDB adapter for fixed-anchor Class A component reconstruction."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from dvbfixer.ccd import CcdComponent, ccd_to_heavy_atom_graph
from dvbfixer.domain.nonprotein_reconstruction import (
    ComponentClass,
    CoordinateSource,
    ExternalLink,
    ObservedAtom,
    ReconstructionRequest,
    ReconstructionResult,
    ReconstructionStatus,
    UserMappedAuthority,
    ValidationFinding,
)
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.nonprotein_reconstruction import chemical_graph_digest, reconstruct_nonprotein

MAX_PDB_BYTES = 100 * 1024 * 1024
_SAFE_BUNDLE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}\Z")


class PdbComponentError(ValueError):
    """A PDB cannot be mapped or materialized without guessing identity."""


@dataclass(frozen=True, slots=True)
class PdbAtom:
    line_index: int
    model: int
    serial: int
    record: str
    name: str
    alternate_location: str
    residue_name: str
    chain_id: str
    sequence_number: str
    insertion_code: str
    element: str
    position: tuple[float, float, float]
    occurrence: int


@dataclass(frozen=True, slots=True)
class PdbReconstructionBundle:
    result: ReconstructionResult
    input_sha256: str
    output_sha256: str | None
    added_atom_names: tuple[str, ...]
    structure_name: str


def _json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _json_value(asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


def _parse_atom_lines(lines: list[str]) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    model = 1
    previous_key: tuple[int, str, str, str, str] | None = None
    occurrence_by_key: dict[tuple[int, str, str, str, str], int] = {}
    active_occurrence = 1
    for line_index, line in enumerate(lines):
        if line.startswith("MODEL "):
            try:
                model = int(line[10:14].strip())
            except ValueError as exc:
                raise PdbComponentError("MODEL records must contain an integer serial") from exc
            previous_key = None
            continue
        if line.startswith("TER"):
            previous_key = None
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line.rstrip("\r\n")) < 78:
            raise PdbComponentError(f"truncated PDB coordinate record on line {line_index + 1}")
        try:
            serial = int(line[6:11])
            position = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError as exc:
            raise PdbComponentError(f"invalid PDB coordinate record on line {line_index + 1}") from exc
        if not all(math.isfinite(value) for value in position):
            raise PdbComponentError(f"non-finite coordinate on line {line_index + 1}")
        if not 1 <= serial <= 99999:
            raise PdbComponentError(f"atom serial is outside PDB limits on line {line_index + 1}")
        key = (model, line[21], line[22:26].strip(), line[26], line[17:20].strip())
        if key != previous_key:
            active_occurrence = occurrence_by_key.get(key, 0) + 1
            occurrence_by_key[key] = active_occurrence
        previous_key = key
        element = line[76:78].strip().upper()
        if not element:
            raise PdbComponentError(f"atom on line {line_index + 1} lacks an explicit element")
        atoms.append(
            PdbAtom(
                line_index,
                model,
                serial,
                line[:6],
                line[12:16].strip(),
                line[16].strip(),
                line[17:20].strip(),
                line[21],
                line[22:26].strip(),
                line[26].strip(),
                element,
                position,
                active_occurrence,
            )
        )
    serials = [atom.serial for atom in atoms]
    if len({atom.model for atom in atoms}) > 1:
        raise PdbComponentError("multi-model PDB reconstruction is unsupported")
    if len(serials) != len(set(serials)):
        raise PdbComponentError("PDB atom serials must be unique")
    return atoms


def _component_ref(atom: PdbAtom) -> ComponentInstanceRef:
    return ComponentInstanceRef(
        atom.model,
        atom.chain_id,
        atom.sequence_number,
        atom.insertion_code,
        atom.alternate_location,
        atom.occurrence,
    )


def _selected_atoms(
    atoms: list[PdbAtom],
    component: CcdComponent,
    *,
    model: int,
    chain_id: str,
    sequence_number: str,
    insertion_code: str,
    occurrence: int,
) -> list[PdbAtom]:
    selected = [
        atom
        for atom in atoms
        if atom.model == model
        and atom.chain_id == chain_id
        and atom.sequence_number == sequence_number
        and atom.insertion_code == insertion_code
        and atom.residue_name == component.component_id
        and atom.occurrence == occurrence
    ]
    if not selected:
        raise PdbComponentError(
            "no exact component instance matches model/chain/residue/insertion-code/occurrence"
        )
    if any(atom.record != "HETATM" for atom in selected):
        raise PdbComponentError("component reconstruction accepts HETATM residues only")
    if any(atom.alternate_location for atom in selected):
        raise PdbComponentError("alternate-location component reconstruction is unsupported")
    return selected


def _external_links(
    lines: list[str], atoms: list[PdbAtom], selected: list[PdbAtom]
) -> tuple[ExternalLink, ...]:
    selected_serials = {atom.serial for atom in selected}
    by_serial = {atom.serial: atom for atom in atoms}
    selected_identity = _component_ref(selected[0])
    links: dict[tuple[str, int], ExternalLink] = {}
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        values: list[int] = []
        for start in range(6, len(line.rstrip("\r\n")), 5):
            field = line[start : start + 5].strip()
            if not field:
                continue
            try:
                values.append(int(field))
            except ValueError as exc:
                raise PdbComponentError("invalid CONECT atom serial") from exc
        if len(values) < 2:
            continue
        source = values[0]
        for target in values[1:]:
            selected_serial = source if source in selected_serials else target if target in selected_serials else None
            other_serial = target if selected_serial == source else source
            if selected_serial is None or other_serial in selected_serials:
                continue
            selected_atom = by_serial.get(selected_serial)
            other_atom = by_serial.get(other_serial)
            if selected_atom is None or other_atom is None:
                raise PdbComponentError("selected component has a dangling external CONECT record")
            links[(selected_atom.name, other_serial)] = ExternalLink(
                ExactAtomRef(selected_identity, selected_atom.name),
                ExactAtomRef(_component_ref(other_atom), other_atom.name),
            )
    selected_key = (
        selected[0].chain_id,
        selected[0].sequence_number,
        selected[0].insertion_code,
        selected[0].residue_name,
    )
    for line in lines:
        if not line.startswith("LINK  "):
            continue
        if len(line.rstrip("\r\n")) < 57:
            raise PdbComponentError("truncated LINK record")
        endpoints = (
            (line[12:16].strip(), line[16].strip(), line[17:20].strip(), line[21], line[22:26].strip(), line[26].strip()),
            (line[42:46].strip(), line[46].strip(), line[47:50].strip(), line[51], line[52:56].strip(), line[56].strip()),
        )
        matches = [
            (chain, sequence, icode, residue) == selected_key
            for _name, _altloc, residue, chain, sequence, icode in endpoints
        ]
        if matches == [False, False] or matches == [True, True]:
            continue
        selected_endpoint = endpoints[0] if matches[0] else endpoints[1]
        other_endpoint = endpoints[1] if matches[0] else endpoints[0]
        selected_name, selected_altloc, _residue, _chain, _sequence, _icode = selected_endpoint
        if selected_altloc:
            raise PdbComponentError("alternate-location LINK records are unsupported")
        other_name, other_altloc, other_residue, other_chain, other_sequence, other_icode = other_endpoint
        candidates = [
            atom
            for atom in atoms
            if atom.name == other_name
            and atom.alternate_location == other_altloc
            and atom.residue_name == other_residue
            and atom.chain_id == other_chain
            and atom.sequence_number == other_sequence
            and atom.insertion_code == other_icode
            and atom.model == selected[0].model
        ]
        if len(candidates) != 1:
            raise PdbComponentError("LINK endpoint cannot be mapped to one exact atom")
        other_atom = candidates[0]
        links[(selected_name, other_atom.serial)] = ExternalLink(
            ExactAtomRef(selected_identity, selected_name),
            ExactAtomRef(_component_ref(other_atom), other_atom.name),
        )
    return tuple(links[key] for key in sorted(links))


def build_reconstruction_request(
    pdb_content: bytes,
    component: CcdComponent,
    *,
    model: int,
    chain_id: str,
    sequence_number: str,
    insertion_code: str = "",
    occurrence: int = 1,
    ph: float = 7.0,
) -> tuple[ReconstructionRequest, list[str], list[PdbAtom], list[PdbAtom]]:
    """Map one exact PDB HETATM instance and pinned CCD record into the core service."""
    if not pdb_content or len(pdb_content) > MAX_PDB_BYTES:
        raise PdbComponentError("PDB input is empty or exceeds the 100 MiB limit")
    if component.provenance.mode != "local-cif":
        raise PdbComponentError(
            "reconstruction requires an explicitly selected local CCD CIF; online lookup is read-only"
        )
    try:
        text = pdb_content.decode("ascii")
    except UnicodeDecodeError as exc:
        raise PdbComponentError("PDB input must be ASCII") from exc
    lines = text.splitlines(keepends=True)
    atoms = _parse_atom_lines(lines)
    selected = _selected_atoms(
        atoms,
        component,
        model=model,
        chain_id=chain_id,
        sequence_number=sequence_number,
        insertion_code=insertion_code,
        occurrence=occurrence,
    )
    component_ref = _component_ref(selected[0])
    heavy_observed = tuple(
        ObservedAtom(ExactAtomRef(component_ref, atom.name), atom.element, atom.position)
        for atom in selected
        if atom.element != "H"
    )
    graph = ccd_to_heavy_atom_graph(component)
    authority = UserMappedAuthority(
        graph,
        f"wwPDB CCD CIF {component.component_id} sha256:{component.provenance.sha256}",
        chemical_graph_digest(graph),
    )
    request = ReconstructionRequest(
        component_ref,
        heavy_observed,
        authority,
        ph,
        external_links=_external_links(lines, atoms, selected),
    )
    return request, lines, atoms, selected


def _format_atom_name(name: str, element: str) -> str:
    if not 1 <= len(name) <= 4 or not name.isascii():
        raise PdbComponentError(f"atom name {name!r} cannot be represented in PDB")
    if len(name) == 4:
        return name
    if len(element) == 1 and not name[0].isdigit():
        return f" {name:<3}"
    return f"{name:<4}"


def _format_generated_atom(
    *,
    serial: int,
    name: str,
    element: str,
    position: tuple[float, float, float],
    residue_name: str,
    chain_id: str,
    sequence_number: str,
    insertion_code: str,
) -> str:
    if serial > 99999 or len(residue_name) > 3 or len(chain_id) != 1:
        raise PdbComponentError("generated identity exceeds fixed-column PDB limits")
    try:
        residue_integer = int(sequence_number)
    except ValueError as exc:
        raise PdbComponentError("PDB residue number must be an integer") from exc
    if not -999 <= residue_integer <= 9999 or len(insertion_code) > 1:
        raise PdbComponentError("generated residue identity exceeds fixed-column PDB limits")
    if len(element) not in {1, 2}:
        raise PdbComponentError(f"element {element!r} cannot be represented in PDB")
    coordinate_fields = [f"{value:8.3f}" for value in position]
    if any(len(field) != 8 for field in coordinate_fields):
        raise PdbComponentError("generated coordinate exceeds fixed-column PDB limits")
    atom_field = _format_atom_name(name, element)
    return (
        f"HETATM{serial:5d} {atom_field} {residue_name:>3} {chain_id}"
        f"{residue_integer:4d}{insertion_code or ' '}   {''.join(coordinate_fields)}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _validate_environment(
    result: ReconstructionResult, atoms: list[PdbAtom], selected: list[PdbAtom]
) -> None:
    selected_serials = {atom.serial for atom in selected}
    environment = [atom for atom in atoms if atom.serial not in selected_serials]
    generated = [
        coordinate
        for coordinate in result.coordinates
        if coordinate.source is CoordinateSource.SOURCE_TEMPLATE
    ]
    for coordinate in generated:
        for atom in environment:
            distance = math.dist(coordinate.position, atom.position)
            if distance < 0.70:
                raise PdbComponentError(
                    f"generated atom {coordinate.identity.name} clashes with "
                    f"{atom.chain_id}:{atom.sequence_number}:{atom.name} at {distance:.3f} A"
                )


def materialize_reconstruction(
    lines: list[str],
    atoms: list[PdbAtom],
    selected: list[PdbAtom],
    component: CcdComponent,
    result: ReconstructionResult,
) -> tuple[bytes, tuple[str, ...]]:
    """Add only generated heavy atoms; every source coordinate line remains byte-identical."""
    if result.status is not ReconstructionStatus.SUCCEEDED or result.resolved_graph is None:
        raise PdbComponentError("only successful reconstruction geometry can be materialized")
    _validate_environment(result, atoms, selected)
    generated = [
        coordinate
        for coordinate in result.coordinates
        if coordinate.source is CoordinateSource.SOURCE_TEMPLATE
    ]
    if not generated:
        return "".join(lines).encode("ascii"), ()
    next_serial = max((atom.serial for atom in atoms), default=0) + 1
    serial_by_name = {atom.name: atom.serial for atom in selected}
    generated_lines: list[str] = []
    for coordinate in generated:
        serial_by_name[coordinate.identity.name] = next_serial
        generated_lines.append(
            _format_generated_atom(
                serial=next_serial,
                name=coordinate.identity.name,
                element=coordinate.element,
                position=coordinate.position,
                residue_name=component.component_id,
                chain_id=result.component.chain_id,
                sequence_number=result.component.sequence_number,
                insertion_code=result.component.insertion_code,
            )
        )
        next_serial += 1

    insert_at = max(atom.line_index for atom in selected) + 1
    while insert_at < len(lines) and lines[insert_at].startswith("ANISOU"):
        insert_at += 1
    if insert_at and not lines[insert_at - 1].endswith(("\n", "\r")):
        raise PdbComponentError("cannot insert atoms after a source record without a line ending")
    output = lines[:insert_at] + generated_lines + lines[insert_at:]
    generated_names = {coordinate.identity.name for coordinate in generated}
    conect: list[str] = []
    adjacency: dict[str, list[str]] = {name: [] for name in generated_names}
    for bond in result.resolved_graph.bonds:
        if bond.atom1 in generated_names:
            adjacency[bond.atom1].append(bond.atom2)
        if bond.atom2 in generated_names:
            adjacency[bond.atom2].append(bond.atom1)
    for name in sorted(adjacency, key=lambda item: serial_by_name[item]):
        neighbors = sorted(
            {serial_by_name[item] for item in adjacency[name]},
        )
        for start in range(0, len(neighbors), 4):
            values = [serial_by_name[name], *neighbors[start : start + 4]]
            conect.append("CONECT" + "".join(f"{value:5d}" for value in values) + "\n")
    end_index = next(
        (index for index in range(len(output) - 1, -1, -1) if output[index].startswith("END")),
        len(output),
    )
    output[end_index:end_index] = conect
    return "".join(output).encode("ascii"), tuple(coordinate.identity.name for coordinate in generated)


def reconstruct_pdb_component(
    pdb_content: bytes,
    component: CcdComponent,
    *,
    model: int,
    chain_id: str,
    sequence_number: str,
    insertion_code: str = "",
    occurrence: int = 1,
    ph: float = 7.0,
) -> tuple[PdbReconstructionBundle, bytes | None]:
    """Reconstruct one exact instance and return bytes only after all validation passes."""
    if component.provenance.mode != "local-cif":
        raise PdbComponentError(
            "reconstruction requires an explicitly selected local CCD CIF; online lookup is read-only"
        )
    if not pdb_content or len(pdb_content) > MAX_PDB_BYTES:
        raise PdbComponentError("PDB input is empty or exceeds the 100 MiB limit")
    provisional_class = component.provisional_class
    if provisional_class in {"C", "D"}:
        component_class = ComponentClass(provisional_class)
        ambiguous = provisional_class == "D" and component.ambiguous
        status = ReconstructionStatus.AMBIGUOUS if ambiguous else ReconstructionStatus.UNSUPPORTED
        code = "ambiguous-authority" if ambiguous else f"class-{provisional_class.lower()}-unsupported"
        message = (
            "CCD marks this component as chemically ambiguous"
            if ambiguous
            else "component chemistry is outside the isolated organic Class A reconstruction contract"
        )
        result = ReconstructionResult(
            status,
            component_class,
            ComponentInstanceRef(
                model,
                chain_id,
                sequence_number,
                insertion_code,
                "",
                occurrence,
            ),
            findings=(ValidationFinding(code, message),),
        )
        return (
            PdbReconstructionBundle(
                result,
                hashlib.sha256(pdb_content).hexdigest(),
                None,
                (),
                "structure.pdb",
            ),
            None,
        )
    request, lines, atoms, selected = build_reconstruction_request(
        pdb_content,
        component,
        model=model,
        chain_id=chain_id,
        sequence_number=sequence_number,
        insertion_code=insertion_code,
        occurrence=occurrence,
        ph=ph,
    )
    result = reconstruct_nonprotein(request)
    if result.status is not ReconstructionStatus.SUCCEEDED:
        return (
            PdbReconstructionBundle(
                result,
                hashlib.sha256(pdb_content).hexdigest(),
                None,
                (),
                "structure.pdb",
            ),
            None,
        )
    output, added = materialize_reconstruction(lines, atoms, selected, component, result)
    return (
        PdbReconstructionBundle(
            result,
            hashlib.sha256(pdb_content).hexdigest(),
            hashlib.sha256(output).hexdigest(),
            added,
            "structure.pdb",
        ),
        output,
    )


def publish_pdb_bundle(
    bundle: PdbReconstructionBundle,
    output: bytes,
    root: Path,
    name: str,
) -> Path:
    """Atomically publish structure and provenance as one new directory."""
    if bundle.result.status is not ReconstructionStatus.SUCCEEDED:
        raise PdbComponentError("a non-successful reconstruction cannot be published")
    if bundle.output_sha256 != hashlib.sha256(output).hexdigest():
        raise PdbComponentError("materialized structure digest does not match provenance")
    root = root.expanduser()
    if root.is_symlink():
        raise PdbComponentError("publication root must be a non-symlink directory")
    root = root.resolve(strict=True)
    if not root.is_dir():
        raise PdbComponentError("publication root must be a non-symlink directory")
    if _SAFE_BUNDLE_NAME.fullmatch(name) is None:
        raise PdbComponentError("bundle name must be one portable path component")
    destination = root / name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"publication destination already exists: {destination}")
    staging = Path(tempfile.mkdtemp(prefix=".dvbfixer-component-", dir=root))
    committed = False
    try:
        structure_path = staging / bundle.structure_name
        report_path = staging / "provenance.json"
        structure_path.write_bytes(output)
        report_path.write_text(
            json.dumps(_json_value(bundle), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        for path in (structure_path, report_path):
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
        directory_fd = os.open(staging, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        os.replace(staging, destination)
        committed = True
        root_fd = os.open(root, os.O_RDONLY)
        try:
            os.fsync(root_fd)
        finally:
            os.close(root_fd)
        return destination
    except Exception:
        if committed:
            shutil.rmtree(destination, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)
