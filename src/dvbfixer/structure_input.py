"""Normalize CIF structure inputs to the PDB representation used internally.

DVBfixer's scientific pipelines intentionally remain PDB based: a substantial
part of their behaviour depends on PDB records such as SEQRES, LINK, SSBOND and
REMARK 350.  This module is the single format boundary.  It converts CIF input
once, validates that it is representable as PDB, and leaves every downstream
pipeline unchanged.
"""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote

CIF_EXTENSIONS = {".cif", ".mmcif"}
PDB_CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
_ORIGINAL_INPUT_ENV = "DVBFIXER_ORIGINAL_STRUCTURE_INPUT"


class StructureInputError(ValueError):
    """A structure cannot be converted to dvbfixer's internal PDB format."""


@dataclass(frozen=True)
class NormalizedStructure:
    original_path: Path
    pdb_path: Path
    chain_map: dict[str, str]
    dialect: str

    @property
    def changed_chain_ids(self) -> dict[str, str]:
        return {old: new for old, new in self.chain_map.items() if old != new}


def is_cif_path(path: str | Path) -> bool:
    return Path(path).suffix.lower() in CIF_EXTENSIONS


def display_structure_path(internal_path: str | Path) -> Path:
    """Return the original CIF path while a normalized command is running."""
    original = os.environ.get(_ORIGINAL_INPUT_ENV)
    return Path(original) if original else Path(internal_path)


def _cif_dialect(path: Path) -> str:
    import gemmi

    try:
        document = gemmi.cif.read(str(path))
    except Exception as exc:
        raise StructureInputError(f"cannot parse CIF {path}: {exc}") from exc
    if len(document) != 1:
        raise StructureInputError(
            f"CIF must contain exactly one data block; found {len(document)} in {path}"
        )
    block = document[0]
    if block.find_values("_atom_site.Cartn_x") or block.find_values("_atom_site_Cartn_x"):
        return "pdbx/mmCIF"
    if block.find_values("_atom_site_fract_x") or block.find_values("_atom_site.fract_x"):
        return "small-molecule CIF"
    raise StructureInputError(
        "unrecognized CIF dialect: expected PDBx/mmCIF Cartesian _atom_site.Cartn_x "
        "or crystallographic _atom_site_fract_x coordinates"
    )


def _chain_mapping(chain_names: Sequence[str]) -> dict[str, str]:
    unique = list(dict.fromkeys(chain_names))
    preserved = {
        chain for chain in unique
        if len(chain) == 1 and chain in PDB_CHAIN_IDS
    }
    available = iter(character for character in PDB_CHAIN_IDS if character not in preserved)
    mapping: dict[str, str] = {}
    for chain in unique:
        if chain in preserved:
            mapping[chain] = chain
            continue
        try:
            mapping[chain] = next(available)
        except StopIteration as exc:
            raise StructureInputError(
                f"CIF contains {len(unique)} chains; PDB supports at most "
                f"{len(PDB_CHAIN_IDS)} unique chain IDs"
            ) from exc
    return mapping


def _validate_pdb_text(text: str, source: Path) -> None:
    atoms = 0
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms += 1
        if len(line) < 78:
            raise StructureInputError(
                f"CIF conversion produced a truncated PDB atom record at line {line_number}"
            )
        serial = line[6:11].strip()
        residue = line[22:26].strip()
        if not serial.isdigit() or int(serial) > 99999:
            raise StructureInputError(f"{source} exceeds the PDB atom-serial limit (99999)")
        if not residue or not residue.lstrip("-").isdigit():
            raise StructureInputError(
                f"residue identifier at converted PDB line {line_number} is not representable"
            )
        if len(line[17:20].strip()) > 3 or len(line[12:16].strip()) > 4:
            raise StructureInputError(
                f"atom or residue name at converted PDB line {line_number} exceeds PDB limits"
            )
    if atoms == 0:
        raise StructureInputError(f"CIF contains no convertible atom sites: {source}")
    if atoms > 99999:
        raise StructureInputError(f"{source} has {atoms} atoms; PDB supports at most 99999")


def _validate_mmcif_for_pdb(structure: Any, source: Path) -> None:
    atom_count = 0
    for model in structure:
        if model.num < 1 or model.num > 9999:
            raise StructureInputError(
                f"model number {model.num} in {source} is outside the PDB range 1..9999"
            )
        for chain in model:
            for residue in chain:
                if len(residue.name) > 3:
                    raise StructureInputError(
                        f"residue name {residue.name!r} in chain {chain.name!r} exceeds 3 PDB characters"
                    )
                if residue.seqid.num < -999 or residue.seqid.num > 9999:
                    raise StructureInputError(
                        f"residue {chain.name}/{residue.seqid} is outside the PDB residue-number range"
                    )
                if len(residue.seqid.icode.strip()) > 1:
                    raise StructureInputError(
                        f"residue {chain.name}/{residue.seqid} has a multi-character insertion code"
                    )
                for atom in residue:
                    atom_count += 1
                    if len(atom.name.strip()) > 4:
                        raise StructureInputError(
                            f"atom name {atom.name!r} in {chain.name}/{residue.seqid} exceeds 4 PDB characters"
                        )
                    if atom.altloc not in {"\x00", " ", "."} and len(atom.altloc) > 1:
                        raise StructureInputError(f"atom {atom.name} has a multi-character alternate-location ID")
                    for coordinate in (atom.pos.x, atom.pos.y, atom.pos.z):
                        if not math.isfinite(coordinate) or coordinate < -999.999 or coordinate > 9999.999:
                            raise StructureInputError(
                                f"coordinate {coordinate:g} for atom {chain.name}/{residue.seqid}/{atom.name} "
                                "is outside the PDB 8.3 field range"
                            )
                    for label, value in (("occupancy", atom.occ), ("B factor", atom.b_iso)):
                        if not math.isfinite(value) or value < -99.99 or value > 999.99:
                            raise StructureInputError(
                                f"{label} {value:g} for atom {chain.name}/{residue.seqid}/{atom.name} "
                                "is outside the PDB 6.2 field range"
                            )
    if atom_count > 99999:
        raise StructureInputError(f"{source} has {atom_count} atoms; PDB supports at most 99999")


def _mapping_remarks(mapping: dict[str, str]) -> list[str]:
    remarks: list[str] = []
    for old, new in mapping.items():
        if old == new:
            continue
        encoded = quote(old or "<blank>", safe="._-")
        prefix = f"REMARK 999 DVBFIXER CIF_CHAIN_MAP {new} "
        width = 80 - len(prefix)
        if len(encoded) <= width:
            remarks.append(prefix + encoded + "\n")
            continue
        chunks = [encoded[index:index + 30] for index in range(0, len(encoded), 30)]
        for index, chunk in enumerate(chunks, 1):
            remarks.append(
                f"REMARK 999 DVBFIXER CIF_CHAIN_MAP_PART {new} {index}/{len(chunks)} {chunk}\n"
            )
    return remarks


def _insert_mapping_remarks(path: Path, mapping: dict[str, str]) -> None:
    remarks = _mapping_remarks(mapping)
    if not remarks or not path.is_file() or path.suffix.lower() != ".pdb":
        return
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    if any("DVBFIXER CIF_CHAIN_MAP" in line for line in lines):
        return
    insert_at = next(
        (index for index, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM", "MODEL "))),
        0,
    )
    path.write_text("".join(lines[:insert_at] + remarks + lines[insert_at:]))


def _convert_mmcif(source: Path, destination: Path) -> dict[str, str]:
    import gemmi

    try:
        structure = gemmi.read_structure(str(source), format=gemmi.CoorFormat.Mmcif)
    except Exception as exc:
        raise StructureInputError(f"cannot read PDBx/mmCIF {source}: {exc}") from exc
    _validate_mmcif_for_pdb(structure, source)
    chains = [chain.name for model in structure for chain in model]
    mapping = _chain_mapping(chains)
    for model in structure:
        for chain in model:
            chain.name = mapping[chain.name]

    for connection in structure.connections:
        connection.partner1.chain_name = mapping.get(
            connection.partner1.chain_name, connection.partner1.chain_name
        )
        connection.partner2.chain_name = mapping.get(
            connection.partner2.chain_name, connection.partner2.chain_name
        )

    # Assembly generators may refer to author chain names.  Gemmi translates
    # mmCIF assembly categories to this representation and its PDB writer then
    # emits the corresponding REMARK 350 records.
    for assembly in structure.assemblies:
        for generator in assembly.generators:
            generator.chains[:] = [mapping.get(chain, chain) for chain in generator.chains]

    options = gemmi.PdbWriteOptions()
    options.seqres_records = True
    options.ssbond_records = True
    options.link_records = True
    options.conect_records = True
    options.numbered_ter = False
    try:
        text = structure.make_pdb_string(options)
    except Exception as exc:
        raise StructureInputError(f"cannot represent {source} as PDB: {exc}") from exc
    _validate_pdb_text(text, source)
    destination.write_text("".join(_mapping_remarks(mapping)) + text)
    return mapping


def _convert_small_cif(source: Path, destination: Path) -> dict[str, str]:
    import gemmi

    try:
        from openbabel import openbabel
    except ImportError as exc:
        raise StructureInputError(
            "small-molecule CIF conversion requires Open Babel; install the documented "
            "dvbfixer environment or openbabel-wheel"
        ) from exc

    conversion = openbabel.OBConversion()
    if not conversion.SetInAndOutFormats("cif", "pdb"):
        raise StructureInputError("Open Babel CIF or PDB format plugin is unavailable")
    conversion.AddOption("B", openbabel.OBConversion.INOPTIONS)
    molecule = openbabel.OBMol()
    if not conversion.ReadFile(molecule, str(source)) or molecule.NumAtoms() == 0:
        raise StructureInputError(f"Open Babel could not read small-molecule CIF {source}")
    if any(bond.GetBondOrder() > 1 for bond in openbabel.OBMolBondIter(molecule)):
        print(
            "WARNING: small-molecule CIF bond orders are reduced to PDB CONECT "
            "connectivity for downstream processing",
            file=sys.stderr,
        )
    if not conversion.WriteFile(molecule, str(destination)):
        raise StructureInputError(f"Open Babel could not convert {source} to PDB")
    conversion.CloseOutFile()
    block = gemmi.cif.read(str(source))[0]
    labels = list(block.find_values("_atom_site_label"))
    if not labels:
        labels = list(block.find_values("_atom_site.label"))
    if labels and len(labels) != molecule.NumAtoms():
        raise StructureInputError(
            "small-molecule CIF atom labels do not match the converted asymmetric unit"
        )
    if any(len(label) > 4 for label in labels):
        offending = next(label for label in labels if len(label) > 4)
        raise StructureInputError(
            f"small-molecule atom label {offending!r} exceeds the 4-character PDB limit"
        )
    lines = destination.read_text(errors="replace").splitlines(keepends=True)
    normalized_lines: list[str] = []
    atom_index = 0
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            line = line.ljust(80)
            if labels:
                line = line[:12] + f"{labels[atom_index]:>4}" + line[16:]
            atom_index += 1
            line = line[:17] + "LIG" + line[20:21] + "A" + line[22:]
            line = line.rstrip() + "\n"
        normalized_lines.append(line)
    text = "".join(normalized_lines)
    destination.write_text(text)
    _validate_pdb_text(text, source)
    return {"A": "A"}


def normalize_structure(source: str | Path, destination: str | Path) -> NormalizedStructure:
    source_path = Path(source).expanduser().resolve()
    destination_path = Path(destination)
    if not source_path.is_file():
        raise StructureInputError(f"structure input does not exist: {source_path}")
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    dialect = _cif_dialect(source_path)
    if dialect == "pdbx/mmCIF":
        mapping = _convert_mmcif(source_path, destination_path)
    else:
        mapping = _convert_small_cif(source_path, destination_path)
    return NormalizedStructure(source_path, destination_path, mapping, dialect)


def _split_template_reference(value: str) -> tuple[Path, str | None] | None:
    path = Path(value)
    if path.is_file() and is_cif_path(path):
        return path, None
    if ":" not in value:
        return None
    candidate, chain = value.rsplit(":", 1)
    path = Path(candidate)
    if path.is_file() and is_cif_path(path):
        return path, chain
    return None


def _output_argument(argv: Sequence[str]) -> Path | None:
    for index, token in enumerate(argv[:-1]):
        if token in {"-o", "--output"}:
            return Path(argv[index + 1]).expanduser()
    return None


def _option_value(argv: Sequence[str], option: str) -> tuple[int, str] | None:
    try:
        index = argv.index(option)
    except ValueError:
        return None
    if index + 1 >= len(argv):
        return None
    return index + 1, argv[index + 1]


def _translate_fasta(source: Path, destination: Path, mapping: dict[str, str]) -> None:
    output: list[str] = []
    for raw in source.read_text().splitlines(keepends=True):
        if not raw.startswith(">"):
            output.append(raw)
            continue
        newline = "\n" if raw.endswith("\n") else ""
        header = raw[1:].rstrip("\r\n")
        identifier, separator, description = header.partition(" ")
        translated = mapping.get(identifier)
        if translated is None and identifier.lower().startswith("chain_"):
            chain = identifier[6:]
            translated = "chain_" + mapping.get(chain, chain)
        if translated is None and "_" in identifier:
            prefix, chain = identifier.rsplit("_", 1)
            translated = prefix + "_" + mapping.get(chain, chain)
        identifier = translated or identifier
        output.append(">" + identifier + (separator + description if separator else "") + newline)
    destination.write_text("".join(output))


def _plan_contains_cif(path: Path) -> bool:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False

    def visit(value: object) -> bool:
        if isinstance(value, str):
            return value.lower().endswith((".cif", ".mmcif"))
        if isinstance(value, list):
            return any(visit(item) for item in value)
        if isinstance(value, dict):
            return any(visit(item) for item in value.values())
        return False

    return visit(payload)


def _translate_value(value: str, mapping: dict[str, str]) -> str:
    if value in mapping:
        return mapping[value]
    # Most residue/atom selectors begin with CHAIN:, while comma-separated
    # lists are used by a few topology/protonation options.
    parts = value.split(",")
    translated: list[str] = []
    for part in parts:
        head, separator, tail = part.partition(":")
        translated.append(mapping.get(head, head) + (separator + tail if separator else ""))
    return ",".join(translated)


def _translate_colon_positions(
    value: str, mapping: dict[str, str], positions: set[int],
) -> str:
    fields = value.split(":")
    for position in positions:
        if position < len(fields):
            fields[position] = mapping.get(fields[position], fields[position])
    return ":".join(fields)


_DEFAULT_FILE_SUFFIX = {
    "split": "_split.pdb", "renumber": "_renum.pdb", "model": "_model.pdb",
    "pull": "_pulled.pdb", "prepare": "_prepared.pdb", "minimize": "_minimized.pdb",
    "protonate": "_prot.pdb", "rename": "_renamed.pdb", "puppet": "_puppet.pdb",
    "conect": "_conect.pdb", "transplant": "_transplant.pdb",
    "zbs": "_zbs.pdb",
}


@contextmanager
def normalized_command_inputs(command: str, argv: Sequence[str]) -> Iterator[list[str]]:
    """Rewrite CIF path arguments to temporary PDBs for one CLI invocation."""
    original = list(argv)
    non_input_path_options = {
        "-o", "--output", "--pdb", "--dat", "--fit-dir", "--gromacs",
        "--postflight-report",
    }
    references = [
        None if index and original[index - 1] in non_input_path_options
        else _split_template_reference(token)
        for index, token in enumerate(original)
    ]
    cif_refs = [(index, ref) for index, ref in enumerate(references) if ref is not None]
    plan_option = _option_value(original, "--template-plan")
    plan_path = Path(plan_option[1]).resolve() if plan_option else None
    has_cif_plan = bool(plan_path and plan_path.is_file() and _plan_contains_cif(plan_path))
    if not cif_refs and not has_cif_plan:
        yield original
        return

    explicit_output = _output_argument(original)
    primary = cif_refs[0][1][0] if cif_refs else plan_path
    assert primary is not None
    if (
        explicit_output is not None
        and explicit_output.suffix.lower() in CIF_EXTENSIONS
        and command in {*_DEFAULT_FILE_SUFFIX, "convert", "transplant"}
    ):
        raise StructureInputError(
            "CIF output is not supported; choose a .pdb output path"
        )
    top_pdb_option = _option_value(original, "--pdb") if command == "top" else None
    if top_pdb_option and Path(top_pdb_option[1]).suffix.lower() in CIF_EXTENSIONS:
        raise StructureInputError("CIF output is not supported; --pdb must name a .pdb file")
    if explicit_output is None and command == "convert":
        suffix = "_charmm.pdb" if "--to-charmm" in original else "_amber.pdb"
        explicit_output = primary.with_name(primary.stem + suffix)
        original.extend(["-o", str(explicit_output)])
    if explicit_output is None and command == "parametrize":
        explicit_output = Path.cwd() / primary.stem
        original.extend(["-o", primary.stem])
    if explicit_output is None and command in _DEFAULT_FILE_SUFFIX:
        explicit_output = primary.with_name(primary.stem + _DEFAULT_FILE_SUFFIX[command])
        original.extend(["-o", str(explicit_output)])
    if command == "minimize" and "--dat" not in original:
        original_dat = primary.with_suffix(".dat")
        if original_dat.is_file():
            original.extend(["--dat", str(original_dat)])
    work_parent = (explicit_output.parent if explicit_output else Path.cwd()).resolve()
    work_parent.mkdir(parents=True, exist_ok=True)
    before_pdbs = {path.resolve() for path in work_parent.glob("*.pdb")}

    with tempfile.TemporaryDirectory(prefix=".dvbfixer_cif_", dir=work_parent) as temp_name:
        temp_dir = Path(temp_name)
        rewritten = list(original)
        combined_mapping: dict[str, str] = {}
        cache: dict[Path, NormalizedStructure] = {}

        def converted(source: Path) -> NormalizedStructure:
            resolved = source.resolve()
            normalized = cache.get(resolved)
            if normalized is not None:
                return normalized
            destination = temp_dir / f"{len(cache):03d}" / f"{source.stem}.pdb"
            destination.parent.mkdir(parents=True, exist_ok=True)
            normalized = normalize_structure(source, destination)
            cache[resolved] = normalized
            for old, new in normalized.chain_map.items():
                previous = combined_mapping.get(old)
                if previous is not None and previous != new:
                    raise StructureInputError(
                        f"chain {old!r} maps inconsistently across CIF inputs: "
                        f"{previous!r} and {new!r}; use explicit one-character chain IDs"
                    )
                combined_mapping[old] = new
            print(
                f"[cif] converted {source} ({normalized.dialect}) to internal PDB",
                file=sys.stderr,
            )
            if normalized.changed_chain_ids:
                rendered = ", ".join(
                    f"{old or '<blank>'!r}->{new!r}"
                    for old, new in normalized.changed_chain_ids.items()
                )
                print(
                    f"WARNING: CIF chain IDs mapped for PDB compatibility: {rendered}",
                    file=sys.stderr,
                )
            return normalized

        for index, (source, requested_chain) in cif_refs:
            normalized = converted(source)
            replacement = str(normalized.pdb_path)
            if requested_chain is not None:
                replacement += ":" + normalized.chain_map.get(requested_chain, requested_chain)
            rewritten[index] = replacement

        if has_cif_plan and plan_path and plan_option:
            payload = json.loads(plan_path.read_text())

            def rewrite_plan(value: object, local_mapping: dict[str, str] | None = None) -> object:
                if isinstance(value, list):
                    return [rewrite_plan(item, local_mapping) for item in value]
                if not isinstance(value, dict):
                    return value
                result = dict(value)
                item_mapping = local_mapping
                for key in ("file", "path", "fittedFile", "structure"):
                    candidate = result.get(key)
                    if not isinstance(candidate, str) or not candidate.lower().endswith((".cif", ".mmcif")):
                        continue
                    source = Path(candidate)
                    if not source.is_absolute():
                        source = plan_path.parent / source
                    normalized = converted(source)
                    result[key] = str(normalized.pdb_path)
                    item_mapping = normalized.chain_map
                if item_mapping:
                    for key in ("chain", "sourceChain"):
                        if isinstance(result.get(key), str):
                            result[key] = item_mapping.get(result[key], result[key])
                return {key: rewrite_plan(item, item_mapping) for key, item in result.items()}

            rewritten_plan = temp_dir / "template-plan.json"
            rewritten_plan.write_text(json.dumps(rewrite_plan(payload), indent=2) + "\n")
            rewritten[plan_option[0]] = str(rewritten_plan)

        fasta_option = _option_value(rewritten, "--fasta")
        if fasta_option and combined_mapping:
            fasta_source = Path(fasta_option[1])
            if fasta_source.is_file():
                rewritten_fasta = temp_dir / f"{fasta_source.stem}.mapped.fasta"
                _translate_fasta(fasta_source, rewritten_fasta, combined_mapping)
                rewritten[fasta_option[0]] = str(rewritten_fasta)

        # Translate structured CHAIN:... selectors.  Bare values are translated
        # only for options whose grammar is explicitly a chain ID; changing all
        # bare tokens would corrupt unrelated values such as ligand --name.
        bare_chain_options = {"--cap-chain"}
        for index, token in enumerate(rewritten):
            if token.startswith("-") or Path(token).exists():
                continue
            previous = rewritten[index - 1] if index else ""
            if previous == "--ss":
                rewritten[index] = _translate_colon_positions(
                    token, combined_mapping, {0, 2},
                )
            elif previous == "--align":
                rewritten[index] = _translate_colon_positions(
                    token, combined_mapping, {0, 1},
                )
            elif ":" in token or ("," in token and any(part.partition(":")[0] in combined_mapping for part in token.split(","))):
                rewritten[index] = _translate_value(token, combined_mapping)
            elif previous in bare_chain_options:
                rewritten[index] = combined_mapping.get(token, token)

        previous_original = os.environ.get(_ORIGINAL_INPUT_ENV)
        os.environ[_ORIGINAL_INPUT_ENV] = str(primary.resolve())
        try:
            yield rewritten
        finally:
            if previous_original is None:
                os.environ.pop(_ORIGINAL_INPUT_ENV, None)
            else:
                os.environ[_ORIGINAL_INPUT_ENV] = previous_original
            if combined_mapping:
                candidates = set(work_parent.glob("*.pdb")) - before_pdbs
                if explicit_output and explicit_output.suffix.lower() == ".pdb":
                    candidates.add(explicit_output)
                if command == "top":
                    pdb_option = _option_value(rewritten, "--pdb")
                    candidates.add(
                        Path(pdb_option[1]) if pdb_option else work_parent / "conf.pdb"
                    )
                for output in candidates:
                    _insert_mapping_remarks(output, combined_mapping)


def run_with_normalized_inputs(
    command: str,
    command_main: Callable[[list[str]], object],
    argv: Sequence[str],
) -> object:
    """Invoke a command through the shared CIF-to-PDB boundary."""
    try:
        with normalized_command_inputs(command, argv) as normalized:
            return command_main(normalized)
    except StructureInputError as exc:
        raise SystemExit(f"CIF input error: {exc}") from exc
