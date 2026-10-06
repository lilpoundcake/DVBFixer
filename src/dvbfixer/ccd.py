"""Bounded wwPDB Chemical Component Dictionary information adapter."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, BinaryIO, cast

from dvbfixer.domain.nonprotein_reconstruction import ChemicalAtom, ChemicalBond, ChemicalGraph
from dvbfixer.domain.parameterization import COMPLEX_COFACTORS

MAX_CCD_BYTES = 5 * 1024 * 1024
CCD_DOWNLOAD_ROOT = "https://files.rcsb.org/ligands/download"
_COMPONENT_ID = re.compile(r"[A-Z0-9]{1,8}\Z")
_METALS = frozenset(
    {"LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "AL", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG"}
)
_CLASS_A_ELEMENTS = frozenset(
    {"H", "B", "C", "N", "O", "F", "SI", "P", "S", "CL", "SE", "BR", "I"}
)


class CcdError(ValueError):
    """CCD authority could not be read or validated safely."""


@dataclass(frozen=True, slots=True)
class CcdAtom:
    name: str
    alternate_name: str
    element: str
    formal_charge: int | None
    aromatic: bool
    leaving: bool
    stereo: str | None
    ideal_position: tuple[float, float, float] | None


@dataclass(frozen=True, slots=True)
class CcdBond:
    atom1: str
    atom2: str
    order: str
    aromatic: bool
    stereo: str | None


@dataclass(frozen=True, slots=True)
class CcdProvenance:
    mode: str
    source: str
    sha256: str
    retrieved_from: str | None = None


@dataclass(frozen=True, slots=True)
class CcdComponent:
    component_id: str
    name: str
    component_type: str
    formula: str
    formal_charge: int | None
    formula_weight: float | None
    initial_date: str | None
    modified_date: str | None
    ambiguous: bool
    atoms: tuple[CcdAtom, ...]
    bonds: tuple[CcdBond, ...]
    provenance: CcdProvenance

    @property
    def heavy_atoms(self) -> tuple[CcdAtom, ...]:
        return tuple(atom for atom in self.atoms if atom.element != "H")

    @property
    def provisional_class(self) -> str:
        """Classify intrinsic CCD chemistry before structure-specific links."""
        elements = {atom.element for atom in self.atoms}
        if self.component_id in COMPLEX_COFACTORS or elements & _METALS:
            return "C"
        if self.ambiguous or "C" not in elements or not elements <= _CLASS_A_ELEMENTS:
            return "D"
        return "A"

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["atom_count"] = len(self.atoms)
        result["heavy_atom_count"] = len(self.heavy_atoms)
        result["hydrogen_atom_count"] = len(self.atoms) - len(self.heavy_atoms)
        result["bond_count"] = len(self.bonds)
        result["provisional_class"] = self.provisional_class
        return result


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise CcdError("CCD download redirects are disabled")


def normalize_component_id(value: str) -> str:
    component_id = value.strip().upper()
    if _COMPONENT_ID.fullmatch(component_id) is None:
        raise CcdError("component ID must contain 1-8 ASCII letters or digits")
    return component_id


def _optional(value: str | None) -> str | None:
    if value is None:
        return None
    import gemmi

    decoded = gemmi.cif.as_string(value)
    return decoded or None


def _required(block: Any, tag: str) -> str:
    value = _optional(block.find_value(tag))
    if value is None:
        raise CcdError(f"CCD record is missing {tag}")
    return value


def _flag(value: str) -> bool:
    if value == "Y":
        return True
    if value == "N":
        return False
    raise CcdError(f"invalid CCD Y/N flag: {value!r}")


def _position(values: tuple[str, str, str]) -> tuple[float, float, float] | None:
    if any(_optional(value) is None for value in values):
        return None
    try:
        result = (float(values[0]), float(values[1]), float(values[2]))
        if not all(math.isfinite(value) for value in result):
            raise ValueError
        return result
    except ValueError as exc:
        raise CcdError("CCD ideal coordinates must be finite numbers") from exc


def parse_ccd_component(
    content: bytes,
    *,
    expected_component_id: str,
    provenance: CcdProvenance,
) -> CcdComponent:
    """Parse one bounded CCD component CIF without interpreting coordinates."""
    if not content or len(content) > MAX_CCD_BYTES:
        raise CcdError("CCD response is empty or exceeds the 5 MiB limit")
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CcdError("CCD response is not UTF-8") from exc
    try:
        import gemmi

        document = gemmi.cif.read_string(text)
        if len(document) != 1:
            raise CcdError("CCD component CIF must contain exactly one data block")
        block = document.sole_block()
    except CcdError:
        raise
    except Exception as exc:
        raise CcdError(f"invalid CCD CIF: {exc}") from exc

    component_id = normalize_component_id(_required(block, "_chem_comp.id"))
    if component_id != normalize_component_id(expected_component_id):
        raise CcdError(
            f"CCD response component {component_id!r} does not match request "
            f"{expected_component_id!r}"
        )
    atom_table = block.find([
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.alt_atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_aromatic_flag",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
        "_chem_comp_atom.pdbx_stereo_config",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
    ])
    if not atom_table:
        raise CcdError("CCD component contains no atoms")
    atoms: list[CcdAtom] = []
    names: set[str] = set()
    for row in atom_table:
        name = _optional(str(row[0])) or ""
        if not name or name in names:
            raise CcdError("CCD atom names must be non-empty and unique")
        names.add(name)
        charge_text = _optional(str(row[3]))
        try:
            charge = int(charge_text) if charge_text is not None else None
        except ValueError as exc:
            raise CcdError(f"CCD atom {name!r} has an invalid formal charge") from exc
        stereo = _optional(str(row[6]))
        atoms.append(CcdAtom(
            name=name,
            alternate_name=_optional(str(row[1])) or "",
            element=(_optional(str(row[2])) or "").upper(),
            formal_charge=charge,
            aromatic=_flag(_optional(str(row[4])) or ""),
            leaving=_flag(_optional(str(row[5])) or ""),
            stereo=stereo if stereo not in {"N"} else None,
            ideal_position=_position(
                (
                    _optional(str(row[7])) or "",
                    _optional(str(row[8])) or "",
                    _optional(str(row[9])) or "",
                )
            ),
        ))

    bond_table = block.find([
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
        "_chem_comp_bond.pdbx_stereo_config",
    ])
    bonds: list[CcdBond] = []
    seen_bonds: set[frozenset[str]] = set()
    for row in bond_table:
        atom1 = _optional(str(row[0])) or ""
        atom2 = _optional(str(row[1])) or ""
        pair = frozenset((atom1, atom2))
        if atom1 == atom2 or not {atom1, atom2} <= names or pair in seen_bonds:
            raise CcdError("CCD bond endpoints must be distinct, known, and unique")
        seen_bonds.add(pair)
        stereo = _optional(str(row[4]))
        bonds.append(CcdBond(
            atom1=atom1,
            atom2=atom2,
            order=(_optional(str(row[2])) or "").upper(),
            aromatic=_flag(_optional(str(row[3])) or ""),
            stereo=stereo if stereo not in {"N"} else None,
        ))

    def optional_float(tag: str) -> float | None:
        value = _optional(block.find_value(tag))
        if value is None:
            return None
        try:
            result = float(value)
            if not math.isfinite(result):
                raise ValueError
            return result
        except ValueError as exc:
            raise CcdError(f"CCD field {tag} must be numeric") from exc

    formal_charge_text = _optional(block.find_value("_chem_comp.pdbx_formal_charge"))
    try:
        formal_charge = int(formal_charge_text) if formal_charge_text is not None else None
    except ValueError as exc:
        raise CcdError("CCD component formal charge must be an integer or unknown") from exc
    known_atom_charges = [atom.formal_charge for atom in atoms if atom.formal_charge is not None]
    if (
        formal_charge is not None
        and len(known_atom_charges) == len(atoms)
        and sum(known_atom_charges) != formal_charge
    ):
        raise CcdError("CCD atom formal charges do not sum to component charge")
    return CcdComponent(
        component_id=component_id,
        name=_required(block, "_chem_comp.name"),
        component_type=_required(block, "_chem_comp.type"),
        formula=_required(block, "_chem_comp.formula"),
        formal_charge=formal_charge,
        formula_weight=optional_float("_chem_comp.formula_weight"),
        initial_date=_optional(block.find_value("_chem_comp.pdbx_initial_date")),
        modified_date=_optional(block.find_value("_chem_comp.pdbx_modified_date")),
        ambiguous=_flag(_required(block, "_chem_comp.pdbx_ambiguous_flag")),
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        provenance=provenance,
    )


def load_ccd_component(
    path: Path,
    component_id: str,
    *,
    expected_sha256: str | None = None,
) -> CcdComponent:
    """Load one explicitly selected local CCD component CIF."""
    path = path.expanduser()
    if path.is_symlink():
        raise CcdError("local CCD input must be a regular non-symlink file")
    path = path.resolve()
    if not path.is_file():
        raise CcdError("local CCD input must be a regular non-symlink file")
    if path.stat().st_size > MAX_CCD_BYTES:
        raise CcdError("local CCD input exceeds the 5 MiB limit")
    content = path.read_bytes()
    if len(content) > MAX_CCD_BYTES:
        raise CcdError("local CCD input exceeds the 5 MiB limit")
    digest = hashlib.sha256(content).hexdigest()
    if expected_sha256 is not None:
        if re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
            raise CcdError("expected CCD SHA-256 must be 64 lowercase hexadecimal characters")
        if digest != expected_sha256:
            raise CcdError("local CCD SHA-256 does not match the required pin")
    return parse_ccd_component(
        content,
        expected_component_id=component_id,
        provenance=CcdProvenance("local-cif", str(path), digest),
    )


def _read_bounded(stream: BinaryIO) -> bytes:
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = stream.read(min(65536, MAX_CCD_BYTES + 1 - total))
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
        if total > MAX_CCD_BYTES:
            raise CcdError("CCD response exceeds the 5 MiB limit")
    return b"".join(chunks)


def fetch_ccd_component(
    component_id: str,
    *,
    cache_dir: Path | None = None,
    timeout_seconds: float = 15.0,
) -> CcdComponent:
    """Fetch one CCD CIF from the fixed HTTPS endpoint and cache by digest."""
    component_id = normalize_component_id(component_id)
    if not 0 < timeout_seconds <= 120:
        raise CcdError("CCD timeout must be in the interval (0, 120]")
    url = f"{CCD_DOWNLOAD_ROOT}/{component_id}.cif"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "dvbfixer-component-info/1"},
        method="GET",
    )
    opener = urllib.request.build_opener(_NoRedirect())
    try:
        with opener.open(request, timeout=timeout_seconds) as response:
            status = getattr(response, "status", 200)
            if status != 200:
                raise CcdError(f"CCD server returned HTTP {status}")
            content_length = response.headers.get("Content-Length")
            if content_length is not None and int(content_length) > MAX_CCD_BYTES:
                raise CcdError("CCD response exceeds the 5 MiB limit")
            content = _read_bounded(response)
    except CcdError:
        raise
    except (OSError, urllib.error.URLError, ValueError) as exc:
        raise CcdError(f"CCD download failed for {component_id}: {exc}") from exc
    digest = hashlib.sha256(content).hexdigest()
    root = (
        cache_dir.expanduser().resolve()
        if cache_dir is not None
        else Path.home() / ".cache" / "dvbfixer" / "ccd"
    )
    destination = root / component_id / f"{digest}.cif"
    component = parse_ccd_component(
        content,
        expected_component_id=component_id,
        provenance=CcdProvenance("online-ccd", str(destination), digest, url),
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.parent.is_symlink() or not destination.parent.is_dir():
        raise CcdError("CCD cache component path must be a non-symlink directory")
    if destination.is_symlink():
        raise CcdError("CCD cache destination must not be a symlink")
    if destination.exists() and destination.read_bytes() != content:
        raise CcdError("CCD digest cache collision")
    if not destination.exists():
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{digest}.", suffix=".tmp", dir=destination.parent
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
    return component


def component_json(component: CcdComponent) -> str:
    return json.dumps(component.to_dict(), ensure_ascii=False, indent=2, sort_keys=True) + "\n"


_BOND_ORDERS = {
    "SING": 1.0,
    "DOUB": 2.0,
    "TRIP": 3.0,
    "AROM": 1.5,
    "DELO": 1.5,
}


def ccd_to_heavy_atom_graph(component: CcdComponent) -> ChemicalGraph:
    """Convert one unambiguous CCD record to a locked heavy-atom geometry graph.

    This does not infer a pH-dependent state. The CCD record's declared formal
    charges, stereochemistry, bond orders, and ideal coordinates are used as one
    exact microspecies. Hydrogen placement remains outside this adapter.
    """
    if component.ambiguous:
        raise CcdError("CCD marks this component as chemically ambiguous")
    heavy_atoms = component.heavy_atoms
    if not heavy_atoms:
        raise CcdError("CCD component contains no heavy atoms")
    if any(atom.ideal_position is None for atom in heavy_atoms):
        raise CcdError("every CCD heavy atom needs ideal coordinates for reconstruction")
    if component.formal_charge is None or any(atom.formal_charge is None for atom in heavy_atoms):
        raise CcdError("CCD formal charges must be explicit for reconstruction")
    names = {atom.name for atom in heavy_atoms}
    atoms = tuple(
        ChemicalAtom(
            atom.name,
            atom.element,
            cast(int, atom.formal_charge),
            atom.ideal_position,
            aromatic=atom.aromatic,
            stereo=atom.stereo,
        )
        for atom in heavy_atoms
        if atom.ideal_position is not None
    )
    positions = {atom.name: atom.ideal_position for atom in atoms}
    bonds: list[ChemicalBond] = []
    for bond in component.bonds:
        if bond.atom1 not in names or bond.atom2 not in names:
            continue
        order = _BOND_ORDERS.get(bond.order)
        if order is None:
            raise CcdError(f"unsupported CCD bond order {bond.order!r}")
        first = positions[bond.atom1]
        second = positions[bond.atom2]
        length = math.dist(first, second)
        if not math.isfinite(length) or length <= 0:
            raise CcdError(f"invalid ideal bond geometry for {bond.atom1}-{bond.atom2}")
        bonds.append(
            ChemicalBond(
                bond.atom1,
                bond.atom2,
                order,
                length,
                aromatic=bond.aromatic or bond.order in {"AROM", "DELO"},
            )
        )
    net_charge = sum(atom.formal_charge for atom in atoms)
    if net_charge != component.formal_charge:
        raise CcdError("CCD heavy-atom charges do not preserve the component formal charge")
    return ChemicalGraph(
        component_id=component.component_id,
        state_id=f"ccd-{component.provenance.sha256[:16]}",
        atoms=atoms,
        bonds=tuple(bonds),
        net_charge=net_charge,
    )
