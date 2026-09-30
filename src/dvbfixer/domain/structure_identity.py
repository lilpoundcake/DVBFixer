"""Pure structure-identity concepts; no file-format or toolkit dependencies."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum


class ComponentKind(StrEnum):
    POLYMER = "polymer"
    HETEROGEN = "heterogen"
    SOLVENT = "solvent"
    ION = "ion"


@dataclass(frozen=True, order=True)
class ResidueRef:
    chain_id: str
    sequence_number: int
    insertion_code: str = " "
    name: str = ""


@dataclass(frozen=True, order=True)
class AtomRef:
    serial: int
    residue: ResidueRef
    name: str


@dataclass(frozen=True, order=True)
class ComponentInstanceRef:
    """Exact identity of one structural component occurrence."""

    model: int
    chain_id: str
    sequence_number: str
    insertion_code: str = ""
    alternate_location: str = ""
    occurrence: int = 1

    def __post_init__(self) -> None:
        if self.model < 1:
            raise ValueError("model must be at least 1")
        if not self.chain_id:
            raise ValueError("chain_id must be explicit")
        if not self.sequence_number:
            raise ValueError("sequence_number must be explicit")
        if self.occurrence < 1:
            raise ValueError("occurrence must be at least 1")


@dataclass(frozen=True, order=True)
class ExactAtomRef:
    """Atom identity that preserves MODEL, altloc, and component occurrence."""

    component: ComponentInstanceRef
    name: str

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("atom name must be explicit")


@dataclass(frozen=True)
class MolecularComponent:
    identifier: str
    kind: ComponentKind
    residues: tuple[ResidueRef, ...]
    atom_serials: tuple[int, ...]


def allocate_chain_ids(
    count: int, *, alphabet: Sequence[str], reserved: Iterable[str]
) -> tuple[str, ...]:
    """Allocate unique one-character IDs without reusing reserved identities."""
    unavailable = {value for value in reserved if value and value != " "}
    available = [value for value in alphabet if value not in unavailable]
    if count > len(available):
        raise ValueError(
            f"need {count} unique molecule chain IDs, but only {len(available)} "
            "of the 62 PDB chain IDs remain; use mmCIF output or remove molecules"
        )
    return tuple(available[:count])
