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
