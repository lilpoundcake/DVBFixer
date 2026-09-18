"""Compatibility adapter for force-field naming conversion."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from dvbfixer.domain.force_field_naming import (
    AMBER_VARIANTS,
    CHARMM_UNIVERSAL_ATOM_RENAMES,
    CHARMM_VARIANTS,
    GROMACS_AMBER_ATOM_RENAMES,
    GROMACS_AMBER_LYN_ATOM_RENAME,
    GROMACS_AMBER_NA_ATOM_RENAMES,
    GROMACS_CHARMM_CAP_ATOM_RENAMES,
    PROTONATION_AMBER_TO_CHARMM,
    PROTONATION_ATOM_RENAME_TO_AMBER,
    PROTONATION_ATOM_RENAME_TO_CHARMM,
    PROTONATION_CHARMM_TO_AMBER,
    ForceFieldTarget,
    NamingProfile,
    ResidueIdentity,
    VariantOverride,
)
from dvbfixer.force_field_naming import NamingConversionRequest, convert_force_field_naming

__all__ = [
    "AMBER_VARIANTS",
    "CHARMM_UNIVERSAL_ATOM_RENAMES",
    "CHARMM_VARIANTS",
    "GROMACS_AMBER_ATOM_RENAMES",
    "GROMACS_AMBER_LYN_ATOM_RENAME",
    "GROMACS_AMBER_NA_ATOM_RENAMES",
    "GROMACS_CHARMM_CAP_ATOM_RENAMES",
    "PROTONATION_AMBER_TO_CHARMM",
    "PROTONATION_ATOM_RENAME_TO_AMBER",
    "PROTONATION_ATOM_RENAME_TO_CHARMM",
    "PROTONATION_CHARMM_TO_AMBER",
    "apply_variants_to_pdb_text",
]


def _split_key(key: object) -> ResidueIdentity:
    if isinstance(key, tuple) and len(key) == 2:
        return ResidueIdentity(str(key[0]), str(key[1]), "")
    if isinstance(key, tuple) and len(key) == 3:
        return ResidueIdentity(str(key[0]), str(key[1]), str(key[2] or ""))
    raise TypeError(f"Unexpected key shape for amber_renames: {key!r}")


def _atomic_replace(path: Path, data: bytes) -> None:
    destination = path.resolve(strict=True) if path.is_symlink() else path
    mode = destination.stat().st_mode
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", dir=destination.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary_name, mode)
        os.replace(temporary_name, destination)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def apply_variants_to_pdb_text(
    pdb_path: str | Path,
    amber_renames: dict,
    target_ff: str = "amber",
    include_gromacs_shifts: bool = True,
    verbose: bool = False,
) -> int:
    """Rewrite a PDB in place and return changed ATOM/HETATM line count."""
    if target_ff not in ("amber", "charmm"):
        raise ValueError(f"target_ff must be 'amber' or 'charmm', got {target_ff!r}")
    overrides = tuple(
        VariantOverride(_split_key(key), value)
        for key, value in (amber_renames or {}).items()
    )
    path = Path(pdb_path)
    source = path.read_bytes()
    request = NamingConversionRequest(
        pdb_text=source.decode("latin-1"),
        target=ForceFieldTarget(target_ff),
        profile=NamingProfile.GROMACS if include_gromacs_shifts else NamingProfile.STANDARD,
        variant_overrides=overrides,
    )
    result = convert_force_field_naming(request)
    if verbose:
        for change in result.changes:
            print(
                f"  [ff_names] {change.chain}/{change.source_residue_name}"
                f"{change.residue_number}:{change.source_atom_name} -> "
                f"{change.target_residue_name}:{change.target_atom_name}"
            )
    converted = result.pdb_text.encode("latin-1")
    if converted != source:
        _atomic_replace(path, converted)
    return result.summary.changed_atoms
