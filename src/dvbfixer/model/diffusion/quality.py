"""Diffusion-only conformational quality gates."""

from __future__ import annotations

import math
from importlib.metadata import distribution
from typing import Any

import numpy as np

from dvbfixer.diagnose.chemistry import (
    _bonded_atom_pairs,
    _dihedral_deg,
    _find_atom,
    _pos,
    _res_loc,
)
from dvbfixer.diagnose.report import Finding, Severity

_RAMACHANDRAN_BIN_DEGREES = 4.0
_RAMACHANDRAN_NEIGHBORHOOD_BINS = 2
_RAMACHANDRAN_REFERENCE_PATH = distribution("MDAnalysis").locate_file(
    "MDAnalysis/analysis/data/rama_ref_data.npy"
)
_RAMACHANDRAN_REFERENCE: np.ndarray | None = None

_SIDECHAIN_CHI12_BIN_DEGREES = 6.0
_SIDECHAIN_CHI12_NEIGHBORHOOD_BINS = 2
_SIDECHAIN_CHI12_REFERENCE_PATH = distribution("MDAnalysis").locate_file(
    "MDAnalysis/analysis/data/janin_ref_data.npy"
)
_SIDECHAIN_CHI12_REFERENCE: np.ndarray | None = None
_SIDECHAIN_CHI12_ATOMS: dict[str, tuple[str, str, str, str, str]] = {
    "ARG": ("N", "CA", "CB", "CG", "CD"),
    "ASN": ("N", "CA", "CB", "CG", "OD1"),
    "ASP": ("N", "CA", "CB", "CG", "OD1"),
    "GLN": ("N", "CA", "CB", "CG", "CD"),
    "GLU": ("N", "CA", "CB", "CG", "CD"),
    "HIS": ("N", "CA", "CB", "CG", "ND1"),
    "ILE": ("N", "CA", "CB", "CG1", "CD1"),
    "LEU": ("N", "CA", "CB", "CG", "CD1"),
    "LYS": ("N", "CA", "CB", "CG", "CD"),
    "MET": ("N", "CA", "CB", "CG", "SD"),
    "PHE": ("N", "CA", "CB", "CG", "CD1"),
    "TRP": ("N", "CA", "CB", "CG", "CD1"),
    "TYR": ("N", "CA", "CB", "CG", "CD1"),
}
_SIDECHAIN_CHI2_SYMMETRIC = frozenset({"ASP", "LEU", "PHE", "TYR"})
_CANONICAL_AMINO_ACIDS = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
)


def _ramachandran_reference() -> np.ndarray:
    global _RAMACHANDRAN_REFERENCE
    if _RAMACHANDRAN_REFERENCE is None:
        _RAMACHANDRAN_REFERENCE = np.load(str(_RAMACHANDRAN_REFERENCE_PATH))
    return _RAMACHANDRAN_REFERENCE


def _ramachandran_reference_populated(phi: float, psi: float) -> bool:
    reference = _ramachandran_reference()
    phi_bin = int(math.floor((phi + 180.0) / _RAMACHANDRAN_BIN_DEGREES)) % 90
    psi_bin = int(math.floor((psi + 180.0) / _RAMACHANDRAN_BIN_DEGREES)) % 90
    radius = _RAMACHANDRAN_NEIGHBORHOOD_BINS
    return any(
        reference[(psi_bin + psi_offset) % 90, (phi_bin + phi_offset) % 90] >= 1.0
        for psi_offset in range(-radius, radius + 1)
        for phi_offset in range(-radius, radius + 1)
    )


def check_ramachandran(topology: Any, positions: Any) -> list[Finding]:
    """Report gross general-residue Ramachandran outliers."""
    findings: list[Finding] = []
    bonds = _bonded_atom_pairs(topology)
    for chain in topology.chains():
        residues = list(chain.residues())
        for previous, residue, following in zip(residues, residues[1:], residues[2:]):
            if (
                residue.name not in _CANONICAL_AMINO_ACIDS
                or residue.name in {"GLY", "PRO"}
                or following.name == "PRO"
            ):
                continue
            previous_c = _find_atom(previous, "C")
            n = _find_atom(residue, "N")
            ca = _find_atom(residue, "CA")
            c = _find_atom(residue, "C")
            following_n = _find_atom(following, "N")
            atoms = (previous_c, n, ca, c, following_n)
            if any(atom is None for atom in atoms):
                continue
            assert previous_c is not None
            assert n is not None
            assert ca is not None
            assert c is not None
            assert following_n is not None
            atoms = (previous_c, n, ca, c, following_n)
            required_bonds = ((previous_c, n), (n, ca), (ca, c), (c, following_n))
            if any(
                tuple(sorted((first.index, second.index))) not in bonds
                for first, second in required_bonds
            ):
                continue
            phi = _dihedral_deg(*(_pos(positions, atom) for atom in atoms[:4]))
            psi = _dihedral_deg(*(_pos(positions, atom) for atom in atoms[1:]))
            if phi is None or psi is None or _ramachandran_reference_populated(phi, psi):
                continue
            chain_id, resid = _res_loc(residue)
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    category="ramachandran_outlier",
                    chain=chain_id,
                    resid=resid,
                    resname=residue.name,
                    atom="N-CA-C",
                    message=(
                        f"gross Ramachandran outlier (phi={phi:.1f} deg, "
                        f"psi={psi:.1f} deg)"
                    ),
                    fix_hint="manual: rebuild or refine the local backbone geometry",
                    extra={"phi_degrees": phi, "psi_degrees": psi},
                )
            )
    return findings


def _sidechain_chi12_reference() -> np.ndarray:
    global _SIDECHAIN_CHI12_REFERENCE
    if _SIDECHAIN_CHI12_REFERENCE is None:
        _SIDECHAIN_CHI12_REFERENCE = np.load(str(_SIDECHAIN_CHI12_REFERENCE_PATH))
    return _SIDECHAIN_CHI12_REFERENCE


def _sidechain_chi12_reference_populated(chi1: float, chi2: float) -> bool:
    reference = _sidechain_chi12_reference()
    chi1_bin = int(math.floor((chi1 % 360.0) / _SIDECHAIN_CHI12_BIN_DEGREES)) % 60
    chi2_bin = int(math.floor((chi2 % 360.0) / _SIDECHAIN_CHI12_BIN_DEGREES)) % 60
    radius = _SIDECHAIN_CHI12_NEIGHBORHOOD_BINS
    return any(
        reference[(chi2_bin + dy) % 60, (chi1_bin + dx) % 60] >= 1.0
        for dy in range(-radius, radius + 1)
        for dx in range(-radius, radius + 1)
    )


def check_sidechain_chi12(topology: Any, positions: Any) -> list[Finding]:
    """Report gross pooled chi1/chi2 side-chain torsion outliers."""
    findings: list[Finding] = []
    bonds = _bonded_atom_pairs(topology)
    for residue in topology.residues():
        names = _SIDECHAIN_CHI12_ATOMS.get(residue.name)
        if names is None:
            continue
        atoms_by_name: dict[str, list[Any]] = {}
        for atom in residue.atoms():
            atoms_by_name.setdefault(atom.name, []).append(atom)
        if any(len(atoms_by_name.get(name, ())) != 1 for name in names):
            continue
        atoms = tuple(atoms_by_name[name][0] for name in names)
        if any(
            tuple(sorted((first.index, second.index))) not in bonds
            for first, second in zip(atoms, atoms[1:])
        ):
            continue
        chi1 = _dihedral_deg(*(_pos(positions, atom) for atom in atoms[:4]))
        chi2 = _dihedral_deg(*(_pos(positions, atom) for atom in atoms[1:]))
        if chi1 is None or chi2 is None:
            continue
        populated = _sidechain_chi12_reference_populated(chi1, chi2)
        if residue.name in _SIDECHAIN_CHI2_SYMMETRIC:
            populated = populated or _sidechain_chi12_reference_populated(chi1, chi2 + 180.0)
        if populated:
            continue
        chain_id, resid = _res_loc(residue)
        findings.append(
            Finding(
                severity=Severity.ERROR,
                category="sidechain_chi12_outlier",
                chain=chain_id,
                resid=resid,
                resname=residue.name,
                atom="-".join(names),
                message=f"gross side-chain chi1/chi2 outlier (chi1={chi1:.1f}, chi2={chi2:.1f})",
                fix_hint="manual: rebuild or refine the side-chain conformation",
                extra={"chi1_degrees": chi1, "chi2_degrees": chi2},
            )
        )
    return findings
