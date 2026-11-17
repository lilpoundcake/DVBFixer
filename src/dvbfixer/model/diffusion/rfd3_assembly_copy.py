"""Private same-assembly missing-proposal transfer, never experimental truth.

An explicitly approved assembly operation and independently verified observed
frames are caller prerequisites. This pure helper only enforces coordinate and
identity ownership; no donor approval or scientific/public acceptance is granted.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError


def transfer_missing_assembly_copy(
    coordinates: Mapping[AtomIdentity, tuple[float, float, float]],
    *, donor_chain: str, recipient_chain: str,
    donor_atoms: frozenset[AtomIdentity], recipient_atoms: frozenset[AtomIdentity],
    observed: frozenset[AtomIdentity], rotation: Sequence[Sequence[float]],
    translation_angstrom: Sequence[float],
) -> dict[AtomIdentity, tuple[float, float, float]]:
    """Use exact chain/residue/icode/name pairing; only recipient missing atoms move.

    Deposited BIOMT matrices have decimal rounding; permit 1e-5 orthogonality /
    determinant noise, never reflections or general scale/shear transformations.
    The donor MUST also be an unobserved proposal, explicitly identified as such.
    """
    if not donor_chain or not recipient_chain or donor_chain == recipient_chain:
        raise CustodyError("assembly-copy chains must be distinct and case-sensitive")
    if (not donor_atoms or not recipient_atoms or donor_atoms & recipient_atoms
            or (donor_atoms | recipient_atoms) & observed or not observed <= coordinates.keys()
            or not donor_atoms | recipient_atoms <= coordinates.keys()):
        raise CustodyError("assembly-copy donor/recipient ownership is invalid")
    if any(a.chain != donor_chain for a in donor_atoms) or any(a.chain != recipient_chain for a in recipient_atoms):
        raise CustodyError("assembly-copy atom selection has a different exact chain")
    expected = frozenset(AtomIdentity(donor_chain, a.residue_number, a.insertion_code, a.atom_name) for a in recipient_atoms)
    if expected != donor_atoms:
        raise CustodyError("assembly-copy donor/recipient identities do not match exactly")
    matrix, offset = np.asarray(rotation, dtype=float), np.asarray(translation_angstrom, dtype=float)
    if matrix.shape != (3, 3) or offset.shape != (3,) or not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(offset)):
        raise CustodyError("assembly-copy transform is not a finite3x3/3-vector")
    if abs(float(np.linalg.det(matrix)) - 1.) > 1e-5 or float(np.max(np.abs(matrix.T @ matrix - np.eye(3)))) > 1e-5:
        raise CustodyError("assembly-copy transform is not a proper rounded rigid operation")
    if any(len(xyz) != 3 or not all(math.isfinite(v) for v in xyz) for xyz in coordinates.values()):
        raise CustodyError("assembly-copy coordinates are nonfinite or not three-vectors")
    result = dict(coordinates)
    for recipient in recipient_atoms:
        donor = AtomIdentity(donor_chain, recipient.residue_number, recipient.insertion_code, recipient.atom_name)
        result[recipient] = tuple(float(v) for v in matrix @ coordinates[donor] + offset)
    return result
