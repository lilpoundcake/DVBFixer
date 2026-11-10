"""Research-only observed-anchor displacement measurements, not policy changes.

Fits exclude the measured residue. Experimental conformer differences combine
uncertainty, physical motion and model bias; they are not direct error radii.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.geometry import weighted_kabsch

CALIBRATION_REVISION = "leave-anchor-out-experimental-comparison-v1"
PROPOSED_LIMITS = (
    ("current", AnchorMovementLimits()),
    ("candidate-075", AnchorMovementLimits(0.75, 1.5, 0.75)),
    ("candidate-100", AnchorMovementLimits(1.0, 2.0, 1.0)),
    ("candidate-150", AnchorMovementLimits(1.5, 2.5, 1.5)),
)


@dataclass(frozen=True)
class AnchorComparison:
    residue: ResidueIdentity
    matched_heavy_atoms: int
    fit_atoms: int
    fit_rmsd_angstrom: float
    backbone_rmsd_angstrom: float
    ca_displacement_angstrom: float
    maximum_heavy_displacement_angstrom: float
    contact_pairs_reference: int
    contact_pairs_lost: int
    contact_pairs_gained: int
    limits_passed: tuple[str, ...]
    revision: str = CALIBRATION_REVISION


def compare_anchor(
    reference: Mapping[AtomIdentity, tuple[float, float, float]],
    mobile: Mapping[AtomIdentity, tuple[float, float, float]],
    residue: ResidueIdentity,
    fit_atoms: tuple[AtomIdentity, ...],
    *,
    contact_atoms: tuple[AtomIdentity, ...] = (),
    contact_cutoff_angstrom: float = 4.0,
) -> AnchorComparison:
    """Exact-name matched, leave-anchor-out local fit with explicit contact axis.

    Caller must establish sequence/chemical equivalence and heavy-only axes;
    unsupported mappings/altlocs must not be guessed by coordinate proximity.
    """
    def belongs(atom):
        return (atom.chain, atom.residue_number, atom.insertion_code) == (
            residue.chain, residue.residue_number, residue.insertion_code)

    if (len(fit_atoms) < 6 or len(set(fit_atoms)) != len(fit_atoms)
            or any(belongs(a) for a in fit_atoms)):
        raise ValueError("calibration fit requires six unique atoms outside the anchor")
    if (len(set(contact_atoms)) != len(contact_atoms) or any(belongs(a) for a in contact_atoms)
            or not np.isfinite(contact_cutoff_angstrom) or contact_cutoff_angstrom <= 0):
        raise ValueError("calibration contact axis/cutoff is invalid")
    anchors = tuple(sorted(a for a in reference if belongs(a)))
    if set(anchors) != {a for a in mobile if belongs(a)} or not anchors:
        raise ValueError("calibration requires an identical complete heavy-atom anchor axis")
    if not {"N", "CA", "C", "O"} <= {a.atom_name for a in anchors}:
        raise ValueError("calibration anchor backbone is incomplete")

    def points(mapping, atoms):
        try:
            result = np.asarray([mapping[a] for a in atoms], dtype=float)
        except KeyError as exc:
            raise ValueError("calibration atom axis is not shared") from exc
        if result.shape != (len(atoms), 3) or not np.isfinite(result).all():
            raise ValueError("calibration coordinates are nonfinite or malformed")
        return result

    fixed_fit, moving_fit = points(reference, fit_atoms), points(mobile, fit_atoms)
    if min(np.linalg.matrix_rank(p - p.mean(axis=0), tol=1e-8) for p in (fixed_fit, moving_fit)) < 2:
        raise ValueError("calibration fitting frame is degenerate")
    transform = weighted_kabsch(moving_fit, fixed_fit)
    fixed_anchor = points(reference, anchors)
    moved_anchor = transform.apply(points(mobile, anchors))
    displacements = np.linalg.norm(moved_anchor - fixed_anchor, axis=1)
    backbone = np.array([a.atom_name in {"N", "CA", "C", "O"} for a in anchors])
    rmsd = float(np.sqrt(np.mean(displacements[backbone] ** 2)))
    ca = float(displacements[next(i for i, a in enumerate(anchors) if a.atom_name == "CA")])
    maximum = float(displacements.max())
    if contact_atoms:
        ref_contacts = np.linalg.norm(fixed_anchor[:, None] - points(reference, contact_atoms)[None, :], axis=2) <= contact_cutoff_angstrom
        mobile_contacts = np.linalg.norm(moved_anchor[:, None] - transform.apply(points(mobile, contact_atoms))[None, :], axis=2) <= contact_cutoff_angstrom
        contact_counts = (int(ref_contacts.sum()), int((ref_contacts & ~mobile_contacts).sum()),
                          int((mobile_contacts & ~ref_contacts).sum()))
    else:
        contact_counts = (0, 0, 0)
    return AnchorComparison(residue, len(anchors), len(fit_atoms),
        float(np.sqrt(np.mean(np.sum((transform.apply(moving_fit) - fixed_fit) ** 2, axis=1)))),
        rmsd, ca, maximum, *contact_counts,
        tuple(name for name, limits in PROPOSED_LIMITS
              if rmsd <= limits.backbone_rmsd_angstrom_max and ca <= limits.ca_displacement_angstrom_max
              and maximum <= limits.heavy_displacement_angstrom_max))
