"""Private terminal-only transport/derivative custody, not a public capability.

The caller must approve the chemical target and check actual native topology and
global L geometry. Native nonactive fidelity is checked BEFORE coordinate-only
materialization; copying observations back is not a substitute for this check.
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, read_records
from dvbfixer.model.diffusion.rfd3_refinement import RefinementCustody, audit_refined_bytes


@dataclass(frozen=True)
class NativeTerminalTransport:
    active_atoms: int
    nonactive_atoms: int
    max_nonactive_residual_angstrom: float
    max_active_displacement_angstrom: float
    scientific_acceptance: bool = False


def audit_native_terminal_transport(
    before: Mapping[AtomIdentity, tuple[float, float, float]],
    after: Mapping[AtomIdentity, tuple[float, float, float]],
    *, active: frozenset[AtomIdentity], max_displacement_angstrom: float,
    fixed_tolerance_angstrom: float = 1e-4,
) -> NativeTerminalTransport:
    """Require complete native input/output axes; no fit, clamp or atom exclusion."""
    for value in (max_displacement_angstrom, fixed_tolerance_angstrom):
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise CustodyError("terminal transport limits must be finite positive numbers")
    if before.keys() != after.keys() or not active or not active < before.keys():
        raise CustodyError("terminal native transport inventory/active mask changed")
    distances = {}
    for atom, xyz in before.items():
        if (len(xyz) != 3 or len(after[atom]) != 3
                or not all(math.isfinite(v) for v in (*xyz, *after[atom]))):
            raise CustodyError("terminal native coordinates are not finite three-vectors")
        distances[atom] = math.dist(xyz, after[atom])
    nonactive_maximum = max(v for atom, v in distances.items() if atom not in active)
    active_maximum = max(distances[atom] for atom in active)
    if nonactive_maximum > fixed_tolerance_angstrom:
        raise CustodyError("native terminal restoration moved a nonactive atom")
    if active_maximum > max_displacement_angstrom:
        raise CustodyError("native terminal restoration exceeded declared movement scope")
    return NativeTerminalTransport(len(active), len(before) - len(active), nonactive_maximum, active_maximum)


def audit_terminal_derivative(
    request: DiffusionRequest, source: bytes, context: bytes, raw: bytes, final: bytes,
    *, active: frozenset[AtomIdentity], raw_sha256: str, context_sha256: str,
    max_displacement_angstrom: float,
) -> RefinementCustody:
    """Protect all inactive proposal atoms as well as experimental observations."""
    if not active or not active <= set(request.generated_atoms):
        raise CustodyError("terminal derivative mask is not a nonempty missing-atom subset")
    report = audit_refined_bytes(request, source, context, raw, final,
        raw_sha256=raw_sha256, context_sha256=context_sha256,
        max_displacement_angstrom=max_displacement_angstrom)
    before, _, _, _ = read_records(raw.decode())
    after, _, _, _ = read_records(final.decode())
    final_records = {atom.identity: atom.record for atom in after}
    if any(atom.record != final_records[atom.identity] for atom in before if atom.identity not in active):
        raise CustodyError("terminal derivative changed an inactive atom record")
    return report
