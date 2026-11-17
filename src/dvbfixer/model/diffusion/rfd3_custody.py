"""Private, offline canonical-proposal adapter. No runner/public registration.

Raw fidelity is checked before generated-only materialization. Callers must
retain raw inputs before invoking this pure service, including on exceptions.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    SourceLedger,
    build_ledger,
    digest,
    expected_names,
    read_records,
)


@dataclass(frozen=True)
class RawProposal:
    atom_axis: tuple[AtomIdentity, ...]
    residue_names: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    fixed_mask: tuple[bool, ...]
    sequence_fixed_mask: tuple[bool, ...]


@dataclass(frozen=True)
class RawPrecision:
    """Explicit caller declaration; no production/default tolerance is inferred."""
    max_observed_residual_angstrom: float

    def __post_init__(self):
        if (isinstance(self.max_observed_residual_angstrom, bool)
                or not math.isfinite(self.max_observed_residual_angstrom) or self.max_observed_residual_angstrom <= 0):
            raise CustodyError("raw precision must be explicit, finite and positive")


@dataclass(frozen=True)
class RawAudit:
    raw_sha256: str
    max_observed_residual_angstrom: float
    per_observation_residuals: tuple[tuple[AtomIdentity, float], ...]
    rotation: tuple[tuple[float, ...], ...]
    translation: tuple[float, ...]
    precision: RawPrecision


@dataclass(frozen=True)
class MaterializedCandidate:
    text: str
    sha256: str
    raw_audit: RawAudit
    scientific_acceptance: bool = False


def audit_raw(request: DiffusionRequest, ledger: SourceLedger, raw: RawProposal, precision: RawPrecision):
    """One proper full-observation fit, no reinjection and no per-chain frames."""
    names = expected_names(request)
    if ledger.canonical_request_sha256 != digest(request.to_json().encode()):
        raise CustodyError("frozen request changed")
    axis = raw.atom_axis
    passthrough = {atom for atom in request.fixed_atoms if atom.atom_name == "OXT"}
    expected = set(names) - passthrough
    if len(axis) != len(set(axis)) or set(axis) != expected:
        raise CustodyError("raw canonical inventory changed; only observed OXT is preserved-only")
    if len(raw.residue_names) != len(axis) or any(name != names[atom] for atom, name in zip(axis, raw.residue_names)):
        raise CustodyError("raw known sequence changed")
    coordinates = np.asarray(raw.coordinates, dtype=np.float64)
    fixed, seq = np.asarray(raw.fixed_mask), np.asarray(raw.sequence_fixed_mask)
    if coordinates.shape != (len(axis), 3) or not np.isfinite(coordinates).all():
        raise CustodyError("raw coordinates must have a finite complete atom axis")
    if fixed.dtype.kind != "b" or fixed.shape != (len(axis),) or seq.dtype.kind != "b" or seq.shape != fixed.shape or not seq.all():
        raise CustodyError("raw masks must be exact boolean axes and all sequence fixed")
    if not np.array_equal(fixed, [atom in request.fixed_atoms for atom in axis]):
        raise CustodyError("raw observation ownership mask changed")
    source = {atom.identity: atom.xyz for atom in ledger.atoms}
    observed = tuple(atom for atom in axis if atom in request.fixed_atoms)
    mobile = coordinates[fixed]
    target = np.array([source[atom] for atom in observed])
    if len(target) < 4 or any(np.linalg.matrix_rank(points - points.mean(axis=0)) < 3 for points in (mobile, target)):
        raise CustodyError("raw/source observed frame is degenerate")
    left, _, right = np.linalg.svd((mobile - mobile.mean(axis=0)).T @ (target - target.mean(axis=0)))
    if np.linalg.det(left @ right) <= 0:
        raise CustodyError("raw frame requires a reflected transform")
    transform = weighted_kabsch(mobile, target)
    aligned = transform.apply(coordinates)
    residual = np.linalg.norm(aligned[fixed] - target, axis=1)
    maximum = float(residual.max())
    if maximum > precision.max_observed_residual_angstrom:
        raise CustodyError(f"raw fixed observations drift before restoration: {maximum:.9g} Angstrom")
    payload = {"axis": [(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in axis],
               "names": raw.residue_names, "coordinates": coordinates.tolist(),
               "fixed_mask": fixed.tolist(), "sequence_fixed_mask": seq.tolist()}
    raw_digest = digest(json.dumps(payload, sort_keys=True, allow_nan=False).encode())
    audit = RawAudit(raw_digest, maximum, tuple(zip(observed, map(float, residual))),
                     tuple(tuple(map(float, row)) for row in transform.rotation), tuple(map(float, transform.translation)), precision)
    return audit, {atom: tuple(map(float, aligned[i])) for i, atom in enumerate(axis)}


def materialize_private_candidate(
    request: DiffusionRequest, ledger: SourceLedger, source: bytes, context: bytes,
    raw: RawProposal, *, precision: RawPrecision,
) -> MaterializedCandidate:
    if build_ledger(request, source, context) != ledger:
        raise CustodyError("source/request/context differs from frozen ledger")
    audit, aligned = audit_raw(request, ledger, raw, precision)
    names = expected_names(request)
    # Deliberately exclude EVERY observation, including preserved-only OXT.
    axis = request.generated_atoms
    text = materialize_candidate_pdb(source.decode("utf-8"), axis,
        tuple(aligned[atom] for atom in axis), tuple(names[atom] for atom in axis), request,
        elements=tuple(atom.atom_name[0] for atom in axis))
    final, link_records, links, auxiliary = read_records(text)
    current = {atom.identity: atom.record for atom in final}
    if (any(current.get(atom.identity) != atom.record for atom in ledger.atoms)
            or link_records != ledger.link_records or links != ledger.links or auxiliary != ledger.auxiliary_records):
        raise CustodyError("materializer changed an observation, auxiliary record or required explicit link")
    return MaterializedCandidate(text, digest(text.encode()), audit)
