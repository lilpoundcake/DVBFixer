"""Private derivative custody for a retained, already audited RFD3 candidate.

This does not alter native-output fidelity or grant scientific acceptance. Only
missing-atom coordinate columns may differ from the retained raw materialization.
Crystallographic waters enter a disposable FF input, never rewritten context.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    SourceLedger,
    build_ledger,
    digest,
    expected_names,
    read_records,
)
from dvbfixer.model.diffusion.rfd3_water_context import _water_records, assess_water_context


@dataclass(frozen=True)
class RefinementInput:
    text: str
    water_serial_mapping: tuple[tuple[int, int], ...]
    water_atoms: int


@dataclass(frozen=True)
class RefinementCustody:
    raw_candidate_sha256: str
    refined_candidate_sha256: str
    source_sha256: str
    context_sha256: str
    canonical_request_sha256: str
    observations_byte_identical: int
    retained_link_pairs: int
    changed_generated_atoms: int
    generated_rms_movement_angstrom: float
    generated_max_movement_angstrom: float
    movement_limit_angstrom: float
    per_atom_displacements: tuple[tuple[AtomIdentity, float], ...]
    scientific_acceptance: bool = False


def _candidate(request: DiffusionRequest, ledger: SourceLedger, data: bytes):
    atoms, records, links, auxiliary = read_records(data.decode("utf-8"))
    names = expected_names(request)
    current = {atom.identity: atom for atom in atoms}
    if current.keys() != names.keys() or any(
        atom.residue_name != names[atom.identity] or atom.element != atom.identity.atom_name[0] for atom in atoms
    ):
        raise CustodyError("refinement candidate inventory/sequence/element changed")
    if any(current[atom.identity].record != atom.record for atom in ledger.atoms):
        raise CustodyError("refinement changed a fixed source record")
    if records != ledger.link_records or links != ledger.links or auxiliary != ledger.auxiliary_records:
        raise CustodyError("refinement changed source links or auxiliary records")
    return current


def build_refinement_input(
    request: DiffusionRequest, source: bytes, context: bytes, raw: bytes,
    *, context_sha256: str, raw_sha256: str,
) -> RefinementInput:
    """Attach immutable original water oxygens to FF-only input, with unique serials.

    The refiner freezes input heavy atoms outside the exact missing-atom mask.
    Newly introduced temporary FF hydrogens are disposable, not observations.
    """
    if digest(raw) != raw_sha256 or digest(context) != context_sha256:
        raise CustodyError("retained raw/full-context binding changed")
    ledger = build_ledger(request, source, context)
    candidate = _candidate(request, ledger, raw)
    water = assess_water_context(request, source, context, context_sha256=context_sha256, candidate=raw)
    if not water.source_context_clear:
        raise CustodyError("severe immutable source/water defect blocks refinement")
    waters, _ = _water_records(context)
    identities = {AtomIdentity(w.chain, w.residue_number, w.insertion_code, w.atom_name) for w in waters}
    if len(identities) != len(waters) or identities & candidate.keys():
        raise CustodyError("ambiguous protein/water FF atom identity")
    if not raw.endswith(b"END\n"):
        raise CustodyError("private refinement input requires a terminal END newline")
    first_serial = max(atom.serial for atom in candidate.values()) + 1
    if first_serial + len(waters) - 1 > 99999:
        raise CustodyError("temporary FF atom serial capacity exceeded")
    lines, mapping = [], []
    for serial, atom in enumerate(waters, first_serial):
        lines.append(atom.record[:6] + f"{serial:5d}" + atom.record[11:])
        mapping.append((atom.serial, serial))
    text = raw[:-4].decode("utf-8") + "".join(lines) + "END\n"
    # Canonical strict reader checks all temporary serials/identities; no link or
    # spatial proximity is used to reinterpret protein-water pairs as covalent.
    atoms, _, _, _ = read_records(text)
    if len(atoms) != len(candidate) + len(waters):
        raise CustodyError("temporary FF context lost atoms")
    return RefinementInput(text, tuple(mapping), len(waters))


def audit_refined_bytes(
    request: DiffusionRequest, source: bytes, context: bytes, raw: bytes, refined: bytes,
    *, raw_sha256: str, context_sha256: str, max_displacement_angstrom: float,
) -> RefinementCustody:
    """Separate derivative audit, never pretending refined atoms equal native raw."""
    if (isinstance(max_displacement_angstrom, bool) or not math.isfinite(max_displacement_angstrom)
            or max_displacement_angstrom <= 0):
        raise CustodyError("refinement movement limit must be explicit finite positive")
    if digest(raw) != raw_sha256 or digest(context) != context_sha256:
        raise CustodyError("retained raw/full-context binding changed")
    ledger = build_ledger(request, source, context)
    before, after = _candidate(request, ledger, raw), _candidate(request, ledger, refined)
    initial, final = raw.splitlines(keepends=True), refined.splitlines(keepends=True)
    if len(initial) != len(final):
        raise CustodyError("refinement changed physical record layout")
    for old, new in zip(initial, final, strict=True):
        identity = None
        if old.startswith((b"ATOM  ", b"HETATM")):
            identity = AtomIdentity(old[21:22].decode(), old[22:26].strip().decode(),
                                    old[26:27].strip().decode(), old[12:16].strip().decode())
        if identity in request.generated_atoms:
            if old[:30] != new[:30] or old[54:] != new[54:]:
                raise CustodyError("refinement changed missing-atom metadata")
        elif old != new:
            raise CustodyError("refinement changed a non-mobile record")
    distances = tuple((atom, math.dist(before[atom].xyz, after[atom].xyz)) for atom in request.generated_atoms)
    maximum = max(distance for _, distance in distances)
    if maximum > max_displacement_angstrom:
        raise CustodyError(f"refinement exceeded declared missing-atom movement: {maximum:.6g} Angstrom")
    rms = math.sqrt(sum(distance * distance for _, distance in distances) / len(distances))
    return RefinementCustody(digest(raw), digest(refined), ledger.source_sha256, ledger.context_sha256,
                             ledger.canonical_request_sha256, len(ledger.atoms), len(ledger.links),
                             sum(distance > 0 for _, distance in distances), rms, maximum,
                             max_displacement_angstrom, distances)
