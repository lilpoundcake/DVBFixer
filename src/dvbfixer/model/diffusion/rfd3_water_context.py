"""Private immutable water-obstacle assessment, not native water conditioning.

Protein/water cross pairs are noncovalent by this narrow contract. Short contacts
must not disappear through generic proximity-based covalent-bond exclusions.
Canonical radii, clash presets and heavy H-bond classification are reused.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from dvbfixer.diagnose._common import is_water
from dvbfixer.diagnose.geometry import VDW_RADII_ANGSTROM
from dvbfixer.diagnose.report import Severity
from dvbfixer.diagnose.steric import CLASH_MODE_PRESETS, DEFAULT_CLASH_MODE, _is_heavy_hbond_pair
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.rfd3_custody_records import (
    CustodyError,
    SourceAtom,
    digest,
    expected_names,
    read_records,
)

WATER_CONTEXT_REVISION = "private-fixed-water-obstacles-v1"


@dataclass(frozen=True)
class WaterAtom:
    chain: str
    residue_number: str
    insertion_code: str
    residue_name: str
    atom_name: str
    altloc: str
    serial: int
    occupancy: float
    xyz: tuple[float, float, float]
    record: str


@dataclass(frozen=True)
class WaterContact:
    protein_atom: AtomIdentity
    water: WaterAtom
    distance_angstrom: float
    overlap_angstrom: float
    severity: Severity


@dataclass(frozen=True)
class WaterAssessment:
    revision: str
    context_sha256: str
    water_atoms: int
    nonunit_water_occupancies: int
    water_altlocs: int
    baseline_contacts: tuple[WaterContact, ...]
    restored_contacts: tuple[WaterContact, ...] | None
    native_water_conditioning: bool = False
    scientific_acceptance: bool = False

    @property
    def source_context_clear(self) -> bool:
        return not any(contact.severity == Severity.ERROR for contact in self.baseline_contacts)

    @property
    def restored_region_clear(self) -> bool | None:
        if self.restored_contacts is None:
            return None
        return not any(contact.severity == Severity.ERROR for contact in self.restored_contacts)


def _water_records(context: bytes):
    waters, protein, serials = [], [], set()
    water_identities = set()
    for line in context.decode("utf-8").splitlines(keepends=True):
        kind = line[:6].strip()
        if kind in {"MODEL", "ENDMDL"}:
            raise CustodyError("multi-model context is not admitted")
        if kind in {"LINK", "CONECT"}:
            raise CustodyError("LINK/CONECT-bearing full context requires a separate component-link audit")
        if kind not in {"ATOM", "HETATM"}:
            continue
        serial = int(line[6:11])
        if serial in serials:
            raise CustodyError("duplicate full-context atom serial")
        serials.add(serial)
        name = line[17:20].strip()
        if not is_water(name):
            protein.append(line)
            continue
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        occupancy = float(line[54:60])
        element = line[76:78].strip() or line[12:16].strip()[:1]
        if element != "O" or line[12:16].strip() not in {"O", "OW"}:
            raise CustodyError("this context contract requires oxygen-only crystallographic waters")
        if not all(math.isfinite(v) for v in (*xyz, occupancy)) or not 0 <= occupancy <= 1:
            raise CustodyError("invalid water coordinates or occupancy")
        identity = (line[21], line[22:26].strip(), line[26].strip(), name, line[12:16].strip(), line[16].strip())
        if identity in water_identities:
            raise CustodyError("ambiguous duplicate water identity")
        water_identities.add(identity)
        waters.append(WaterAtom(line[21], line[22:26].strip(), line[26].strip(), name,
                                line[12:16].strip(), line[16].strip(), serial, occupancy, xyz, line))
    if not waters:
        raise CustodyError("declared water context is empty")
    return tuple(waters), "".join(protein)


def _contacts(protein: tuple[SourceAtom, ...], waters: tuple[WaterAtom, ...]) -> tuple[WaterContact, ...]:
    if not protein:
        return ()
    warning, error = CLASH_MODE_PRESETS[DEFAULT_CLASH_MODE]
    water_points = np.array([water.xyz for water in waters])
    tree = cKDTree(water_points)
    # 3 A contains every warning for the admitted C/N/O/S versus water-O pairs.
    found = []
    oxygen = SimpleNamespace(element=SimpleNamespace(symbol="O"))
    for atom in protein:
        if atom.element not in VDW_RADII_ANGSTROM:
            raise CustodyError("unknown protein element for water assessment")
        radius = VDW_RADII_ANGSTROM[atom.element] + VDW_RADII_ANGSTROM["O"]
        if radius - warning > 3.0:
            raise CustodyError("water neighbor cutoff no longer covers canonical warning range")
        probe = SimpleNamespace(element=SimpleNamespace(symbol=atom.element))
        for index in sorted(tree.query_ball_point(atom.xyz, r=3.0)):
            water = waters[index]
            distance = float(np.linalg.norm(np.asarray(atom.xyz) - water_points[index]))
            overlap = radius - distance
            if overlap < warning - 1e-12 or _is_heavy_hbond_pair(probe, oxygen, distance):
                continue
            # Never infer a covalent protein-water bond or exclude a short pair.
            found.append(WaterContact(atom.identity, water, distance, overlap,
                                      Severity.ERROR if overlap >= error - 1e-12 else Severity.WARNING))
    return tuple(found)


def assess_water_context(
    request: DiffusionRequest, source: bytes, context: bytes, *, context_sha256: str,
    candidate: bytes | None = None,
) -> WaterAssessment:
    if digest(source) != request.normalized_pdb.sha256 or digest(context) != context_sha256:
        raise CustodyError("source/full-context binding changed")
    source_atoms, _, _, _ = read_records(source.decode("utf-8"))
    source_map = {atom.identity: atom for atom in source_atoms}
    names = expected_names(request)
    if source_map.keys() != set(request.fixed_atoms):
        raise CustodyError("source fixed inventory differs from frozen request")
    if any(atom.residue_name != names[atom.identity] or atom.element != atom.identity.atom_name[0] for atom in source_atoms):
        raise CustodyError("source sequence/element changed during water audit")
    waters, protein_text = _water_records(context)
    original_atoms, _, _, _ = read_records(protein_text)
    originals = {atom.identity: atom for atom in original_atoms}
    if originals.keys() != source_map.keys() or any(
        originals[key].xyz != atom.xyz or originals[key].residue_name != atom.residue_name
        or originals[key].element != atom.element for key, atom in source_map.items()
    ):
        raise CustodyError("omitted context is not waters-only in the unchanged source frame")
    restored = None
    if candidate is not None:
        final_atoms, _, _, _ = read_records(candidate.decode("utf-8"))
        final = {atom.identity: atom for atom in final_atoms}
        if final.keys() != names.keys() or any(atom.residue_name != names[atom.identity] or atom.element != atom.identity.atom_name[0] for atom in final_atoms):
            raise CustodyError("candidate inventory/sequence changed during water audit")
        if any(final[key].record != atom.record for key, atom in source_map.items()):
            raise CustodyError("candidate changed fixed source records during water audit")
        restored = _contacts(tuple(atom for atom in final_atoms if atom.identity in request.generated_atoms), waters)
    return WaterAssessment(WATER_CONTEXT_REVISION, digest(context), len(waters),
                           sum(water.occupancy != 1 for water in waters), sum(bool(water.altloc) for water in waters),
                           _contacts(source_atoms, waters), restored)
