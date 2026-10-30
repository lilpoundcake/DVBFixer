"""Validated optional Probe 2.18 adapter; never mutate the source structure."""
from __future__ import annotations

import math
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from dvbfixer.diagnose._common import is_water
from dvbfixer.diagnose.geometry import AtomIdentity
from dvbfixer.diagnose.report import Finding, Severity

SUPPORTED_VERSION = "probe.2.18.211005"
CONTACT_TYPES = frozenset({"wc", "cc", "so", "bo", "hb"})


@dataclass(frozen=True)
class ProbeAtom:
    identity: AtomIdentity
    resname: str
    emitted_resname: str | None = None


@dataclass
class ProbeResult:
    findings: list[Finding] | None = None
    version: str | None = None
    status: str = "failed"
    reason: str | None = None


@dataclass
class ContactPair:
    overlap: float = 0.0
    types: set[str] = field(default_factory=set)


def _serial_key(value: str) -> str:
    value = value.strip()
    try:
        return str(int(value))
    except ValueError:
        return value


def normalize_input(text: str, include_water: bool) -> tuple[str, dict[AtomIdentity, ProbeAtom]]:
    """Align atom fields, select waters, and build an injective identity map.

    Alternate conformations and repeated identities need explicit selection
    before Probe can analyse the same state as the caller's topology.
    """
    from dvbfixer.ffutils.variants import VARIANT_TO_PARENT

    lines = text.splitlines(keepends=True)
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        raise ValueError("select one MODEL before using Probe")
    inventory: dict[AtomIdentity, ProbeAtom] = {}
    aliases: set[AtomIdentity] = set()
    removed_serials: set[str] = set()
    normalized = []
    for line in lines:
        if not include_water and line.startswith("TER") and is_water(line[17:20]):
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            normalized.append(line)
            continue
        resname = line[17:20].strip()
        if not include_water and is_water(resname):
            serial = _serial_key(line[6:11])
            if serial:
                removed_serials.add(serial)
            continue
        if len(line.rstrip("\r\n")) < 78:
            raise ValueError("Probe requires explicit PDB element columns")
        if line[16].strip():
            raise ValueError("select alternate conformations before using Probe")
        if not line[12:27].isascii():
            raise ValueError("Probe atom identifiers must fit ASCII PDB fields")
        name, element = line[12:16].strip(), line[76:78].strip()
        if not name or not resname or not element or not element.isalpha() or len(element) > 2:
            raise ValueError("missing or invalid atom name, residue name or element")
        if not all(math.isfinite(float(line[start:start + 8])) for start in (30, 38, 46)):
            raise ValueError("nonfinite PDB coordinates")
        identity = AtomIdentity(line[21], str(int(line[22:26])), line[26].strip(), name)
        emitted = AtomIdentity(identity.chain, identity.resid, identity.icode, name.upper())
        if identity in inventory or emitted in aliases:
            raise ValueError("Probe atom-name uppercasing makes input identities ambiguous")
        aliases.add(emitted)
        emitted_resname = None
        if line[20].strip():
            # CHARMM ASPP/GLUP occupy the PDB separator column. Probe interprets
            # that fourth letter as part of a two-character chain ID. Normalize
            # only known variants in this private input, preserving all atoms,
            # coordinates, protonation hydrogens and the original display name.
            original_resname = line[17:21].strip()
            emitted_resname = VARIANT_TO_PARENT.get(original_resname)
            if emitted_resname is None or len(emitted_resname) > 3:
                raise ValueError("unsupported residue/chain separator for Probe; unknown four-character residue")
            resname = original_resname
            line = line[:17] + f"{emitted_resname:3s} " + line[21:]
        inventory[identity] = ProbeAtom(identity, resname, emitted_resname)
        atom_field = f" {name:<3}" if len(element) == 1 and len(name) < 4 and not name[0].isdigit() else f"{name:<4}"
        normalized.append(line[:12] + atom_field + line[16:])
    if not inventory:
        raise ValueError("no atoms selected for Probe")
    # Never leave dangling serials after selecting out waters. Drop a complete
    # record if it references an omitted atom, without changing retained serials.
    selected = []
    for line in normalized:
        if line.startswith("CONECT") and removed_serials:
            serials = {_serial_key(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)}
            if serials & removed_serials:
                continue
        selected.append(line)
    return "".join(selected), inventory


def parse_contacts(
    raw: str,
    inventory: dict[AtomIdentity, ProbeAtom],
    version: str,
    radii: str,
    warn_a: float,
    error_a: float,
) -> list[Finding]:
    aliases = {
        AtomIdentity(key.chain, key.resid, key.icode, key.atom.upper()): atom
        for key, atom in inventory.items()
    }
    if len(aliases) != len(inventory):
        raise ValueError("ambiguous Probe identity map")
    pairs: dict[tuple[AtomIdentity, AtomIdentity], ContactPair] = {}
    for line in raw.splitlines():
        if not line.strip():
            continue
        fields = line.split(":")
        if len(fields) != 20 or fields[0] or fields[1] != "1->1" or fields[2] not in CONTACT_TYPES:
            raise ValueError("unsupported Probe condensed output")
        if int(fields[5]) <= 0:
            raise ValueError("invalid Probe dot count")
        numeric = [float(fields[i]) for i in (*range(6, 13), *range(15, 20))]
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError("nonfinite Probe contact data")
        identities = []
        for atom_field in fields[3:5]:
            if len(atom_field) != 16 or atom_field[15].strip():
                raise ValueError("unsupported Probe atom identifier")
            emitted = AtomIdentity(atom_field[:2].strip() or " ", str(int(atom_field[2:6])),
                                   atom_field[6].strip(), atom_field[11:15].strip())
            atom = aliases.get(emitted)
            if atom is None or (atom.emitted_resname or atom.resname).upper() != atom_field[7:10].strip():
                raise ValueError(f"unresolved Probe atom identifier: {atom_field!r}")
            identities.append(atom.identity)
        if identities[0] == identities[1]:
            raise ValueError("Probe emitted a self-contact")
        first, second = sorted(identities)
        pair = pairs.setdefault((first, second), ContactPair())
        pair.types.add(fields[2])
        if fields[2] in {"bo", "so"}:
            pair.overlap = max(pair.overlap, -numeric[0])
    findings = []
    for (first, second), pair in sorted(pairs.items()):
        if pair.overlap < warn_a - 1e-12:
            continue
        left, right = inventory[first], inventory[second]
        left_label = f"{first.chain.strip() or '_'}/{left.resname}{first.residue.display_resid}:{first.atom}"
        right_label = f"{second.chain.strip() or '_'}/{right.resname}{second.residue.display_resid}:{second.atom}"
        findings.append(Finding(
            severity=Severity.ERROR if pair.overlap >= error_a - 1e-12 else Severity.WARNING,
            category="clash", chain=first.chain, resid=first.residue.display_resid,
            resname=left.resname, atom=first.atom,
            message=f"{left_label} clashes with partner {right_label} — overlap {pair.overlap:.2f} Å (via MolProbity probe)",
            fix_hint="dvbfixer minimize",
            extra={"geometry_schema": "dvbfixer.geometry.v1", "metric": "steric_probe_overlap.v1",
                   "engine": {"name": "molprobity-probe", "version": version},
                   "radii_convention": radii, "overlap_angstrom": pair.overlap,
                   "warning_threshold_angstrom": warn_a, "error_threshold_angstrom": error_a,
                   "contact_types": sorted(pair.types), "atom_fields_normalized": True,
                   "atoms": [first.to_dict(), second.to_dict()],
                   "clash_partner": {"chain": second.chain, "resid": second.residue.display_resid,
                                     "icode": second.icode, "resname": right.resname, "atom": second.atom}},
        ))
    return findings


def analyze(
    path: str | Path, binary: str, *, include_water: bool, radii: str, warn_a: float, error_a: float,
) -> ProbeResult:
    if radii not in {"electron", "nuclear"}:
        raise ValueError(f"unknown Probe radii convention: {radii}")
    if not (math.isfinite(warn_a) and math.isfinite(error_a) and 0 < warn_a <= error_a):
        raise ValueError("Probe clash thresholds must be finite, positive and ordered")
    result = ProbeResult()
    try:
        normalized, inventory = normalize_input(Path(path).read_text(), include_water)
    except (OSError, ValueError, UnicodeError) as exc:
        result.status, result.reason = "unsupported_input", str(exc)
        return result
    try:
        version = subprocess.run([binary, "-version"], capture_output=True, text=True, timeout=10, check=False)
        result.version = version.stdout.strip()
        if version.returncode or "ERROR" in version.stderr.upper() or result.version != SUPPORTED_VERSION:
            raise ValueError(f"unsupported or failed Probe version preflight: {result.version!r}; expected {SUPPORTED_VERSION}")
        with tempfile.TemporaryDirectory(prefix="dvbfixer-probe-") as directory:
            temporary = Path(directory) / "input.pdb"
            temporary.write_text(normalized)
            command = [binary, *(["-nuclear"] if radii == "nuclear" else []),
                       "-unformated", "-condensed", "-Q", "-mc", "-het",
                       *(["-waters", "-wat2wat"] if include_water else ["-nowaters"]),
                       "-self", "ALL", str(temporary)]
            native = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
        if native.returncode or "ERROR" in native.stderr.upper():
            raise ValueError(f"Probe failed (exit {native.returncode}): {native.stderr.strip()[:400]}")
        result.findings = parse_contacts(native.stdout, inventory, result.version, radii, warn_a, error_a)
        result.status = "complete"
    except (OSError, ValueError, UnicodeError, subprocess.TimeoutExpired) as exc:
        result.reason = f"{type(exc).__name__}: {exc}"
    return result
