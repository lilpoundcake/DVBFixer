"""Private single-target PDB adapter; no production routing or filesystem writes.

Canonical protein/water only. Source observations are parsed in full, never
accepted from a caller. Ambiguous layouts and linked chemistry fail closed.
"""
from __future__ import annotations

import hashlib
import math
from collections.abc import Callable
from dataclasses import dataclass, replace

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.dat import AddedAtom, DatRecord
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS
from dvbfixer.pdb_component_reconstruction import (
    PdbAtom,
    _component_ref,
    _format_generated_atom,
    _parse_atom_lines,
)
from dvbfixer.protein_completion.contracts import (
    Addition,
    Candidate,
    CompletionRequest,
    Observation,
)
from dvbfixer.protein_completion.geometry import cross, dot, sub
from dvbfixer.protein_completion.links import parse_declared_links
from dvbfixer.protein_completion.rotamers import SingleChiLibrary
from dvbfixer.protein_completion.validation import (
    admission_reason,
    audit_candidate,
    sidechain_bonds,
)

# The separate branched_source entry explicitly opts into VAL/THR. Never enroll
# these default source/cohort/publication consumers when the kernel expands.
SOURCE_RESIDUES = frozenset({"ALA", "SER", "CYS"})


class SourceCompletionError(ValueError):
    """Source binding, supported scope or materialized geometry failed."""


@dataclass(frozen=True)
class MaterializedCompletion:
    pdb: bytes
    source_sha256: str
    output_sha256: str
    request_digest: str
    generated: tuple[Addition, ...]  # read back from rounded PDB, not raw proposal

    def dat_record(self) -> DatRecord:
        """Fresh canonical sidecar; no observed atom can acquire generated ownership."""
        name = self._residue_name()
        atoms: list[AddedAtom] = [
            {"chain": a.atom.component.chain_id, "resid": a.atom.component.sequence_number,
             "icode": a.atom.component.insertion_code, "resname": name,
             "atom": a.atom.name, "element": a.element} for a in self.generated
        ]
        return DatRecord(description="Private canonical sidechain completion; generated atoms only",
                         added_atoms=atoms)

    def _residue_name(self) -> str:
        # Derive from independently parsed output, not a second caller authority.
        atoms = _parse_atom_lines(self.pdb.decode("ascii").splitlines(keepends=True))
        target = self.generated[0].atom.component if self.generated else None
        return next((a.residue_name for a in atoms if _component_ref(a) == target), "")


def _layout(source: bytes) -> tuple[list[str], list[PdbAtom], set[int]]:
    if not isinstance(source, bytes) or not 0 < len(source) <= 16*1024*1024:
        raise SourceCompletionError("bounded nonempty PDB bytes required")
    try:
        text = source.decode("ascii")
    except UnicodeDecodeError as exc:
        raise SourceCompletionError("ASCII PDB required") from exc
    lines = text.splitlines(keepends=True)
    if any(not line.endswith("\n") or "\r" in line for line in lines):
        raise SourceCompletionError("LF-terminated source records required")
    ended = False
    for line in lines:
        record = line[:6].strip()
        if ((record in {"ATOM", "HETATM", "ANISOU", "CONECT"}
             and not line.startswith({"ATOM": "ATOM  ", "HETATM": "HETATM", "ANISOU": "ANISOU", "CONECT": "CONECT"}[record]))
                or (line.startswith("TER") and record != "TER")):
            raise SourceCompletionError("malformed structural record prefix")
        if any(line.startswith(prefix) and record != prefix for prefix in ("LINK", "SSBOND")):
            raise SourceCompletionError("malformed structural record prefix")
        if ended and line.strip():
            raise SourceCompletionError("records after END")
        if record in {"MODEL", "ENDMDL", "MASTER"}:
            raise SourceCompletionError(f"unsupported source record: {record}")
        if record == "END":
            ended = True
    atoms = _parse_atom_lines(lines)
    if not atoms:
        raise SourceCompletionError("no source atoms")
    serials = {a.serial for a in atoms}
    by_serial = {a.serial: a for a in atoms}
    seen_anisou: set[int] = set()
    for line in lines:
        record = line[:6].strip()
        if record == "TER" and line[6:11].strip():
            try:
                serial = int(line[6:11])
            except ValueError as exc:
                raise SourceCompletionError("invalid TER serial") from exc
            if not 1 <= serial <= 99999 or serial in serials:
                raise SourceCompletionError("duplicate/out-of-range TER serial")
            serials.add(serial)
        if record == "ANISOU":
            try:
                serial = int(line[6:11])
                values = [int(line[i:i+7]) for i in range(28, 70, 7)]
            except ValueError as exc:
                raise SourceCompletionError("invalid ANISOU") from exc
            atom = by_serial.get(serial)
            if (atom is None or serial in seen_anisou or len(values) != 6
                    or line[12:27] != lines[atom.line_index][12:27]
                    or line[76:78] != lines[atom.line_index][76:78]):
                raise SourceCompletionError("orphan/duplicate/mismatched ANISOU")
            seen_anisou.add(serial)
    seen: set[ExactAtomRef] = set()
    names: dict[ComponentInstanceRef, str] = {}
    for atom in atoms:
        ref = _component_ref(atom)
        identity = ExactAtomRef(ref, atom.name)
        if atom.occurrence != 1 or identity in seen or names.get(ref, atom.residue_name) != atom.residue_name:
            raise SourceCompletionError("ambiguous/split residue or atom identity")
        seen.add(identity)
        names[ref] = atom.residue_name
        try:
            occupancy = float(lines[atom.line_index][54:60])
            bfactor = float(lines[atom.line_index][60:66])
            number = int(atom.sequence_number)
        except ValueError as exc:
            raise SourceCompletionError("invalid occupancy/B-factor/residue number") from exc
        if (atom.alternate_location or occupancy != 1 or not math.isfinite(bfactor)
                or str(number) != atom.sequence_number or not -999 <= number <= 9999):
            raise SourceCompletionError("altloc/partial occupancy or noncanonical residue number")
        if atom.residue_name not in CANONICAL_HEAVY_ATOMS and atom.residue_name != "HOH":
            raise SourceCompletionError("source scope is canonical protein and water only")
        allowed = (*CANONICAL_HEAVY_ATOMS[atom.residue_name], "OXT") if atom.residue_name != "HOH" else ("O",)
        if atom.name not in allowed or atom.element != atom.name[0]:
            raise SourceCompletionError("source scope requires canonical heavy-atom names/elements")
    return lines, atoms, serials


def _global_ca(atoms: list[PdbAtom], *, allow_missing_cb: ComponentInstanceRef | None,
               additional_missing_cb: frozenset[ComponentInstanceRef] = frozenset()) -> None:
    residues: dict[ComponentInstanceRef, list[PdbAtom]] = {}
    for atom in atoms:
        residues.setdefault(_component_ref(atom), []).append(atom)
    for ref, group in residues.items():
        if group[0].residue_name in {"GLY", "HOH"}:
            continue
        xyz = {a.name: a.position for a in group}
        needed = {"N", "CA", "C", "CB"}
        if (ref == allow_missing_cb or ref in additional_missing_cb) and "CB" not in xyz:
            needed.remove("CB")
        if not needed <= xyz.keys():
            raise SourceCompletionError(f"undefined global protein CA: {ref}")
        if "CB" in xyz:
            ca = xyz["CA"]
            volume = dot(cross(sub(xyz["N"], ca), sub(xyz["C"], ca)), sub(xyz["CB"], ca))
            if volume <= .5:
                raise SourceCompletionError(f"non-L/nearplanar global protein CA: {ref}")


def plan_source(source: bytes, residue: ComponentInstanceRef) -> CompletionRequest:
    """Rebuild complete observations, exact missing set and observed external links."""
    return _plan_source(source, residue)


def _plan_source(source: bytes, residue: ComponentInstanceRef, *,
                 other_missing_cb: frozenset[ComponentInstanceRef] = frozenset(),
                 target_capability: frozenset[str] = SOURCE_RESIDUES,
                 validate_admission: Callable[[CompletionRequest], str | None] = admission_reason) -> CompletionRequest:
    """Internal policy-specific admission; default entry points never opt in implicitly."""
    lines, atoms, _ = _layout(source)
    selected = [a for a in atoms if _component_ref(a) == residue]
    if not selected or any(a.record != "ATOM  " for a in selected):
        raise SourceCompletionError("exact ATOM target not found")
    name = selected[0].residue_name
    if name not in target_capability:
        raise SourceCompletionError("source adapter target capability remains ALA/SER/free CYS"
                                    if target_capability == SOURCE_RESIDUES else
                                    "branched source target capability is VAL/THR only"
                                    if target_capability == frozenset({"VAL", "THR"}) else
                                    "target outside explicit source capability")
    expected = CANONICAL_HEAVY_ATOMS.get(name)
    if expected is None:
        raise SourceCompletionError("target is not canonical protein")
    observations = tuple(Observation(ExactAtomRef(_component_ref(a), a.name), a.element, a.position)
                         for a in atoms)
    by_serial = {a.serial: a for a in atoms}
    pairs: dict[frozenset[ExactAtomRef], tuple[ExactAtomRef, ExactAtomRef]] = {}

    def link(a: PdbAtom, b: PdbAtom) -> None:
        left, right = ExactAtomRef(_component_ref(a), a.name), ExactAtomRef(_component_ref(b), b.name)
        if left == right:
            raise SourceCompletionError("self CONECT bond")
        if left.component != right.component:
            pairs.setdefault(frozenset((left, right)), (left, right))

    # File-order neighbouring residues only, with TER and chain-case barriers.
    groups: list[list[PdbAtom]] = []
    for atom in atoms:
        if not groups or _component_ref(groups[-1][0]) != _component_ref(atom):
            groups.append([])
        groups[-1].append(atom)
    for left, right in zip(groups, groups[1:]):
        a, b = left[0], right[0]
        if (a.chain_id != b.chain_id or a.residue_name not in CANONICAL_HEAVY_ATOMS
                or b.residue_name not in CANONICAL_HEAVY_ATOMS
                or any(l[:6].strip() == "TER" for l in lines[left[-1].line_index+1:b.line_index])):
            continue
        c = next((x for x in left if x.name == "C"), None)
        n = next((x for x in right if x.name == "N"), None)
        if c and n and 1 <= math.dist(c.position, n.position) <= 1.8:
            link(c, n)
    by_identity = {ExactAtomRef(_component_ref(a), a.name): a for a in atoms}
    disulfide_partners: dict[ExactAtomRef, ExactAtomRef] = {}

    def declared_bond(a: PdbAtom, b: PdbAtom) -> None:
        left, right = ExactAtomRef(_component_ref(a), a.name), ExactAtomRef(_component_ref(b), b.name)
        if left == right:
            raise SourceCompletionError("self declared bond")
        if left.component == right.component:
            if a.residue_name not in SOURCE_RESIDUES | target_capability:
                raise SourceCompletionError("unsupported source internal graph capability")
            allowed = {frozenset(edge) for edge in sidechain_bonds(a.residue_name)}
            if frozenset((a.name, b.name)) not in allowed:
                raise SourceCompletionError("noncanonical target/context bond")
            # Geometry of the target is checked independently by audit_candidate.
            # Metadata cannot add a new edge or remove an implicit canonical edge.
            return
        pair = frozenset((left, right))
        if pair in pairs:
            return
        if a.residue_name == b.residue_name == "CYS" and a.name == b.name == "SG":
            if residue in {left.component, right.component}:
                raise SourceCompletionError("disulfide-linked target requires coupled reconstruction")
            if not 1.8 <= math.dist(a.position, b.position) <= 2.3:
                raise SourceCompletionError("invalid observed disulfide distance")
            for x, y in ((left, right), (right, left)):
                if x in disulfide_partners and disulfide_partners[x] != y:
                    raise SourceCompletionError("multiple disulfide partners")
                disulfide_partners[x] = y
            pairs[pair] = (left, right)
            return
        raise SourceCompletionError("nonpeptide/nonlocal declared bond requires link-aware support")

    for declared in parse_declared_links(lines, atoms):
        a, b = by_identity[declared.left], by_identity[declared.right]
        if (declared.reported_distance is not None
                and abs(math.dist(a.position, b.position)-declared.reported_distance) > .1):
            raise SourceCompletionError("LINK/SSBOND distance disagrees with observed atoms")
        declared_bond(a, b)
    for line in lines:
        if line[:6] != "CONECT":
            continue
        try:
            fields = line[6:].rstrip("\n")
            values = [int(fields[i:i+5]) for i in range(0, len(fields), 5) if fields[i:i+5].strip()]
        except ValueError as exc:
            raise SourceCompletionError("invalid CONECT") from exc
        if len(values) < 2 or any(v not in by_serial for v in values):
            raise SourceCompletionError("dangling/empty CONECT")
        if len(set(values[1:])) != len(values[1:]):
            raise SourceCompletionError("CONECT bond multiplicity requires bond-order support")
        for other in values[1:]:
            a, b = by_serial[values[0]], by_serial[other]
            declared_bond(a, b)
    missing = tuple(n for n in expected if n not in {a.name for a in selected})
    request = CompletionRequest(hashlib.sha256(source).hexdigest(), residue, name,
                                observations, missing, tuple(pairs.values()))
    reason = validate_admission(request)
    if reason:
        raise SourceCompletionError(reason)
    if name == "CYS":
        ca = next(a.position for a in selected if a.name == "CA")
        # Account for two possible CA--SG reaches even if the other SG is absent;
        # coordinate proximity must not accidentally classify a disulfide as free.
        if any(a.residue_name == "CYS" and a.name == "CA" and _component_ref(a) != residue
               and math.dist(a.position, ca) < 10 for a in atoms):
            raise SourceCompletionError("nearby cysteine requires disulfide-aware admission")
    _global_ca(atoms, allow_missing_cb=residue, additional_missing_cb=other_missing_cb)
    return request


def materialize_source(source: bytes, request: CompletionRequest, candidate: Candidate, *,
                       library: SingleChiLibrary | None = None) -> MaterializedCompletion:
    """Rebind source, audit raw proposal, insert only new atoms, then audit readback."""
    authoritative = plan_source(source, request.residue)
    if authoritative != request:
        raise SourceCompletionError("source request/context/link binding mismatch")
    audit = audit_candidate(authoritative, candidate, library=library)
    if not audit.accepted:
        raise SourceCompletionError(f"raw candidate rejected: {audit.findings}")
    output = _insert_candidates(source, (request,), (candidate,))
    added = audit_materialized_source(source, request, candidate, output, library=library)
    return MaterializedCompletion(output, request.source_sha256, hashlib.sha256(output).hexdigest(),
                                  request.digest, added)


def audit_materialized_source(source: bytes, request: CompletionRequest, candidate: Candidate,
                              output: bytes, *, library: SingleChiLibrary | None = None) -> tuple[Addition, ...]:
    """Independent readback: exact old records, exact new identities and bounded rounding."""
    if plan_source(source, request.residue) != request:
        raise SourceCompletionError("source request/context/link binding mismatch")
    if not audit_candidate(request, candidate, library=library).accepted:
        raise SourceCompletionError("raw candidate rejected")
    added = _readback_additions(source, output, candidate.additions)
    rounded = replace(candidate, additions=added)
    final_audit = audit_candidate(request, rounded, library=library, pdb_roundtrip=True)
    if not final_audit.accepted:
        raise SourceCompletionError(f"rounded PDB rejected: {final_audit.findings}")
    _, readback, _ = _layout(output)
    _global_ca(readback, allow_missing_cb=None)
    return added


def _insert_candidates(source: bytes, requests: tuple[CompletionRequest, ...],
                       candidates: tuple[Candidate, ...]) -> bytes:
    """One insertion pass; allocate serials in output order, never intermediate files."""
    lines, atoms, used = _layout(source)
    insertions: dict[int, tuple[CompletionRequest, Candidate]] = {}
    for request, candidate in zip(requests, candidates, strict=True):
        insertion = max(a.line_index for a in atoms if _component_ref(a) == request.residue)+1
        while insertion < len(lines) and lines[insertion][:6] == "ANISOU":
            insertion += 1
        if insertion in insertions:
            raise SourceCompletionError("duplicate insertion target")
        insertions[insertion] = (request, candidate)
    result = []
    serial = max(used)
    for index in range(len(lines)+1):
        if index in insertions:
            request, candidate = insertions[index]
            for addition in candidate.additions:
                serial += 1
                line = _format_generated_atom(serial=serial, name=addition.atom.name, element=addition.element,
                                              position=addition.xyz, residue_name=request.residue_name,
                                              chain_id=request.residue.chain_id,
                                              sequence_number=request.residue.sequence_number,
                                              insertion_code=request.residue.insertion_code)
                result.append("ATOM  "+line[6:])
        if index < len(lines):
            result.append(lines[index])
    return "".join(result).encode("ascii")


def _readback_additions(source: bytes, output: bytes, additions: tuple[Addition, ...]) -> tuple[Addition, ...]:
    """Shared physical-record and per-coordinate rounding audit; no geometry waiver."""
    lines, atoms, used = _layout(source)
    out_lines, readback, _ = _layout(output)
    original = {ExactAtomRef(_component_ref(a), a.name) for a in atoms}
    new_indices = {a.line_index for a in readback if ExactAtomRef(_component_ref(a), a.name) not in original}
    if [line for i, line in enumerate(out_lines) if i not in new_indices] != lines:
        raise SourceCompletionError("source records changed or extra non-generated records")
    new_atoms = [a for a in readback if a.line_index in new_indices]
    if any(a.record != "ATOM  " or a.serial != max(used)+i+1
           or out_lines[a.line_index][54:66] != "  1.00  0.00"
           for i, a in enumerate(new_atoms)):
        raise SourceCompletionError("generated record/serial/occupancy/B-factor mismatch")
    added = tuple(Addition(ExactAtomRef(_component_ref(a), a.name), a.element, a.position)
                  for a in readback if ExactAtomRef(_component_ref(a), a.name) not in original)
    proposed = {a.atom: a for a in additions}
    if len(proposed) != len(additions) or len(added) != len(proposed) or {a.atom for a in added} != proposed.keys():
        raise SourceCompletionError("rounded generated ownership mismatch")
    for addition in added:
        raw = proposed[addition.atom]
        if (addition.element != raw.element
                or any(abs(a-b) > .000500001 for a, b in zip(addition.xyz, raw.xyz))):
            raise SourceCompletionError("readback exceeds PDB rounding bound")
    return added
