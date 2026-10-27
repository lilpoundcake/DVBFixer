"""Pure PDB alternate-state analysis and residue-coherent selection.

Coordinates, occupancy and serials are preserved. This is selection, not repair.
The CLI owns CIF normalization and filesystem publication.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any


class AltlocError(ValueError):
    """A selection cannot be materialized without guessing atom identity."""


@dataclass(frozen=True)
class Atom:
    line: int
    model: int
    serial: int
    key: str
    name: str
    label: str
    resname: str
    occupancy: float | None
    heavy: bool


def parse_atoms(text: str) -> tuple[list[Atom], dict[str, dict[str, Any]], dict[int, int]]:
    atoms: list[Atom] = []
    groups: dict[str, dict[str, Any]] = {}
    contexts: dict[int, int] = {}
    model, ordinal = 1, 0
    previous = None
    occurrences: dict[tuple[int, str, str, str], int] = {}
    serials: set[tuple[int, int]] = set()
    for index, line in enumerate(text.splitlines(keepends=True)):
        record = line[:6].strip()
        if record == "MODEL":
            ordinal += 1
            model, previous = ordinal, None
        contexts[index] = model
        if record in {"TER", "ENDMDL"}:
            previous = None
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(line.rstrip("\r\n")) < 54:
            raise AltlocError(f"truncated coordinate record at line {index + 1}")
        try:
            serial = int(line[6:11])
            int(line[22:26])
            coordinates = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            occupancy = float(line[54:60]) if line[54:60].strip() else None
        except ValueError as exc:
            raise AltlocError(f"invalid atom fields at line {index + 1}") from exc
        if not all(math.isfinite(value) for value in coordinates) or (
            occupancy is not None and (not math.isfinite(occupancy) or not 0 <= occupancy <= 1)
        ):
            raise AltlocError(f"invalid coordinate/occupancy at line {index + 1}")
        if (model, serial) in serials:
            raise AltlocError(f"duplicate atom serial {serial} in model {model}")
        if not 1 <= serial <= 99999:
            raise AltlocError(f"atom serial outside PDB limits at line {index + 1}")
        serials.add((model, serial))
        base = (model, line[21], line[22:26].strip(), line[26].strip())
        if base != previous:
            occurrences[base] = occurrences.get(base, 0) + 1
        previous = base
        key = json.dumps([*base, occurrences[base]], separators=(",", ":"))
        group = groups.setdefault(key, {
            "key": key, "model": model, "chain": base[1], "residueNumber": base[2],
            "insertionCode": base[3], "occurrence": occurrences[base],
        })
        name, label, resname = line[12:16].strip(), line[16].strip(), line[17:20].strip()
        if not name or not resname:
            raise AltlocError(f"missing atom/residue name at line {index + 1}")
        element = line[76:78].strip().upper() or name.lstrip("0123456789")[:1].upper()
        atoms.append(Atom(index, model, serial, group["key"], name, label, resname,
                          occupancy, element not in {"H", "D"}))
    if not atoms:
        raise AltlocError("input contains no atoms")
    return atoms, groups, contexts


def select_altlocs(
    text: str, *, source_sha256: str, policy: str = "occupancy",
    prefer: str | None = None, fallback: str = "error",
    overrides: dict[str, str] | None = None,
) -> tuple[str | None, dict[str, Any]]:
    """Analyze and, when all choices resolve, produce a single-state PDB."""
    if policy not in {"manual", "occupancy", "first", "prefer"}:
        raise AltlocError("unknown selection policy")
    if prefer is not None and policy != "prefer":
        raise AltlocError("--prefer requires --policy prefer")
    if fallback not in {"error", "occupancy", "first"}:
        raise AltlocError("unknown preferred-label fallback")
    if policy == "prefer" and (not prefer or len(prefer) != 1 or prefer.isspace()):
        raise AltlocError("prefer policy requires one nonblank --prefer label")
    atoms, groups, contexts = parse_atoms(text)
    by_group: dict[str, list[Atom]] = {}
    for atom in atoms:
        by_group.setdefault(atom.key, []).append(atom)
    overrides = overrides or {}
    alternatives = {key for key, records in by_group.items() if any(atom.label for atom in records)}
    if set(overrides) - alternatives:
        raise AltlocError("selection names a missing or non-alternate residue")
    rows = []
    selected_lines: set[int] = set()
    warnings: list[str] = []
    unresolved = []
    for key, records in by_group.items():
        shared = [atom for atom in records if not atom.label]
        labels = list(dict.fromkeys(atom.label for atom in records if atom.label))
        if len({atom.name for atom in shared}) != len(shared):
            raise AltlocError(f"duplicate shared atom identity in residue {key}")
        if not labels:
            if len({atom.resname for atom in records}) != 1:
                raise AltlocError(f"conflicting unlabelled residue names in {key}")
            selected_lines.update(atom.line for atom in records)
            continue
        states = []
        row_warnings = []
        for label in labels:
            labeled = [atom for atom in records if atom.label == label]
            combined = shared + labeled
            if len({atom.name for atom in combined}) != len(combined):
                raise AltlocError(f"shared/labeled or duplicate atom collision in residue {key}, state {label}")
            if len({atom.resname for atom in combined}) != 1:
                raise AltlocError(f"incompatible shared residue names in {key}, state {label}")
            heavy = [atom for atom in labeled if atom.heavy]
            values = [atom.occupancy for atom in heavy]
            score = sum(values) / len(values) if values and all(v is not None for v in values) else None
            states.append({
                "label": label, "residueName": labeled[0].resname,
                "occupancy": score, "atomNames": [atom.name for atom in labeled],
                "serials": [atom.serial for atom in labeled], "atomCount": len(labeled),
            })
        if len({tuple(sorted(state["atomNames"])) for state in states}) > 1:
            row_warnings.append("State atom inventories differ; selection does not rebuild missing atoms.")
        if len({state["residueName"] for state in states}) > 1:
            row_warnings.append("Microheterogeneity: state residue names differ; inspect sequence metadata after selection.")
        if any(state["occupancy"] is None for state in states):
            row_warnings.append("Missing heavy-atom occupancy; occupancy policy needs manual selection or explicit first fallback.")
        choice = overrides.get(key)
        if choice is not None and choice not in labels:
            raise AltlocError(f"selected state {choice!r} is absent in residue {key}")
        method = "override" if choice else policy
        effective = policy
        if not choice and policy == "prefer":
            if prefer in labels:
                choice = prefer
            else:
                row_warnings.append(f"Preferred state {prefer!r} is absent; fallback is {fallback}.")
                effective = fallback
        if not choice and effective == "first":
            choice = labels[0]
        if not choice and effective == "occupancy":
            if all(state["occupancy"] is not None for state in states):
                best = max(state["occupancy"] for state in states)
                tied = [state["label"] for state in states if abs(state["occupancy"] - best) < 1e-9]
                choice = tied[0]
                if len(tied) > 1:
                    row_warnings.append("Equal occupancy: chose first encountered state; review manually.")
        if not choice:
            unresolved.append(key)
        else:
            selected_lines.update(atom.line for atom in records if atom.label in {"", choice})
        row = {**groups[key], "states": states, "selected": choice, "method": method,
               "sharedSerials": [atom.serial for atom in shared], "warnings": row_warnings}
        rows.append(row)
        warnings.extend(f"{key}: {message}" for message in row_warnings)
    report: dict[str, Any] = {
        "schemaVersion": 1, "sourceSha256": source_sha256,
        "normalizedSha256": hashlib.sha256(text.encode()).hexdigest(),
        "policy": policy, "prefer": prefer, "fallback": fallback, "overrides": overrides,
        "residues": rows, "warnings": warnings, "unresolved": unresolved,
        "inputAtoms": len(atoms), "removedAtoms": None, "outputSha256": None,
    }
    if unresolved:
        return None, report
    output, record_warnings = _materialize(text, atoms, contexts, selected_lines)
    report["warnings"].extend(record_warnings)
    report["removedAtoms"] = len(atoms) - len(selected_lines)
    report["outputSha256"] = hashlib.sha256(output.encode()).hexdigest()
    return output, report


def _materialize(text: str, atoms: list[Atom], contexts: dict[int, int], selected: set[int]) -> tuple[str, list[str]]:
    source = {(atom.model, atom.serial): atom for atom in atoms}
    retained = {key: atom for key, atom in source.items() if atom.line in selected}
    by_line = {atom.line: atom for atom in atoms}
    models = sorted({atom.model for atom in atoms})
    explicit_models = any(line.startswith("MODEL ") for line in text.splitlines())
    active_model = False
    output = []
    warnings = []
    source_lines = text.splitlines(keepends=True)
    previous_retained: Atom | None = None
    for index, raw in enumerate(source_lines):
        line = raw.rstrip("\r\n")
        record = line[:6].strip()
        model = contexts[index]
        if record == "MODEL":
            active_model = True
        if record == "ENDMDL":
            active_model = False
        applicable = [model] if active_model or not explicit_models else models
        if record in {"MASTER"}:
            continue
        if index in by_line:
            if index not in selected:
                continue
            previous_retained = by_line[index]
            line = line[:16] + " " + line[17:]
        elif record == "TER":
            if previous_retained and previous_retained.model == model and len(line) >= 27:
                atom_line = source_lines[previous_retained.line]
                if line[21:27] == atom_line[21:27]:
                    line = line[:17] + previous_retained.resname.rjust(3) + line[20:]
            previous_retained = None
        elif record in {"ANISOU", "SIGATM", "SIGUIJ"}:
            try:
                atom = source[(model, int(line[6:11]))]
            except (KeyError, ValueError) as exc:
                raise AltlocError(f"unmatched {record} at line {index + 1}") from exc
            original = source_lines[atom.line]
            if line[12:27] != original[12:27]:
                raise AltlocError(f"{record} identity mismatch at line {index + 1}")
            if atom.line not in selected:
                continue
            line = line[:16] + " " + line[17:]
        elif record == "CONECT":
            try:
                values = [int(line[start:start + 5]) for start in range(6, len(line), 5) if line[start:start + 5].strip()]
            except ValueError as exc:
                raise AltlocError(f"malformed CONECT at line {index + 1}") from exc
            variants = []
            for current in applicable:
                if len(values) < 2 or any((current, serial) not in source for serial in values):
                    variants.append(None)
                    continue
                kept = [serial for serial in values[1:] if (current, serial) in retained]
                variants.append(tuple([values[0], *kept]) if (current, values[0]) in retained and kept else None)
            if len(set(variants)) != 1:
                raise AltlocError("global CONECT requires incompatible selections across MODELs")
            if variants[0] is None:
                warnings.append(f"Dropped CONECT at line {index + 1} with removed/dangling endpoints.")
                continue
            line = "CONECT" + "".join(f"{serial:5d}" for serial in variants[0])
        elif record == "LINK":
            if len(line) < 57:
                raise AltlocError(f"truncated LINK at line {index + 1}")
            decisions = []
            for current in applicable:
                valid = True
                for name, label, resname, chain, number, icode in (
                    (line[12:16].strip(), line[16].strip(), line[17:20].strip(), line[21], line[22:26].strip(), line[26].strip()),
                    (line[42:46].strip(), line[46].strip(), line[47:50].strip(), line[51], line[52:56].strip(), line[56].strip()),
                ):
                    matches = [a for a in atoms if a.model == current and a.name == name
                               and a.resname == resname and (not label or a.label == label)
                               and json.loads(a.key)[1:4] == [chain, number, icode] and a.line in selected]
                    if len(matches) > 1:
                        raise AltlocError("LINK endpoint matches repeated residue occurrences")
                    valid &= len(matches) == 1
                decisions.append(valid)
            if len(set(decisions)) != 1:
                raise AltlocError("global LINK requires incompatible selections across MODELs")
            if not decisions[0]:
                warnings.append(f"Dropped LINK at line {index + 1} referencing a removed state.")
                continue
            line = line[:16] + " " + line[17:46] + " " + line[47:]
        elif record == "SSBOND":
            if len(line) < 36:
                raise AltlocError(f"truncated SSBOND at line {index + 1}")
            decisions = []
            for current in applicable:
                valid = True
                for chain, number, icode in ((line[15], line[17:21].strip(), line[21].strip()),
                                              (line[29], line[31:35].strip(), line[35].strip())):
                    matches = [a for a in retained.values() if a.model == current and a.name == "SG"
                               and a.resname in {"CYS", "CYX", "CYM"}
                               and json.loads(a.key)[1:4] == [chain, number, icode]]
                    if len(matches) > 1:
                        raise AltlocError("SSBOND endpoint matches repeated residue occurrences")
                    valid &= len(matches) == 1
                decisions.append(valid)
            if len(set(decisions)) != 1:
                raise AltlocError("global SSBOND requires incompatible selections across MODELs")
            if not decisions[0]:
                warnings.append(f"Dropped SSBOND at line {index + 1} without a retained sulfur endpoint.")
                continue
        output.append(line + "\n")
    return "".join(output), warnings
