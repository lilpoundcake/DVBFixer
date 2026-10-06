"""Read-only CCD component information and structure-completeness reporting."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from dvbfixer.ccd import (
    CcdComponent,
    CcdError,
    fetch_ccd_component,
    load_ccd_component,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dvbfixer component-info",
        description=(
            "Query one authoritative wwPDB CCD component and optionally compare "
            "its heavy-atom inventory with every matching PDB residue."
        ),
    )
    parser.add_argument("component", help="CCD component ID, for example FAD or LBN")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--ccd-cif", help="Explicit local CCD component CIF")
    source.add_argument(
        "--online",
        action="store_true",
        help="Explicitly allow a bounded HTTPS request to files.rcsb.org",
    )
    parser.add_argument("--structure", help="Optional PDB to compare by residue name")
    parser.add_argument("--cache-dir", help="Private online CCD cache directory")
    parser.add_argument("--timeout", type=float, default=15.0,
                        help="Online request timeout in seconds (default: 15)")
    parser.add_argument("--json", dest="json_path", metavar="PATH",
                        help="Write machine-readable JSON; use '-' for stdout")
    from dvbfixer.batch import add_runtime_help

    add_runtime_help(parser)
    return parser.parse_args(argv)


def _structure_inventory(path: Path, component: CcdComponent) -> list[dict[str, Any]]:
    expected = {atom.name: atom.element for atom in component.heavy_atoms}
    residues: dict[tuple[int, str, str, str, int], dict[str, str]] = defaultdict(dict)
    model = 1
    previous_key: tuple[int, str, str, str] | None = None
    active_occurrence = 1
    occurrence_by_key: dict[tuple[int, str, str, str], int] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.startswith("MODEL "):
            try:
                model = int(line[10:14].strip())
            except ValueError as exc:
                raise CcdError(f"invalid MODEL record on line {line_number}") from exc
            previous_key = None
            continue
        if line.startswith("TER"):
            previous_key = None
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 78:
            raise CcdError(f"truncated PDB coordinate record on line {line_number}")
        residue_key = (model, line[21], line[22:26].strip(), line[26].strip())
        if residue_key != previous_key:
            active_occurrence = occurrence_by_key.get(residue_key, 0) + 1
            occurrence_by_key[residue_key] = active_occurrence
        previous_key = residue_key
        if line[17:20].strip() != component.component_id:
            continue
        altloc = line[16].strip()
        if altloc not in {"", "A"}:
            raise CcdError(
                f"component {component.component_id} has unsupported alternate locations"
            )
        key = (*residue_key, active_occurrence)
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if element == "H":
            continue
        if not element:
            raise CcdError(f"component atom {name!r} lacks an explicit PDB element")
        previous = residues[key].setdefault(name, element)
        if previous != element:
            raise CcdError(f"component atom {name!r} has inconsistent elements")
    result: list[dict[str, Any]] = []
    for (model, chain, residue_number, insertion_code, occurrence), observed in sorted(residues.items()):
        missing = sorted(set(expected) - set(observed))
        unexpected = sorted(set(observed) - set(expected))
        element_mismatches = sorted(
            name
            for name in set(expected) & set(observed)
            if expected[name] != observed[name]
        )
        result.append({
            "model": model,
            "chain": chain,
            "residue_number": residue_number,
            "insertion_code": insertion_code,
            "occurrence": occurrence,
            "observed_heavy_atom_count": len(observed),
            "expected_heavy_atom_count": len(expected),
            "missing_heavy_atoms": missing,
            "unexpected_heavy_atoms": unexpected,
            "element_mismatches": element_mismatches,
            "complete": not missing and not unexpected and not element_mismatches,
        })
    return result


def _human_report(component: CcdComponent, instances: list[dict[str, Any]]) -> str:
    lines = [
        f"Component: {component.component_id} — {component.name}",
        f"Type: {component.component_type}",
        f"Formula: {component.formula}",
        f"Formal charge: {component.formal_charge if component.formal_charge is not None else 'unknown'}",
        f"Atoms: {len(component.atoms)} total, {len(component.heavy_atoms)} heavy",
        f"Bonds: {len(component.bonds)}",
        f"Provisional DVBfixer class: {component.provisional_class}",
        f"Authority: {component.provenance.mode} sha256:{component.provenance.sha256}",
    ]
    if component.provisional_class == "A":
        lines.append(
            "Class note: provisional Class A assumes no external covalent links in the structure."
        )
    elif component.provisional_class == "C":
        lines.append(
            "Class note: cofactor/metal chemistry requires explicit state and coordination evidence."
        )
    elif component.ambiguous:
        lines.append("Class note: the CCD record is marked chemically ambiguous.")
    if instances:
        lines.append("Structure instances:")
        for item in instances:
            identity = (
                f"model {item['model']} "
                f"{item['chain'] or '<blank>'}:{item['residue_number']}"
                f"{item['insertion_code']} occurrence {item['occurrence']}"
            )
            lines.append(
                f"  {identity}: {item['observed_heavy_atom_count']}/"
                f"{item['expected_heavy_atom_count']} heavy atoms; "
                + ("complete" if item["complete"] else "incomplete")
            )
            if item["missing_heavy_atoms"]:
                lines.append("    missing: " + ", ".join(item["missing_heavy_atoms"]))
            if item["unexpected_heavy_atoms"]:
                lines.append("    unexpected: " + ", ".join(item["unexpected_heavy_atoms"]))
            if item["element_mismatches"]:
                lines.append("    element mismatch: " + ", ".join(item["element_mismatches"]))
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        component = (
            fetch_ccd_component(
                args.component,
                cache_dir=Path(args.cache_dir) if args.cache_dir else None,
                timeout_seconds=args.timeout,
            )
            if args.online
            else load_ccd_component(Path(args.ccd_cif), args.component)
        )
        instances = (
            _structure_inventory(Path(args.structure).expanduser().resolve(), component)
            if args.structure
            else []
        )
    except (CcdError, OSError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    payload = component.to_dict()
    payload["structure_instances"] = instances
    if args.json_path:
        text = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
        if args.json_path == "-":
            sys.stdout.write(text)
        else:
            Path(args.json_path).write_text(text, encoding="utf-8")
            print(f"Wrote {args.json_path}")
    else:
        sys.stdout.write(_human_report(component, instances))
