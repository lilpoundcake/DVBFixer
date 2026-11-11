#!/usr/bin/env python3
"""Independently bind planned gaps to deposited REMARK 465, not masked ATOMs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def missing_records(text: str) -> dict[tuple[str, str, str], str]:
    records = {}
    for line in text.splitlines():
        if not line.startswith("REMARK 465") or not line[21:26].strip().lstrip("-").isdigit():
            continue
        key = (line[19:20], line[21:26].strip(), line[26:27].strip())
        name = line[15:18].strip()
        if key in records and records[key] != name:
            raise ValueError(f"conflicting missing-residue annotation: {key}")
        records[key] = name
    return records


def observed_residues(text: str) -> dict[str, tuple]:
    chains: dict[str, list] = {}
    seen = set()
    for line in text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        key = (line[21], line[22:26].strip(), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(key[0], []).append((*key[1:], line[17:20].strip()))
    return {chain: tuple(residues) for chain, residues in chains.items()}


def audit(source: str, assembly: str, request: dict) -> dict:
    from dvbfixer.model.cli import AA3TO1
    from dvbfixer.model.pipeline import parse_seqres

    annotations = missing_records(source)
    original_sequences = parse_seqres(source.splitlines())
    assembled_sequences = parse_seqres(assembly.splitlines())
    source_observed, output_observed = observed_residues(source), observed_residues(assembly)
    # A transformed assembly copy has exactly the same residue sequence/author
    # identities as its source chain. Identical source chains need not be guessed:
    # all matching candidates must agree on each selected missing-residue record.
    mapping = {chain: [candidate for candidate, residues in source_observed.items()
        if residues == observed and original_sequences.get(candidate) == assembled_sequences.get(chain)]
        for chain, observed in output_observed.items()}
    gaps = []
    targets = {t["chain"]: t["sequence"] for t in request["target_sequences"]}
    for gap in request["gaps"]:
        chain = gap["chain"]
        interval = gap["target_interval"]
        sequence = targets[chain][interval["start"]:interval["stop"]]
        expected = []
        for residue in gap["generated_residues"]:
            names = {annotations.get((candidate, residue["residue_number"], residue["insertion_code"]))
                for candidate in mapping.get(chain, [])}
            if len(names) != 1 or None in names:
                raise ValueError(f"gap is not unambiguously deposited as missing: {residue}; candidates={mapping.get(chain)}")
            expected.append(AA3TO1[next(iter(names))])
        if "".join(expected) != sequence:
            raise ValueError(f"REMARK 465 and full SEQRES gap sequence differ for {chain}: {sequence}")
        gaps.append({"chain": chain, "source_chains": mapping[chain], "sequence": sequence,
                     "generated_residues": gap["generated_residues"], "kind": gap["gap_kind"]})
    return {"passed": True, "gap_count": len(gaps), "assembly_chain_sources": mapping, "gaps": gaps}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for directory in sorted((args.root / "cases").iterdir()):
        before = directory / "request-before-preparation.json"
        if not before.exists():
            rows.append({"pdb": directory.name, "passed": False, "reason": "request was not constructed"})
            continue
        try:
            report = audit((args.root / "sources" / f"{directory.name}.pdb").read_text(),
                           (directory / "selected.pdb").read_text(), json.loads(before.read_text()))
            report["request_before_sha256"] = hashlib.sha256(before.read_bytes()).hexdigest()
        except (ValueError, KeyError) as exc:
            report = {"passed": False, "reason": str(exc)}
        with (directory / "natural-gap-audit.json").open("x") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
        rows.append({"pdb": directory.name, **report})
    with (args.root / "natural-gap-audit.json").open("x") as stream:
        json.dump(rows, stream, indent=2)
        stream.write("\n")
    print(json.dumps([{k: r[k] for k in ("pdb", "passed", "gap_count", "reason") if k in r} for r in rows], indent=2))


if __name__ == "__main__":
    main()
