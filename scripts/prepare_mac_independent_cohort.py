#!/usr/bin/env python3
"""Prospective independent-structure mask admission only; no reserved-seed runs.

Two named unrelated protein families, full unique protein chain, fixed ordinal
masks. No endpoint trimming, occupancy edits, atom repair or case replacement.
The exact validation policy/seeds must be frozen separately before sampling.
"""
from __future__ import annotations

import argparse
import collections
import json
import traceback
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_mac_comparison_expansion import prepare

ENTRIES = (("1PGA", "protein-G B1 domain"), ("2CI2", "chymotrypsin inhibitor 2"))
MASKS = (("internal4", 9, 13), ("internal8", 23, 31), ("Nterminal4", 0, 4))


def run(output: Path, *, retained_sources: Path | None = None, preserve_terminal_oxt=False):
    output.mkdir(parents=True, exist_ok=False)
    save(output / "protocol.json", {"entries": ENTRIES, "ordinal_masks_zero_based_half_open": MASKS,
        "projection": "full unique ATOM protein chain; no trim-to-pass or masked-error-driven selection",
        "purpose": "input/reference admission, not sampling or policy tuning", "new_sampling_calls": 0,
        "source_acquisition": "RCSB public experimental PDB download; no sequence upload or MSA service",
        "unrelated_to_development": "distinct structural-family labels, not a sequence-cluster or model-training-disjointness proof",
        "observed_terminal_oxt": "preserve exact true-terminal observed OXT; never generated/deleted" if preserve_terminal_oxt else "original strict canonical-internal inventory",
        "script_sha256": digest(Path(__file__)), "qualification": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    sources = output / "sources"
    sources.mkdir()
    rows = []
    project = Path(__file__).resolve().parents[1]
    for pdb, family in ENTRIES:
        source = sources / f"{pdb}.pdb"
        url = f"https://files.rcsb.org/download/{pdb}.pdb"
        try:
            if retained_sources is None:
                with urllib.request.urlopen(url, timeout=60) as response:
                    body = response.read(2*1024**2+1)
                    if len(body) > 2*1024**2:
                        raise ValueError("source exceeds the declared 2 MiB PDB download bound")
                    headers = dict(response.headers)
            else:
                body = (retained_sources / f"{pdb}.pdb").read_bytes()
                headers = {"retained_source": str(retained_sources / f"{pdb}.pdb"), "sha256": digest(retained_sources / f"{pdb}.pdb")}
            source.write_bytes(body)
            lines = body.decode().splitlines(keepends=True)
            if any(line.startswith(("MODEL ", "SSBOND", "LINK  ")) for line in lines):
                raise ValueError("this independent admission slice excludes multiple models and explicit covalent links")
            chains = {line[21] for line in lines if line.startswith("ATOM  ")}
            if len(chains) != 1:
                raise ValueError("independent slice requires a unique protein chain, not an automatically selected chain")
            chain = next(iter(chains))
            numbers = list(dict.fromkeys(int(line[22:26]) for line in lines if line.startswith("ATOM  ")))
            if numbers != list(range(numbers[0], numbers[-1]+1)):
                raise ValueError("full observed chain is not numerically continuous")
            save(sources / f"{pdb}.download.json", {"url": url, "sha256": digest(source), "http_headers": headers,
                "downloaded_utc": datetime.now(UTC).isoformat(), "family_label": family,
                "explicit_projection": "protein-only; all HETATM omitted, no claim about complete entry",
                "omitted_heterogen_atom_counts": dict(collections.Counter(line[17:20].strip() for line in lines if line.startswith("HETATM"))),
                "first": numbers[0], "last": numbers[-1], "chain": chain})
            for label, start, stop in MASKS:
                case_id = f"{pdb}-{chain}{numbers[0]}-{numbers[-1]}-{label}"
                spec = {"id": case_id, "stratum": label, "family_label": family,
                    "fixture": str(source.relative_to(project)), "spans": [[chain, numbers[0], numbers[-1]]],
                    "hidden": [[chain, numbers[0]+start, numbers[0]+stop-1]]}
                try:
                    manifest = prepare(spec, output / case_id, allow_observed_terminal_oxt=preserve_terminal_oxt)
                    rows.append({"case": case_id, "pdb": pdb, "family_label": family, "mask_admitted": True,
                        "manifest": manifest, "reserved_seed_sampling_started": False})
                except Exception:
                    error = traceback.format_exc()
                    save(output / case_id / "refusal.json", {"error": error, "sampling": False})
                    rows.append({"case": case_id, "pdb": pdb, "mask_admitted": False, "error": error})
        except Exception:
            error = traceback.format_exc()
            save(sources / f"{pdb}.refusal.json", {"error": error, "sampling": False})
            rows.append({"pdb": pdb, "entry_admitted": False, "error": error})
    save(output / "results.json", rows)
    print(json.dumps([{k: v for k, v in row.items() if k != "manifest"} for row in rows], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--retained-sources", type=Path)
    parser.add_argument("--preserve-terminal-oxt", action="store_true")
    args = parser.parse_args()
    run(args.output.resolve(), retained_sources=None if args.retained_sources is None else args.retained_sources.resolve(),
        preserve_terminal_oxt=args.preserve_terminal_oxt)
