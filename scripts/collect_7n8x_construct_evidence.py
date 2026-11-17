#!/usr/bin/env python3
"""Retain primary-publication evidence; no coordinates or admission are changed."""
from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from pathlib import Path
from urllib.request import Request, urlopen

import gemmi
from scripts.mac_reconstruction_comparison import digest, save

URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8982864/fullTextXML"
DOI = "10.1371/journal.pbio.3001589"
LIMIT = 8*1024**2


def paragraphs(raw):
    root = ET.fromstring(raw)
    ids = {e.attrib.get("pub-id-type"): e.text for e in root.findall(".//article-id")}
    if ids.get("doi") != DOI:
        raise ValueError("primary publication DOI differs")
    rows = []
    needles = ("7N8X", "7N9Y", "central region", "full-length TcdB", "pHis22-TcdB")
    for element in root.iter("p"):
        text = " ".join("".join(element.itertext()).split())
        if any(needle in text for needle in needles):
            rows.append(text)
    if not any("7N8X" in text for text in rows):
        raise ValueError("publication does not explicitly identify7N8X")
    return rows


def collect(sources, output):
    ledger_path = sources / "source-results.json"
    complete = json.loads((sources / "completed.json").read_text())
    if digest(ledger_path) != complete["source_results_sha256"]:
        raise ValueError("source ledger changed")
    row = next(r for r in json.loads(ledger_path.read_text())["rows"] if r["entry_id"] == "7N8X")
    source = Path(row["source_path"])
    if digest(source) != row["source_sha256"]:
        raise ValueError("original7N8X source changed")
    document = gemmi.cif.read(str(source))
    block = document.sole_block()
    categories = {name: block.get_mmcif_category(name) for name in
                  ("_struct.", "_entity.", "_struct_ref.", "_struct_ref_seq.", "_struct_ref_seq_dif.", "_citation.")}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "declaration.json", {"entry_id": "7N8X", "source_path": str(source), "source_sha256": digest(source),
        "collector_sha256": digest(Path(__file__)), "primary_article_url": URL, "expected_DOI": DOI,
        "maximum_response_bytes": LIMIT, "timeout_seconds": 60, "HTTP_attempts_declared": 1,
        "source_or_request_changes": False, "no_coordinate_reference_or_sampling": True,
        "earlier_discovery": "publisher webfetch timed out; PMC reader succeeded; raw article now retained independently"})
    request = Request(URL, headers={"User-Agent": "dvbfixer-real-gap-construct-review/1"})
    try:
        with urlopen(request, timeout=60) as response:
            raw = response.read(LIMIT+1)
            receipt = {"url": URL, "final_url": response.url, "status": response.status,
                       "received_utc": datetime.now(UTC).isoformat()}
        with (output / "primary-publication.xml").open("xb") as stream:
            stream.write(raw)
        if len(raw) > LIMIT:
            raise ValueError("primary-publication response exceeds declared8MiB ceiling")
        excerpts = paragraphs(raw)
    except Exception as exc:
        save(output / "collection-failure.json", {"type": type(exc).__name__, "message": str(exc)})
        raise
    receipt.update(sha256=digest(output / "primary-publication.xml"), bytes=len(raw))
    save(output / "publication-receipt.json", receipt)
    save(output / "construct-evidence.json", {"entry_id": "7N8X", "source_sha256": digest(source),
        "publication_DOI": DOI, "publication_sha256": receipt["sha256"], "deposited_categories": categories,
        "deposited_label_axis_inventory": row["inventory"], "primary_publication_paragraphs": excerpts,
        "review_status": "evidence collected; scientific reviewer must distinguish sample construct from modeled region",
        "source_modified": False, "scientific_admission_complete": False})
    print(json.dumps({"entry": "7N8X", "primary_publication_retained": True,
                      "extracted_paragraphs": len(excerpts), "source_unchanged": True}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    collect(args.sources.resolve(), args.output.resolve())
