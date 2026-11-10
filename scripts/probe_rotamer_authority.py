"""Retain data-authority availability/terms evidence; never download statistics.

The old public mirror's license is evidence to review, not accepted terms or a
license grant for a future artifact. A page fetch never clears the data gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from urllib.request import urlopen

PAGES = {
    "publisher": "https://dunbrack.fccc.edu/lab/bbdep2010",
    "tutorial": "https://dunbrack.fccc.edu/lab/bbdep2010_tutorial",
    "linked-terms": "https://dunbrack.fccc.edu/lab/bbdep2010_lic",
    "legacy-mirror-terms": "https://graylab.jhu.edu/Sergey/for.Rocco/nobu-3/ExtendedOpt1-5/AcademicLicenseBBDep2010.txt",
    "legacy-mirror-description": "https://graylab.jhu.edu/Sergey/for.Rocco/nobu-3/ExtendedOpt1-5/Description.txt",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)

    def get(item):
        name, url = item
        record = {"url": url}
        try:
            with urlopen(url, timeout=10) as response:
                body = response.read(512*1024+1)
                if len(body) > 512*1024:
                    raise ValueError("authority page exceeds bound")
                record.update({"status": "retrieved", "final_url": response.url,
                               "sha256": hashlib.sha256(body).hexdigest(), "bytes": len(body)})
            (args.output/f"{name}.txt").write_bytes(body)
        except (OSError, ValueError) as exc:
            record.update({"status": "unavailable", "error": f"{type(exc).__name__}: {exc}"})
        return name, record

    with ThreadPoolExecutor(max_workers=len(PAGES)) as pool:
        pages = dict(pool.map(get, PAGES.items()))
    result = {"schema": "dvbfixer.rotamer-authority-probe.v1", "checked_at": datetime.now(UTC).isoformat(),
              "python": platform.python_version(), "pages": pages,
              "statistics_downloaded": False, "license_accepted": False,
              "artifact_redistribution_cleared": False,
              "note": "Publisher search results report CC BY 4.0 from 2019; the old mirror ships restrictive academic terms. No statistics authority selected or relicensed."}
    (args.output/"evidence.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
