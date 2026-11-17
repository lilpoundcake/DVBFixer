"""Exact archived-code resolution for private, source-bound research declarations.

Only Python implementation bytes can resolve through an exact code snapshot.
Structures, requests, checkpoints and executable interpreters never can.
"""
from __future__ import annotations

import json
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save


def freeze_code(plan_path):
    plan = json.loads(plan_path.read_text())
    directory = plan_path.parent / "executed-code"
    directory.mkdir(exist_ok=False)
    rows = []
    for name, expected in plan["source_code_and_input_sha256"].items():
        path = Path(name)
        if path.suffix != ".py":
            continue
        if digest(path) != expected:
            raise ValueError(f"implementation changed before archival: {name}")
        snapshot = directory / f"{expected}.py"
        if not snapshot.exists():
            with snapshot.open("xb") as stream:
                stream.write(path.read_bytes())
        rows.append({"original_path": str(path), "snapshot_path": str(snapshot), "sha256": expected})
    save(directory / "index.json", {"plan_sha256": digest(plan_path), "snapshots": rows})


def verify_bindings(bindings, plans=()):
    aliases = {}
    for plan_path in plans:
        index_path = plan_path.parent / "executed-code/index.json"
        if not index_path.is_file():
            continue
        index = json.loads(index_path.read_text())
        if digest(plan_path) != index["plan_sha256"]:
            raise ValueError("code archive binds a different declaration")
        for row in index["snapshots"]:
            key = row["original_path"], row["sha256"]
            aliases[key] = Path(row["snapshot_path"])
    resolved = {}
    for name, expected in bindings.items():
        path = Path(name)
        if path.is_file() and digest(path) == expected:
            resolved[name] = path
            continue
        archived = aliases.get((name, expected)) if path.suffix == ".py" else None
        if archived is None or not archived.is_file() or digest(archived) != expected:
            raise ValueError(f"declared input changed or disappeared: {name}")
        resolved[name] = archived
    return resolved
