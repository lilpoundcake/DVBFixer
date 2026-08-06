"""Integrity checks for the tracked structural fixture set."""

from __future__ import annotations

import hashlib
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures"
REPO_ROOT = Path(__file__).resolve().parent.parent


def test_fixture_manifest_is_complete_and_matches_checksums() -> None:
    manifest = FIXTURES / "MANIFEST.sha256"
    entries = []
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = REPO_ROOT / relative
        assert path.is_file(), f"tracked fixture missing: {relative}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"fixture checksum mismatch: {relative}"
        entries.append(path.resolve())

    fixture_inputs = {
        path.resolve() for path in FIXTURES.rglob("*")
        if path.is_file() and path.suffix.lower() in {".pdb", ".fasta"}
    }
    assert set(entries) == fixture_inputs, "fixture manifest has missing or stale entries"
