"""Shared pytest fixtures + fixture-file paths for dvbfixer tests.

The repository's ``test/`` directory (singular, historical) holds the raw
input PDBs used as fixtures. The ``tests/`` directory (plural, this one)
holds the pytest suite. Do not confuse the two.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES_ROOT = REPO_ROOT / "test"
GOLDEN_ROOT = Path(__file__).resolve().parent / "golden"


@pytest.fixture(scope="session")
def fixtures_root() -> Path:
    return FIXTURES_ROOT


@pytest.fixture(scope="session")
def golden_root() -> Path:
    return GOLDEN_ROOT


def _fixture_or_skip(path: Path) -> Path:
    """Return ``path`` if it exists, else skip the calling test.

    The ``test/`` directory is intentionally untracked (large binary
    fixtures live outside git). On CI where the fixtures aren't
    provisioned, tests that need them skip cleanly instead of ERROR-ing.
    Locally the tests run against the developer's own ``test/`` tree.
    """
    if not path.exists():
        pytest.skip(f"fixture missing: {path}")
    return path


@pytest.fixture(scope="session")
def small_pdb() -> Path:
    """Single-residue ASN — smallest input we ship."""
    return _fixture_or_skip(FIXTURES_ROOT / "ASN.pdb")


@pytest.fixture(scope="session")
def default_pdb() -> Path:
    """Small default PDB used across smoke tests."""
    return _fixture_or_skip(FIXTURES_ROOT / "default.pdb")


@pytest.fixture(scope="session")
def multistate_pdb() -> Path:
    """Multi-MODEL PDB — 11 states × 3 chains, exercises split's MODEL path."""
    return _fixture_or_skip(FIXTURES_ROOT / "multistate" / "test_multistate.pdb")


@pytest.fixture(scope="session")
def glycam_dir() -> Path:
    """GLYCAM-named glycoprotein fixture."""
    return _fixture_or_skip(FIXTURES_ROOT / "glycosilated_mAb_Amber_Glycam")


@pytest.fixture(scope="session")
def charmm_glycan_dir() -> Path:
    """CHARMM-GUI-named glycoprotein fixture (mirror of glycam_dir)."""
    return _fixture_or_skip(FIXTURES_ROOT / "glycosilated_mAb_Charmm-GUI")


@pytest.fixture(scope="session")
def trastuzumab_dir() -> Path:
    """Antibody fixture — templates + target FASTA for homology tests."""
    return _fixture_or_skip(FIXTURES_ROOT / "trastuzumab")


@pytest.fixture(scope="session")
def shit_dir() -> Path:
    """Curated set of raw PDBs that historically broke the zbs pipeline
    (coincident-atom H's from Modeller, chirality oscillation, etc.).
    Used by the ``@pytest.mark.slow`` end-to-end regression suite."""
    return _fixture_or_skip(FIXTURES_ROOT / "shit")


@pytest.fixture()
def tmp_workdir(tmp_path: Path) -> Path:
    """Per-test scratch dir with a stable name inside pytest's tmp_path."""
    d = tmp_path / "work"
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_pdb_atoms(path: Path) -> list[str]:
    """Return the ATOM/HETATM lines of a PDB, normalising trailing whitespace.

    Used by regression tests to compare tool output against golden files
    without being sensitive to header lines, timestamps, or CONECT ordering
    changes upstream. When a test does care about non-ATOM records, it can
    do its own comparison — this helper is the common case.
    """
    lines = []
    for raw in path.read_text().splitlines():
        if raw.startswith(("ATOM  ", "HETATM")):
            lines.append(raw.rstrip())
    return lines
