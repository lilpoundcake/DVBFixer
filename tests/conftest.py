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


# ---------------------------------------------------------------------------
# Input-class fixtures for the comprehensive integration matrix (0.7.5+).
# Each returns a specific PDB or directory representing one input class.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def pure_protein_small() -> Path:
    """Small pure-protein PDB — smoke fixture for prepare/minimize/protonate."""
    return _fixture_or_skip(FIXTURES_ROOT / "default.pdb")


@pytest.fixture(scope="session")
def pure_protein_antibody() -> Path:
    """Full antibody (heavy + light chains, Kabat numbering)."""
    return _fixture_or_skip(FIXTURES_ROOT / "1HZH.pdb")


@pytest.fixture(scope="session")
def glycoprot_stress() -> Path:
    """Glycoprotein stress-test PDB — trastuzumab with N-glycans on
    the Fc region. Exercises the ASN → NLN + sugar-tree paths."""
    return _fixture_or_skip(FIXTURES_ROOT / "shit" / "trastuzumab.pdb")


@pytest.fixture(scope="session")
def glycoprot_amber_glycam_dir() -> Path:
    """Glycoprotein with GLYCAM-canonical residue names (NLN + linkage
    codes like 0YB/4YB/VMA). Ready-to-load by AMBER+GLYCAM FF."""
    return _fixture_or_skip(FIXTURES_ROOT / "glycosilated_mAb_Amber_Glycam")


@pytest.fixture(scope="session")
def glycoprot_charmm_gui_dir() -> Path:
    """Glycoprotein with CHARMM-GUI 4-char residue names (BGLC, BMAN,
    etc.). Exercises the convert-to-glycam path."""
    return _fixture_or_skip(FIXTURES_ROOT / "glycosilated_mAb_Charmm-GUI")


@pytest.fixture(scope="session")
def ss_bonded_antibody() -> Path:
    """1DQJ antibody with 12 disulfide bonds — SS-preservation regression."""
    return _fixture_or_skip(FIXTURES_ROOT / "shit" / "1DQJ_original.pdb")


@pytest.fixture(scope="session")
def donor_pdb() -> Path:
    """Transplant donor PDB."""
    return _fixture_or_skip(FIXTURES_ROOT / "donor.pdb")


@pytest.fixture(scope="session")
def shit_broken_geom_pdbs() -> list[Path]:
    """The four historically-broken PDBs used by test_zbs_shit_inputs."""
    names = ("1EMV_original.pdb", "1FR2_original.pdb",
             "2VLN_original.pdb", "2VLQ_original.pdb")
    return [_fixture_or_skip(FIXTURES_ROOT / "shit" / n) for n in names]


# ---------------------------------------------------------------------------
# Helpers used by integration tests.
# ---------------------------------------------------------------------------


def count_d_ca_residues(pdb_path: Path) -> int:
    """Return the number of D-Cα residues in a PDB file. Zero means the
    chirality invariant is satisfied — this is what every dvbfixer
    output is REQUIRED to guarantee after minimize's force-reflect."""
    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import find_d_residues
    pdb = PDBFile(str(pdb_path))
    return len(find_d_residues(pdb.topology, pdb.positions))


def count_ss_bonds(pdb_path: Path) -> int:
    """Count SG-SG bonds in the loaded topology."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))
    return sum(1 for b in pdb.topology.bonds()
               if b[0].name == "SG" and b[1].name == "SG")


def count_residues_by_name(pdb_path: Path, names: set[str]) -> int:
    """Count residues whose name is in ``names``."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))
    return sum(1 for r in pdb.topology.residues() if r.name in names)


_SUGAR_NAMES = frozenset({
    "BGL", "BGC", "GLC", "GAL", "BGA", "MAN", "AMA", "BMA",
    "NAG", "NDG", "NGA", "A2G", "FUC", "FUL", "AFU", "SIA",
    "BGLC", "AGLC", "BMAN", "AMAN", "BGAL", "AGAL", "BFUC", "AFUC",
    "BGLCNA", "AGLCNA",
})

_AMBER_VARIANT_NAMES = frozenset({
    "HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN",
})


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
