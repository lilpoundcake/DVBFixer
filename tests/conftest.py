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


@pytest.fixture(scope="session")
def eightcz8_dir() -> Path:
    """8CZ8 fixture — chain E carries several LYS residues truncated to
    backbone+CB by real crystallographic disorder. Exercises PDBFixer's
    addMissingAtoms seeded-retry rebuild (unseeded clash-escape MD could
    previously leave a rebuilt sidechain non-deterministic or D-chiral)."""
    return _fixture_or_skip(FIXTURES_ROOT / "8cz8")


@pytest.fixture(scope="session")
def lipid_dir() -> Path:
    """7X35 fixture — a 3-chain viral capsid (no SEQRES records) with a
    bound palmitic acid (PLM) ligand whose atom-naming convention
    (C1..C9, CA, CB, CC...) collides with standard protein backbone atom
    names. Exercises `model/renumber.py`'s HETATM-vs-gap-filled-protein
    resSeq collision (chain A has an 11-residue gap vs. its full FASTA
    sequence; PLM's naively-assigned resSeq lands exactly where the
    gap-fill needs to go) and the downstream minimize position-restore
    merge that must not cross-contaminate the two residues' coordinates."""
    return _fixture_or_skip(FIXTURES_ROOT / "lipid")


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


@pytest.fixture(scope="session")
def protein_ligand_gaps_pdb() -> Path:
    """Protein with SEQRES gaps + two real ligands (DAN, EPE) —
    exercises heterogen H-addition valence correctness."""
    return _fixture_or_skip(FIXTURES_ROOT / "protein_ligand" / "1VCU.pdb")


@pytest.fixture(scope="session")
def glycoprot_underannotated_conect_pdb() -> Path:
    """Glycoprotein with IUPAC sugar names, full header records
    (SEQRES/HELIX/SHEET/LINK/CRYST1), and CONECT/LINK annotation that's
    genuinely incomplete for one of its three real N-glycosylation
    sites — exercises `convert`'s header pass-through and its
    CONECT-vs-distance glycosidic bond detection."""
    return _fixture_or_skip(FIXTURES_ROOT / "3ry6" / "3ry6.pdb")


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


def count_h_on_atom(pdb_path: Path, resname: str, atom_name: str) -> int:
    """Count H atoms bonded to a specific named heavy atom, summed
    across every residue instance matching `resname` in the file.
    Used to check ligand valence correctness (e.g. a sulfonate O
    should have zero, a ring CH should have exactly one)."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))
    targets = {a.index for r in pdb.topology.residues() if r.name == resname
               for a in r.atoms() if a.name == atom_name}
    count = 0
    for b in pdb.topology.bonds():
        a1, a2 = b[0], b[1]
        if a1.index in targets and a2.element.symbol == "H":
            count += 1
        elif a2.index in targets and a1.element.symbol == "H":
            count += 1
    return count


def count_nonbonded_clashes(
    pdb_path: Path, cutoff_angstrom: float = 1.5,
    only_residue_names: set[str] | None = None,
) -> list[tuple[float, str, str]]:
    """Whole-structure non-bonded contact scan: any two atoms in
    DIFFERENT residues closer than ``cutoff_angstrom`` that aren't
    directly bonded to each other. A real, correctly-minimized
    structure should have none — this is what caught the glycan-tree
    drift (a tree could swing 4-10 Å off its covalent anchor into an
    unrelated, clashing part of the protein surface even though its
    OWN bond graph was perfectly correct) and the non-covalent-ligand
    splice-back clash (a ligand pasted back at pre-minimize coordinates
    into a pocket the protein had since moved into).

    If ``only_residue_names`` is given, only report a pair when at
    least one atom's residue name is in that set — use this to scope
    out the SEPARATE, already-accepted "minor local packing strain"
    tolerance from the chirality invariant's forced-reflect fallback
    (CLAUDE.md: zero D-Cα is non-negotiable, a small residual clash on
    the reflected residue is not), which is unrelated to and shouldn't
    be conflated with a genuine drift/clash regression on a specific
    residue class (e.g. glycans, ligands).
    """
    import numpy as np
    from openmm.app import PDBFile
    from openmm.unit import angstrom
    from scipy.spatial import cKDTree

    pdb = PDBFile(str(pdb_path))
    pos = np.array([list(pdb.positions[i].value_in_unit(angstrom))
                    for i in range(len(pdb.positions))])
    atoms = list(pdb.topology.atoms())
    bonded = set()
    for b in pdb.topology.bonds():
        bonded.add((b[0].index, b[1].index))
        bonded.add((b[1].index, b[0].index))

    tree = cKDTree(pos)
    violations = []
    for i, j in tree.query_pairs(cutoff_angstrom):
        a1, a2 = atoms[i], atoms[j]
        if a1.residue == a2.residue or (i, j) in bonded:
            continue
        if only_residue_names is not None and not (
            a1.residue.name in only_residue_names
            or a2.residue.name in only_residue_names
        ):
            continue
        d = float(np.linalg.norm(pos[i] - pos[j]))
        violations.append((
            d,
            f"{a1.residue.chain.id}/{a1.residue.name}{a1.residue.id}/{a1.name}",
            f"{a2.residue.chain.id}/{a2.residue.name}{a2.residue.id}/{a2.name}",
        ))
    violations.sort()
    return violations


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
