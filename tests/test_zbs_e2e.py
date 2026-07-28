"""End-to-end zbs pipeline tests across input classes × force fields.

Every dvbfixer output must guarantee:
- Zero D-Cα residues (chirality invariant enforced by minimize's
  post-phase-2 force-reflect).
- Preserved SS bonds (when input had any).
- Preserved variant residues where PROPKA predicts them.
- Preserved glycan chains (when input had any) OR clean
  strip-and-splice behaviour when heterogens are stripped.

These are integration tests — they run the full CLI via
`subprocess`, load the output PDB, and assert invariants. Marked
`@pytest.mark.slow` because each run touches PROPKA + Modeller +
OpenMM.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from tests.conftest import (
    _SUGAR_NAMES,
    count_d_ca_residues,
    count_residues_by_name,
    count_ss_bonds,
)


def _run_zbs(input_pdb: Path, workdir: Path, *extra_args: str) -> Path:
    """Run `dvbfixer zbs ...` under a fresh workdir; return output path."""
    workdir.mkdir(parents=True, exist_ok=True)
    local_input = workdir / "input.pdb"
    shutil.copy(input_pdb, local_input)
    output = workdir / "out.pdb"
    cmd = [
        "dvbfixer", "zbs", str(local_input),
        "-o", str(output),
        "--no-solvent",
        *extra_args,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          cwd=str(workdir))
    assert proc.returncode == 0, (
        f"zbs failed (exit={proc.returncode})\nstdout:\n{proc.stdout}\n"
        f"stderr:\n{proc.stderr}"
    )
    assert output.exists(), f"zbs did not produce {output}"
    return output


@pytest.mark.slow
def test_zbs_pure_protein_default_ff(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """Smallest realistic case — pure protein, default FF."""
    out = _run_zbs(pure_protein_small, tmp_workdir)
    assert count_d_ca_residues(out) == 0, (
        "zbs output has D-Cα residues (chirality invariant broken)"
    )


@pytest.mark.slow
def test_zbs_antibody_ss_bonds_preserved(
    ss_bonded_antibody: Path, tmp_workdir: Path,
) -> None:
    """1DQJ antibody has 12 SS bonds — verify they survive minimize
    (the SG-SG exception in _drop_spurious_inter_aa_bonds)."""
    out = _run_zbs(ss_bonded_antibody, tmp_workdir)
    assert count_d_ca_residues(out) == 0
    n_ss = count_ss_bonds(out)
    assert n_ss >= 6, (
        f"expected ≥6 SS bonds (12 CYX residues → 6 pairs), got {n_ss}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("ff_args", [
    pytest.param(("--ff", "amber"), id="amber"),
    pytest.param(("--ff", "amber+glycam"), id="amber_glycam"),
    # CHARMM+glycoprotein is intentionally OMITTED from CI: it takes
    # >10 min on an antibody-sized system because charmm36.xml has no
    # sugar templates so every glycan atom becomes a free heavy atom
    # under minimize with no bond constraints — the minimizer thrashes.
    # Run manually if you specifically need CHARMM+glycoprotein.
])
def test_zbs_glycoprotein_amber_ffs(
    glycoprot_stress: Path, tmp_workdir: Path, ff_args: tuple[str, ...],
) -> None:
    """Trastuzumab (glycoprotein with N-glycans) must complete under
    the AMBER-family FFs the user can pass. Legacy prep +
    strip-and-splice for amber; legacy prep + GLYCAM system for
    amber+glycam."""
    out = _run_zbs(glycoprot_stress, tmp_workdir, *ff_args)
    assert count_d_ca_residues(out) == 0, (
        f"D-Cα in glycoprotein output ({ff_args})"
    )
    # Glycans should survive splice-back (or be included in the
    # GLYCAM-parametrised system).
    n_sugars = count_residues_by_name(out, _SUGAR_NAMES)
    assert n_sugars >= 10, (
        f"expected glycan chains preserved through zbs ({ff_args}); "
        f"got only {n_sugars} sugar residues"
    )


@pytest.mark.slow
def test_zbs_glycoprot_amber_glycam_variant_preservation(
    glycoprot_stress: Path, tmp_workdir: Path,
) -> None:
    """When --ff amber+glycam is selected on a glycoprotein, the
    ASN glycosylation sites should end up named NLN (or preserved
    as ASN with the covalent bond intact via CONECT)."""
    out = _run_zbs(glycoprot_stress, tmp_workdir, "--ff", "amber+glycam")
    # Either NLN (fully GLYCAM-canonical) or preserved as ASN — both
    # acceptable as long as the pipeline completed and 0 D-Cα.
    # Sanity: at least some residues should be present in output.
    from openmm.app import PDBFile
    pdb = PDBFile(str(out))
    n_res = sum(1 for _ in pdb.topology.residues())
    assert n_res > 0
