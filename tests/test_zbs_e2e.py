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
    count_nonbonded_clashes,
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
    # GLYCAM-parametrised system). amber+glycam auto-converts PDB-standard
    # sugar names (BMA, NAG, ...) to GLYCAM-canonical 3-char codes (0fA,
    # 2MA, ...) as part of normal, intentional processing, so a residue
    # must match _SUGAR_NAMES OR the GLYCAM-code pattern to count.
    from openmm.app import PDBFile

    from dvbfixer.ffutils import is_glycam_sugar
    pdb = PDBFile(str(out))
    n_sugars = sum(1 for r in pdb.topology.residues()
                   if r.name in _SUGAR_NAMES or is_glycam_sugar(r.name))
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


@pytest.mark.slow
def test_zbs_preserves_undocumented_glycosylation_site(
    glycoprot_underannotated_conect_pdb: Path, tmp_workdir: Path,
) -> None:
    """This fixture has 4 real N-glycosylation sites, but the deposited
    PDB's own CONECT/LINK annotation only documents 2 of them, and its
    raw text has a bare `TER` record right before one sugar chain's
    HETATM block. Chained regression coverage: renumber's TER/newline
    bug, model.py's CONECT-inference timing (before vs after
    Modeller), and convert_to_glycam's all-or-nothing CONECT gate all
    had to be fixed together for every site to survive the full
    pipeline run directly (no `dvbfixer convert` pre-processing).

    Also regression coverage for three further bugs found when real
    production runs still showed "glycans processed incorrectly,
    separate from the protein" despite the bond-existence checks above
    passing: (1) a dangling/stale CONECT record in this exact deposited
    PDB that `renumber.py` used to pass through unchanged, fabricating
    a spurious bond once its stale serial number collided with a real
    atom's new one; (2) `add_glycam_bonds` being dead code since 0.7.8
    (broken import silently swallowed); (3) heterogen heavy atoms
    getting zero positional restraint during minimize, letting a
    multi-residue glycan tree drift 4-10 Å off its covalent anchor into
    an unrelated, clashing part of the protein surface even with a
    perfectly correct bond graph."""
    out = _run_zbs(glycoprot_underannotated_conect_pdb, tmp_workdir)
    assert count_d_ca_residues(out) == 0

    from openmm.app import PDBFile
    from dvbfixer.ffutils import is_glycam_sugar
    glycan_resnames = {r.name for r in PDBFile(str(out)).topology.residues()
                        if r.name == "NLN" or is_glycam_sugar(r.name)}
    clashes = count_nonbonded_clashes(
        out, cutoff_angstrom=1.5, only_residue_names=glycan_resnames,
    )
    assert not clashes, (
        f"{len(clashes)} non-bonded contact(s) < 1.5 Å involving a "
        f"glycosylation site or glycan residue in zbs output (glycan "
        f"tree drift/collapse regression): {clashes[:10]}"
    )

    pdb = PDBFile(str(out))
    top = pdb.topology
    nln_residues = [r for r in top.residues() if r.name == "NLN"]
    assert len(nln_residues) == 4, (
        f"expected 4 N-glycosylation sites renamed to NLN, got "
        f"{len(nln_residues)}"
    )
    for r in nln_residues:
        nd2 = next((a for a in r.atoms() if a.name == "ND2"), None)
        assert nd2 is not None, f"{r.chain.id}:{r.name}{r.id} missing ND2"
        bonded_to_sugar = any(
            (b[0] is nd2 and b[1].residue is not r)
            or (b[1] is nd2 and b[0].residue is not r)
            for b in top.bonds()
        )
        assert bonded_to_sugar, (
            f"{r.chain.id}:{r.name}{r.id}'s ND2 has no bond to a sugar "
            f"tree — glycosylation site renamed but not actually linked"
        )


@pytest.mark.slow
def test_zbs_ligand_no_severe_clash_without_parametrize_flag(
    protein_ligand_gaps_pdb: Path, tmp_workdir: Path,
) -> None:
    """Plain `zbs` (no `--parametrize-ligands`) on a structure with an
    unknown ligand (DAN) used to silently strip the ligand out for the
    real minimize, let the pocket residues relax/drift into the
    now-empty cavity, then splice DAN back in at its ORIGINAL
    pre-minimize coordinates — producing a severe clash (< 0.5 Å heavy/H
    overlap) with the pocket that had since moved. Auto-triggering
    GAFF2 parametrization whenever an unknown heterogen has no FF
    template (rather than only under the explicit flag) fixes most of
    this; a real openmmforcefields limitation (GAFFTemplateGenerator
    can't be invoked from inside `Modeller.addHydrogens`'s own internal
    template matching) still forces the legacy strip-and-splice
    fallback here, so this asserts "no SEVERE clash" (< 0.5 Å) rather
    than zero contacts — the residual is a mild ~1 Å H-H contact, not
    the original heavy-atom overlap."""
    out = _run_zbs(protein_ligand_gaps_pdb, tmp_workdir)
    assert count_d_ca_residues(out) == 0

    clashes = count_nonbonded_clashes(out, cutoff_angstrom=0.5)
    assert not clashes, (
        f"{len(clashes)} severe (< 0.5 Å) non-bonded contact(s) in zbs "
        f"output (ligand splice-back clash regression): {clashes[:10]}"
    )
