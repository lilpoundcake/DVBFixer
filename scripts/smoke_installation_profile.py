#!/usr/bin/env python3
"""Run real scientific installation smokes and retain per-step logs/evidence.

Run with the NEW environment's Python, without PYTHONPATH. Outputs must be new.
This is acceptance tooling, not an alternative scientific implementation.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def chemistry(out: Path) -> None:
    import numpy as np
    import openmm
    from openbabel import openbabel as ob
    from openff.toolkit import Molecule
    from openmm.app import ForceField
    from openmmforcefields.generators import GAFFTemplateGenerator
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    Chem.MolToMolFile(mol, str(out / "ethanol.sdf"))
    conversion = ob.OBConversion()
    assert conversion.SetInAndOutFormats("sdf", "pdb")
    obmol = ob.OBMol()
    assert conversion.ReadFile(obmol, str(out / "ethanol.sdf"))
    assert obmol.NumAtoms() == 9
    assert conversion.WriteFile(obmol, str(out / "ethanol.pdb"))
    assert conversion.FindFormat("cif") is not None
    ligand = Molecule.from_file(str(out / "ethanol.sdf"))
    generator = GAFFTemplateGenerator(molecules=ligand, forcefield="gaff-2.11",
                                      cache=str(out / "gaff-cache.json"))
    ff = ForceField("amber14-all.xml")
    ff.registerTemplateGenerator(generator.generator)
    topology = ligand.to_topology().to_openmm()
    system = ff.createSystem(topology)
    assert system.getNumParticles() == 9
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(ligand.conformers[0].to_openmm())
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        openmm.unit.kilojoule_per_mole)
    assert np.isfinite(energy)
    print(f"RDKit/Open Babel SDF/PDB and OpenFF/GAFF system: {energy} kJ/mol")


def resp(out: Path) -> None:
    import numpy as np

    from dvbfixer.parametrize import _compute_resp_charges_pyscf

    path = out / "water.pdb"
    path.write_text("HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
                    "HETATM    2  H1  HOH A   1       0.758   0.000   0.504  1.00  0.00           H\n"
                    "HETATM    3  H2  HOH A   1      -0.758   0.000   0.504  1.00  0.00           H\nEND\n")
    charges = _compute_resp_charges_pyscf(str(path), verbose=True)
    assert len(charges) == 3 and np.isfinite(charges).all()
    assert abs(sum(charges)) < 1e-5 and charges[0] < 0
    print("RESP HF/6-31G* charges:", charges)


def antibody(out: Path) -> None:
    from anarci import anarci

    sequence = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGR"
                "FTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
    numbered, _, _ = anarci([("smoke", sequence)], scheme="imgt")
    assert numbered[0] and len(numbered[0][0][0]) > 100
    print("ANARCI + hmmscan + HMM data: IMGT V-domain numbered")


def analysis(out: Path, root: Path) -> None:
    import MDAnalysis as mda
    import numpy as np
    import plotly.graph_objects as go
    from MDAnalysis.analysis import align

    from dvbfixer.align import kabsch_align_pdb
    from dvbfixer.trajectory_input import load_universe

    source = root / "tests/fixtures/protein_context/1UBQ.pdb"
    reference = mda.Universe(str(source))
    with mda.Writer(str(out / "trajectory.dcd"), reference.atoms.n_atoms) as writer:
        writer.write(reference.atoms)
        reference.atoms.positions += np.array([1., 0., 0.])
        writer.write(reference.atoms)
    mobile = load_universe(source, out / "trajectory.dcd")
    assert len(mobile.trajectory) == 2
    _, rmsd = align.alignto(mobile, mda.Universe(str(source)), select="protein and name CA")
    assert rmsd < 1e-3
    with mda.Writer(str(out / "reference.gro"), n_atoms=reference.atoms.n_atoms) as writer:
        writer.write(reference.atoms)
    gro = load_universe(out / "reference.gro", out / "trajectory.dcd")
    assert gro.atoms.n_atoms == reference.atoms.n_atoms
    gro.trajectory.close()
    with mda.Writer(str(out / "trajectory.nc"), n_atoms=reference.atoms.n_atoms) as writer:
        writer.write(reference.atoms)
    netcdf = load_universe(source, out / "trajectory.nc")
    assert netcdf.atoms.n_atoms == reference.atoms.n_atoms
    assert np.allclose(netcdf.atoms.positions, reference.atoms.positions, atol=1e-4)
    netcdf.trajectory.close()
    _, fitted_rmsd, matched = kabsch_align_pdb(source, source, out / "fitted.pdb")
    assert matched > 50 and fitted_rmsd is not None and fitted_rmsd < 1e-3
    go.Figure(data=go.Scatter(y=[0., rmsd])).write_html(str(out / "plot.html"))
    assert (out / "plot.html").stat().st_size > 100
    print("MDAnalysis DCD/GRO/NetCDF, application structural fitting + Plotly HTML passed")


def formats(out: Path, root: Path) -> None:
    import gemmi

    from dvbfixer.structure_input import normalize_structure

    source = root / "tests/fixtures/protein_context/1UBQ.pdb"
    structure = gemmi.read_structure(str(source))
    structure.make_mmcif_document().write_file(str(out / "protein.cif"))
    normalized = out / "from-mmcif.pdb"
    result = normalize_structure(out / "protein.cif", normalized)
    assert result.dialect == "pdbx/mmCIF"
    def count(path):
        return sum(line.startswith(("ATOM  ", "HETATM")) for line in path.read_text().splitlines())
    assert count(normalized) == count(source)
    small = out / "small.cif"
    small.write_text("data_small\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n"
                     "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                     "_symmetry_space_group_name_H-M 'P 1'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
                     "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
                     "C1 C 0.1 0.2 0.3\nO1 O 0.22 0.2 0.3\n")
    result = normalize_structure(small, out / "from-crystal.pdb")
    assert result.dialect == "small-molecule CIF" and count(out / "from-crystal.pdb") == 2
    print("Gemmi PDBx/mmCIF and Open Babel crystallographic CIF conversion passed")


def glycan_trajectory(out: Path, root: Path) -> None:
    import MDAnalysis as mda

    universe = mda.Universe(str(root / "tests/fixtures/glycosilated_mAb_CHARMM.pdb"))
    with mda.Writer(str(out / "glycan.dcd"), n_atoms=universe.atoms.n_atoms) as writer:
        for _ in range(3):
            writer.write(universe.atoms)
            universe.atoms.positions += [0.1, 0., 0.]
    print("Three-frame glycoprotein analysis input generated")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--with-modeller", action="store_true")
    parser.add_argument("--extended", action="store_true", help="Formats, ligand/topology and real cluster/plotting smokes")
    parser.add_argument("--psi4-env", help="Also exercise full RESP parameterization through this separate QM environment")
    parser.add_argument("--msa-alternatives", action="store_true", help="Exercise MUSCLE 5 and Clustal Omega")
    parser.add_argument("--gmx", help="Optional GROMACS executable for topology preprocessing")
    parser.add_argument("--prepared-input", type=Path, help="Reuse an explicitly selected prepared PDB for targeted --only checks")
    parser.add_argument("--only", action="append", default=[], help="Run only named steps (repeatable)")
    parser.add_argument("--worker", choices=("chemistry", "resp", "antibody", "analysis", "formats", "glycan-trajectory"))
    args = parser.parse_args()
    out = args.output.resolve()
    root = args.root.resolve()
    if args.worker:
        if args.worker in {"analysis", "formats", "glycan-trajectory"}:
            {"analysis": analysis, "formats": formats, "glycan-trajectory": glycan_trajectory}[args.worker](out, root)
        else:
            {"chemistry": chemistry, "resp": resp, "antibody": antibody}[args.worker](out)
        return 0
    out.mkdir(parents=True, exist_ok=False)
    if args.prepared_input:
        shutil.copyfile(args.prepared_input, out / "legacy.pdb")
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("KEY_MODELLER", None)
    env.update(PYTHONNOUSERSITE="1", OMP_NUM_THREADS="2", OPENMM_CPU_THREADS="2")
    steps: list[tuple[str, list[str]]] = []
    for name in ("chemistry", "resp", "antibody", "analysis"):
        steps.append((name, [sys.executable, str(Path(__file__).resolve()), "--output", str(out),
                             "--root", str(root), "--worker", name]))
    source = str(root / "tests/fixtures/protein_context/1UBQ.pdb")
    cli = str(Path(sys.executable).parent / "dvbfixer")
    for backend in ("legacy", "tleap-reduce"):
        steps.append((f"prepare-{backend}", [cli, "prepare", source, "--backend", backend,
                      "-o", str(out / f"{backend}.pdb"), "--atom-naming", "standard", "--strip-heterogens"]))
    steps.extend([
        ("obminimize", [str(Path(sys.executable).parent / "obminimize"), "-ff", "UFF", "-n", "1",
                        str(out / "ethanol.pdb")]),
        ("protonate", [cli, "protonate", str(out / "legacy.pdb"), "-o", str(out / "protonated.pdb")]),
        ("diagnose", [cli, "diagnose", str(out / "legacy.pdb"), "--format", "json",
                      "-o", str(out / "diagnose.json")]),
        ("minimize", [cli, "minimize", str(out / "legacy.pdb"), "--no-solvent", "--max-iter", "10",
                      "--platform", "Reference", "-o", str(out / "minimized.pdb")]),
        ("topology", [cli, "top", str(out / "legacy.pdb"), "-o", str(out / "topol.top")]),
    ])
    fasta = out / "msa.fasta"
    fasta.write_text(">a\nACDEFGHIKLMNPQRSTVWY\n>b\nACDEFGHIKLMNPQRSAVWY\n")
    steps.append(("mafft", [cli, "msa", str(fasta), "-o", str(out / "aligned.fasta")]))
    if args.msa_alternatives:
        for engine in ("muscle", "clustalo"):
            steps.append((engine, [cli, "msa", str(fasta), "--engine", engine, "-o", str(out / f"{engine}.fasta")]))
    if args.extended:
        for worker in ("formats", "glycan-trajectory"):
            steps.append((worker, [sys.executable, str(Path(__file__).resolve()), "--worker", worker,
                                   "--root", str(root), "--output", str(out)]))
        steps.extend([
            ("cif-cli", [cli, "diagnose", str(out / "protein.cif"), "--format", "json", "-o", str(out / "cif-report.json")]),
            ("antibody-cli", [cli, "renumber", str(root / "tests/fixtures/regressions/1DQJ.pdb"),
                               "--scheme", "imgt", "-o", str(out / "antibody.pdb")]),
            ("cluster-plot", [cli, "cluster", str(root / "tests/fixtures/glycosilated_mAb_CHARMM.pdb"),
                              str(out / "glycan.dcd"), "--plot", "-o", str(out / "clusters")]),
        ])
        for method, flags in (("bcc", []), ("pyscf", ["-c", "resp", "--qm-engine", "pyscf", "--qm-nthreads", "1"])):
            directory = out / method
            directory.mkdir()
            steps.append((f"parametrize-{method}", [cli, "parametrize", str(out / "ethanol.sdf"), "-n", "LIG",
                          "-o", str(directory / "ligand"), *flags]))
        if args.psi4_env:
            directory = out / "psi4"
            directory.mkdir()
            steps.append(("parametrize-psi4", [cli, "parametrize", str(out / "ethanol.sdf"), "-n", "LIG",
                          "-o", str(directory / "ligand"), "-c", "resp", "--qm-engine", "psi4",
                          "--psi4-env", args.psi4_env, "--qm-nthreads", "1", "--qm-memory", "1GB"]))
        for name, ff, extra in (("acpype", "amber", ["--acpype"]),
                                ("custom-ff", "amber", ["--ff-dir", str(root / "FF/amber99sb-ildn-lipid21.ff")]),
                                ("amber19sb", "amber19sb", []),
                                ("charmm", "charmm", [])):
            directory = out / name
            directory.mkdir()
            steps.append((name, [cli, "top", str(out / "legacy.pdb"), "--ff", ff, *extra,
                                "-o", str(directory / "topol.top"), "--pdb", str(directory / "conf.pdb")]))
        if args.gmx:
            directory = out / "amber19sb"
            mdp = directory / "acceptance.mdp"
            mdp.write_text("integrator=md\nnsteps=0\ncontinuation=yes\npbc=xyz\ncutoff-scheme=Verlet\n"
                           "comm-mode=none\nverlet-buffer-tolerance=-1\nnstlist=10\nrlist=1.0\n"
                           "coulombtype=Cut-off\nrcoulomb=1.0\nvdwtype=Cut-off\nrvdw=1.0\n")
            gmx = shutil.which(args.gmx) or args.gmx
            steps.extend([
                ("gmx-box", [gmx, "editconf", "-f", str(directory / "conf.pdb"), "-o", str(directory / "boxed.gro"), "-box", "8", "8", "8"]),
                ("grompp", [gmx, "grompp", "-f", str(mdp), "-p", str(directory / "topol.top"),
                             "-c", str(directory / "boxed.gro"), "-o", str(directory / "acceptance.tpr"), "-maxwarn", "0"]),
            ])
    if args.with_modeller:
        from Bio import SeqIO

        sequence = str(next(SeqIO.parse(source, "pdb-seqres")).seq)
        target = out / "target.fasta"
        target.write_text(">A\n" + sequence + "\n")
        gap_input = out / "gap-input.pdb"
        gap_input.write_text("".join(line for line in Path(source).read_text().splitlines(keepends=True)
                                    if not line.startswith("HETATM") and not (
                                        line.startswith("ATOM  ") and int(line[22:26]) == 10)))
        steps.extend([
            ("model", [cli, "model", str(gap_input), "--fasta", str(target), "-o", str(out / "modeled.pdb"),
                       "--num-models", "1", "--num-loops", "1", "--md-level", "none", "--no-terminal"]),
            ("homology", [cli, "homology", str(target), "--template", source, "-o", str(out / "homology"),
                          "--num-models", "1", "--md-level", "none", "--no-loop-refine"]),
        ])
    results = []
    if args.only:
        unknown = set(args.only) - {name for name, _ in steps}
        if unknown:
            parser.error(f"Unknown/unselected steps: {sorted(unknown)}")
        steps = [(name, command) for name, command in steps if name in args.only]
    for name, command in steps:
        from dvbfixer.command_registry import get_command

        success_codes = get_command(command[1]).success_codes if command[0] == cli else (0,)
        with (out / f"{name}.log").open("w") as stream:
            try:
                workdir = out / "amber19sb" if name in {"gmx-box", "grompp"} else out
                process = subprocess.run(command, cwd=workdir, env=env, stdout=stream,
                                         stderr=subprocess.STDOUT, timeout=600, check=False)
                code = process.returncode
                if name in {"diagnose", "cif-cli"} and code in success_codes:
                    diagnostic = json.loads((out / ("diagnose.json" if name == "diagnose" else "cif-report.json")).read_text())
                    assert isinstance(diagnostic["findings"], list)
                if code == 0 and name.startswith("parametrize-"):
                    import parmed
                    directory = out / name.removeprefix("parametrize-")
                    assert (directory / "LIG.itp").stat().st_size > 100
                    assert len(parmed.load_file(str(directory / "LIG.gro")).atoms) == 9
                if code == 0 and name == "cluster-plot":
                    plots = list(out.glob("clusters*.html"))
                    assert plots and all(path.stat().st_size > 100 for path in plots)
                if code == 0 and name in {"prepare-legacy", "prepare-tleap-reduce", "minimize"}:
                    from openmm.app import PDBFile

                    from dvbfixer.ffutils.geometry import assert_all_l

                    pdb = PDBFile(str(out / ({"prepare-legacy": "legacy.pdb",
                                             "prepare-tleap-reduce": "tleap-reduce.pdb",
                                             "minimize": "minimized.pdb"}[name])))
                    assert pdb.topology.getNumAtoms() > 500
                    assert_all_l(pdb.topology, pdb.positions)
            except subprocess.TimeoutExpired:
                code = 124
        results.append({"name": name, "command": command, "returncode": code,
                        "success_codes": list(success_codes),
                        "state": "passed" if code in success_codes else "failed"})
        print(f"{name}: {results[-1]['state']}", flush=True)
    versions = {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()}
    conda = []
    for path in Path(sys.prefix, "conda-meta").glob("*.json"):
        metadata = json.loads(path.read_text())
        conda.append({k: metadata.get(k) for k in ("name", "version", "build", "channel")})
    report = {"prefix": sys.prefix, "python": sys.version, "versions": versions,
              "prepared_input": str(args.prepared_input.resolve()) if args.prepared_input else None,
              "conda": conda, "steps": results,
              "not_tested": ([] if args.with_modeller else ["MODELLER/Homology (license required)"])
                + ([] if args.extended else ["extended structural/format matrix"])
                + ["other managers/platforms", "external Gaussian execution", "diffusion hardware/checkpoints"]}
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return int(any(result["state"] != "passed" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
