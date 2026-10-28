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

    source = root / "tests/fixtures/protein_context/1UBQ.pdb"
    reference = mda.Universe(str(source))
    with mda.Writer(str(out / "trajectory.dcd"), reference.atoms.n_atoms) as writer:
        writer.write(reference.atoms)
        reference.atoms.positions += np.array([1., 0., 0.])
        writer.write(reference.atoms)
    mobile = mda.Universe(str(source), str(out / "trajectory.dcd"))
    assert len(mobile.trajectory) == 2
    _, rmsd = align.alignto(mobile, mda.Universe(str(source)), select="protein and name CA")
    assert rmsd < 1e-3
    go.Figure(data=go.Scatter(y=[0., rmsd])).write_html(str(out / "plot.html"))
    assert (out / "plot.html").stat().st_size > 100
    print("MDAnalysis trajectory/alignment + Plotly HTML passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--with-modeller", action="store_true")
    parser.add_argument("--only", action="append", default=[], help="Run only named steps (repeatable)")
    parser.add_argument("--worker", choices=("chemistry", "resp", "antibody", "analysis"))
    args = parser.parse_args()
    out = args.output.resolve()
    root = args.root.resolve()
    if args.worker:
        if args.worker == "analysis":
            analysis(out, root)
        else:
            {"chemistry": chemistry, "resp": resp, "antibody": antibody}[args.worker](out)
        return 0
    out.mkdir(parents=True, exist_ok=False)
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
        ("diagnose", [cli, "diagnose", str(out / "legacy.pdb")]),
        ("minimize", [cli, "minimize", str(out / "legacy.pdb"), "--no-solvent", "--max-iter", "10",
                      "--platform", "Reference", "-o", str(out / "minimized.pdb")]),
        ("topology", [cli, "top", str(out / "legacy.pdb"), "-o", str(out / "topol.top")]),
    ])
    fasta = out / "msa.fasta"
    fasta.write_text(">a\nACDEFGHIKLMNPQRSTVWY\n>b\nACDEFGHIKLMNPQRSAVWY\n")
    steps.append(("mafft", [cli, "msa", str(fasta), "-o", str(out / "aligned.fasta")]))
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
        with (out / f"{name}.log").open("w") as stream:
            try:
                process = subprocess.run(command, cwd=out, env=env, stdout=stream,
                                         stderr=subprocess.STDOUT, timeout=600, check=False)
                code = process.returncode
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
                        "state": "passed" if code == 0 else "failed"})
        print(f"{name}: {results[-1]['state']}", flush=True)
    versions = {dist.metadata["Name"]: dist.version for dist in importlib.metadata.distributions()}
    conda = []
    for path in Path(sys.prefix, "conda-meta").glob("*.json"):
        metadata = json.loads(path.read_text())
        conda.append({k: metadata.get(k) for k in ("name", "version", "build", "channel")})
    report = {"prefix": sys.prefix, "python": sys.version, "versions": versions,
              "conda": conda, "steps": results,
              "not_tested": ([] if args.with_modeller else ["MODELLER/Homology (license required)"]) + ["full WP5 structural/format matrix",
                             "other managers/platforms", "diffusion hardware/checkpoints"]}
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return int(any(result["returncode"] for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
