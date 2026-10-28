"""Compare the historical Probe invocation, direct Probe and Python.

Optional input PDBs must already have agreed hydrogens. Synthetic methane
controls include all hydrogens. Direct Probe uses -self ALL and excludes water
to match the Python default; the historical broken command is retained as a control.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
from openmm import Vec3
from openmm.app import Element, PDBFile, Topology
from openmm.unit import Quantity, nanometer

from dvbfixer.diagnose import steric


def probe_identity(field):
    # Probe's 16-column identifier includes a two-character chain field.
    return (field[:2].strip(), field[2:6].strip(), field[6:7].strip(), field[11:15].strip())


def python_identity(atom):
    return tuple(atom[key] for key in ("chain", "resid", "icode", "atom"))


def build_directional_controls(output):
    """Complete methylamine/acetone molecules: same H-O distance, changed donor angle."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from scipy.spatial.transform import Rotation

    def molecule(smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(mol, randomSeed=20261028) != 0:
            raise RuntimeError("Directional-control embedding failed")
        return mol

    donor, acceptor = molecule("CN"), molecule("CC(=O)C")
    nitrogen = next(a for a in donor.GetAtoms() if a.GetSymbol() == "N")
    hydrogen = next(a for a in nitrogen.GetNeighbors() if a.GetSymbol() == "H")
    oxygen = next(a for a in acceptor.GetAtoms() if a.GetSymbol() == "O")
    carbon = oxygen.GetNeighbors()[0]

    def orient(mol, origin, axis, direction, target):
        coords = np.asarray(mol.GetConformer().GetPositions())
        coords = coords - coords[origin]
        vector = coords[axis] / np.linalg.norm(coords[axis])
        rotation, _ = Rotation.align_vectors(np.asarray([direction]), np.asarray([vector]))
        return rotation.apply(coords) + target

    donor_coords = orient(donor, nitrogen.GetIdx(), hydrogen.GetIdx(), [1., 0., 0.], np.zeros(3))
    h_position = donor_coords[hydrogen.GetIdx()]
    controls = {}
    for case, direction in [("linear-hbond", np.array([1., 0., 0.])),
                            ("bent-hbond", np.array([0., 1., 0.]))]:
        acceptor_coords = orient(acceptor, oxygen.GetIdx(), carbon.GetIdx(), direction,
                                 h_position + 1.8 * direction)
        top, coordinates = Topology(), []
        for chain_id, mol, coords in [("D", donor, donor_coords), ("A", acceptor, acceptor_coords)]:
            residue = top.addResidue("LIG", top.addChain(chain_id), "1")
            atoms = [top.addAtom(f"{a.GetSymbol()}{a.GetIdx()+1}", Element.getBySymbol(a.GetSymbol()), residue)
                     for a in mol.GetAtoms()]
            for bond in mol.GetBonds():
                top.addBond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
            coordinates.extend(Vec3(*point) for point in coords / 10)
        path = output / f"{case}.pdb"
        with path.open("w") as handle:
            PDBFile.writeFile(top, Quantity(coordinates, nanometer), handle, keepIds=True)
        controls[case] = {"path": path, "hydrogen_acceptor_distance_angstrom": 1.8,
                          "donor_hydrogen_acceptor_angle_degrees": 180 if case == "linear-hbond" else 90,
                          "target_pair": tuple(sorted((("D", "1", "", f"H{hydrogen.GetIdx()+1}"),
                                                       ("A", "1", "", f"O{oxygen.GetIdx()+1}")))),
                          "scope": "Directional target-pair control; other contacts are not assumed clean"}
    return controls


def build_bonded_controls(output):
    """Known covalent exclusion and known >=4-bond intramolecular collision."""
    from itertools import product

    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    controls = {}
    for case, smiles in [("bonded-ethane", "CC"), ("intramolecular-octane-clash", "CCCCCCCC")]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(mol, randomSeed=20261028) != 0:
            raise RuntimeError("Bonded-control embedding failed")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        target = (0, 1)
        if case == "intramolecular-octane-clash":
            candidates = [(i, j) for i in range(8) for j in range(i + 1, 8)
                          if len(Chem.GetShortestPath(mol, i, j)) >= 5]
            found = False
            nonbonds = [(i, j) for i in range(mol.GetNumAtoms()) for j in range(i+1, mol.GetNumAtoms())
                        if mol.GetBondBetweenAtoms(i, j) is None]
            for angles in product((0., 60., 120., 180., 240., 300.), repeat=4):
                for index, angle in enumerate(angles):
                    rdMolTransforms.SetDihedralDeg(mol.GetConformer(), index, index+1, index+2, index+3, angle)
                coords = np.asarray(mol.GetConformer().GetPositions())
                matches = [(i, j) for i, j in candidates if 2.4 < np.linalg.norm(coords[i] - coords[j]) < 2.9]
                no_extra_inferred_bonds = all(
                    np.linalg.norm(coords[i] - coords[j]) > steric._covalent_cutoff_a(
                        mol.GetAtomWithIdx(i).GetSymbol(), mol.GetAtomWithIdx(j).GetSymbol()
                    ) for i, j in nonbonds
                ) if matches else False
                if matches and no_extra_inferred_bonds:
                    target, found = matches[0], True
                    break
            if not found:
                raise RuntimeError("Cannot materialize the declared nonbonded collision")
        top, mapped = Topology(), {}
        coords = np.asarray(mol.GetConformer().GetPositions())
        ordered = []
        chain = top.addChain("A")
        residues = [top.addResidue("LIG", chain, "1")]
        if case == "bonded-ethane":
            residues.append(top.addResidue("LIG", chain, "2"))
        for residue_index, residue in enumerate(residues):
            for atom in mol.GetAtoms():
                parent = atom.GetIdx() if atom.GetSymbol() == "C" else atom.GetNeighbors()[0].GetIdx()
                if case == "bonded-ethane" and parent != residue_index:
                    continue
                mapped[atom.GetIdx()] = top.addAtom(f"{atom.GetSymbol()}{atom.GetIdx()+1}",
                                                   Element.getBySymbol(atom.GetSymbol()), residue)
                ordered.append(coords[atom.GetIdx()])
        for bond in mol.GetBonds():
            top.addBond(mapped[bond.GetBeginAtomIdx()], mapped[bond.GetEndAtomIdx()])
        path = output / f"{case}.pdb"
        with path.open("w") as handle:
            PDBFile.writeFile(top, Quantity([Vec3(*p) for p in np.asarray(ordered) / 10], nanometer), handle, keepIds=True)
        pair = tuple(sorted((a.residue.chain.id, a.residue.id, "", a.name) for a in (mapped[i] for i in target)))
        controls[case] = {"path": path, "target_pair": pair,
                          "target_expected_clash": case != "bonded-ethane",
                          "target_bond_graph_distance": len(Chem.GetShortestPath(mol, *target)) - 1,
                          "selection_policy": ("MMFF-relaxed covalent ethane" if case == "bonded-ethane" else
                                               "torsion grid, target 2.4-2.9 A, no accidental distance-inferred bonds; no detector findings used for selection"),
                          "target_distance_angstrom": float(np.linalg.norm(coords[target[0]] - coords[target[1]]))}
    return controls


def explain_discordant_pairs(pdb, python_pairs, probe_pairs, contacts):
    """Expose the production Python filters, without redefining their chemistry."""
    from scipy.spatial import cKDTree

    atoms = list(pdb.topology.atoms())
    coords = np.asarray(pdb.positions.value_in_unit(nanometer)) * 10
    nearby = cKDTree(coords).query_pairs(r=steric._NEIGHBOR_CUTOFF_A)
    neighbors = steric._build_neighbor_graph(pdb.topology, atoms, coords, nearby)
    excluded = steric._expand_exclusions(neighbors)
    lookup = {}
    for atom in atoms:
        identity = (atom.residue.chain.id, atom.residue.id, atom.residue.insertionCode.strip(), atom.name)
        lookup.setdefault(identity, []).append(atom)
    result = []
    for pair in sorted(python_pairs ^ probe_pairs):
        row = {"atoms": pair, "reported_by": "python" if pair in python_pairs else "probe"}
        matching = [c for c in contacts if tuple(sorted((probe_identity(c["source"]),
                                                       probe_identity(c["target"])))) == pair]
        row["probe_contact_types"] = sorted({c["type"] for c in matching})
        row["probe_max_overlap_angstrom"] = max((c["overlap_angstrom"] for c in matching), default=None)
        selected = [lookup.get(identity, []) for identity in pair]
        if any(len(matches) != 1 for matches in selected):
            row["identity_status"] = "missing or ambiguous; not adjudicated"
        else:
            a, b = (matches[0] for matches in selected)
            distance = float(np.linalg.norm(coords[a.index] - coords[b.index]))
            filters = []
            if frozenset((a.index, b.index)) in excluded:
                filters.append("bond_graph_1_2_1_3_1_4")
            if a.residue is b.residue:
                filters.append("same_residue")
            if steric._is_hbond_pair(a, b, distance, atoms, neighbors):
                filters.append("polar_hydrogen_envelope")
            if steric._is_heavy_hbond_pair(a, b, distance):
                filters.append("heavy_hbond_envelope")
            row.update({"identity_status": "resolved", "distance_angstrom": distance,
                        "elements": [a.element.symbol, b.element.symbol],
                        "residue_names": [a.residue.name, b.residue.name],
                        "python_overlap_angstrom": sum(steric._VDW_A.get(x.element.symbol.upper(), 1.7)
                                                       for x in (a, b)) - distance,
                        "python_exclusion_filters": filters})
        result.append(row)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input", type=Path, action="append", default=[])
    parser.add_argument("--align-atom-fields", action="store_true",
                        help="Create formatting-only copies using PDB element-aware atom-field alignment")
    parser.add_argument("--probe-nuclear", action="store_true",
                        help="Sensitivity arm: use Probe nuclear-position radii in the direct invocation only")
    args = parser.parse_args()
    binary = str(args.probe.resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "audit-script.py").write_bytes(Path(__file__).read_bytes())
    version = subprocess.run([binary, "-version"], capture_output=True, text=True, timeout=30)
    (args.output / "version.txt").write_text(version.stdout + version.stderr)
    inputs = []
    for case, separation in [("clean-methanes", 7.0), ("clashing-methanes", 2.5)]:
        top = Topology()
        coords = []
        for chain_id, shift in [("A", 0), ("a", separation)]:
            residue = top.addResidue("LIG", top.addChain(chain_id), "1")
            carbon = top.addAtom("C", Element.getBySymbol("C"), residue)
            coords.append(Vec3(shift / 10, 0, 0))
            for i, signs in enumerate([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]):
                hydrogen = top.addAtom(f"H{i+1}", Element.getBySymbol("H"), residue)
                top.addBond(carbon, hydrogen)
                coords.append(Vec3((shift + signs[0] * 0.6293) / 10,
                                   signs[1] * 0.6293 / 10, signs[2] * 0.6293 / 10))
        path = args.output / f"{case}.pdb"
        with path.open("w") as handle:
            PDBFile.writeFile(top, Quantity(coords, nanometer), handle, keepIds=True)
        inputs.append(path)
    inputs.extend(args.input)
    directional_controls = build_directional_controls(args.output)
    inputs.extend(control["path"] for control in directional_controls.values())
    bonded_controls = build_bonded_controls(args.output)
    inputs.extend(control["path"] for control in bonded_controls.values())
    report = {"platform": platform.platform(), "binary": binary,
               "binary_sha256": hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "invocation": sys.argv,
               "probe_direct_radii": "nuclear" if args.probe_nuclear else "electron-cloud (Probe default)",
              "warning_overlap_angstrom": 0.4, "error_overlap_angstrom": 0.5,
              "scope": "contact comparison, not a validated false-positive benchmark", "runs": []}
    for index, path in enumerate(inputs):
        work = args.output / f"{index}-{path.stem}"
        work.mkdir()
        original = path
        if args.align_atom_fields:
            lines = []
            for line in path.read_text().splitlines(keepends=True):
                if line.startswith(("ATOM  ", "HETATM")):
                    name, element = line[12:16].strip(), line[76:78].strip()
                    if len(element) == 1 and len(name) < 4:
                        line = line[:12] + f" {name:<3}" + line[16:]
                lines.append(line)
            path = work / "aligned.pdb"
            path.write_text("".join(lines))
        pdb = PDBFile(str(path))
        record = {"input": str(path.resolve()), "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "source_input": str(original.resolve()),
                  "source_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
                  "atom_fields_aligned": args.align_atom_fields,
                  "atoms": pdb.topology.getNumAtoms(),
                  "hydrogens": sum(a.element.symbol == "H" for a in pdb.topology.atoms())}
        start = time.monotonic()
        findings = steric.clashes_python(pdb.topology, pdb.positions, clash_warn_a=0.4, clash_error_a=0.5)
        record["python_seconds"] = time.monotonic() - start
        record["python_findings"] = [asdict(f) for f in findings]
        for mode, selection in [("historical-production", ["ALL", "ALL"]), ("direct-self", ["-nowaters", "-self", "ALL"])]:
            cmd = [binary, "-unformated", "-condensed", "-Q", "-mc", "-het", *selection, str(path.resolve())]
            if mode == "direct-self" and args.probe_nuclear:
                cmd.insert(1, "-nuclear")
            start = time.monotonic()
            run = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            seconds = time.monotonic() - start
            (work / f"{mode}.stdout").write_text(run.stdout)
            (work / f"{mode}.stderr").write_text(run.stderr)
            contacts = []
            for line in run.stdout.splitlines():
                fields = line.split(":")
                if len(fields) >= 20:
                    contacts.append({"type": fields[2], "source": fields[3], "target": fields[4],
                                     "dot_count": int(fields[5]), "overlap_angstrom": -float(fields[6])})
            record[mode] = {"command": cmd, "returncode": run.returncode, "seconds": seconds,
                            "contacts": contacts,
                            "unique_overlap_pairs_ge_0.4": len({
                                tuple(sorted((c["source"], c["target"]))) for c in contacts
                                if c["type"] in ("bo", "so") and c["overlap_angstrom"] >= 0.4
                            }),
                            "bad_overlap_contacts_ge_0.4": sum(c["type"] in ("bo", "so") and
                                                               c["overlap_angstrom"] >= 0.4 for c in contacts)}
        report["runs"].append(record)
        python_pairs = {tuple(sorted(python_identity(atom) for atom in f.extra["atoms"]))
                        for f in findings}
        probe_pairs = {tuple(sorted((probe_identity(c["source"]), probe_identity(c["target"]))))
                       for c in record["direct-self"]["contacts"]
                       if c["type"] in ("bo", "so") and c["overlap_angstrom"] >= 0.4}
        record["pair_comparison"] = {
            "shared": sorted(python_pairs & probe_pairs),
            "python_only": sorted(python_pairs - probe_pairs),
            "probe_only": sorted(probe_pairs - python_pairs),
            "interpretation": "Different radii/exclusions; discordance alone is not a quality verdict",
        }
        record["discordant_pair_details"] = explain_discordant_pairs(
            pdb, python_pairs, probe_pairs, record["direct-self"]["contacts"]
        )
        if original.stem in directional_controls:
            control = directional_controls[original.stem]
            record["directional_control"] = {k: v for k, v in control.items() if k != "path"}
            target = control["target_pair"]
            record["directional_control"]["python_reports_target_pair"] = target in python_pairs
            record["directional_control"]["probe_target_contact_types"] = sorted({
                c["type"] for c in record["direct-self"]["contacts"]
                if tuple(sorted((probe_identity(c["source"]), probe_identity(c["target"])))) == target
            })
        if original.stem in bonded_controls:
            control = bonded_controls[original.stem]
            record["bonded_control"] = {k: v for k, v in control.items() if k != "path"}
            record["bonded_control"]["python_reports_target_pair"] = control["target_pair"] in python_pairs
            record["bonded_control"]["probe_reports_target_pair"] = control["target_pair"] in probe_pairs
        (args.output / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
        print(path.name, "python", len(findings), "historical production rc", record["historical-production"]["returncode"],
              "direct Probe", record["direct-self"]["bad_overlap_contacts_ge_0.4"], flush=True)
    manifest = {str(path.relative_to(args.output)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(args.output.rglob("*")) if path.is_file()}
    (args.output / "manifest.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
