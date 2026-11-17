#!/usr/bin/env python3
"""Read exact residual atoms/rings from saved endpoints, not a causal replay."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.diagnose.chemistry import _dihedral_deg, check_bond_lengths
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity
from dvbfixer.model.diffusion.covalent_ring_probe import probe_covalent_ring
from dvbfixer.model.diffusion.gap_portfolio import _coordinates


def residue_cycles(atoms, bonds):
    neighbors = {atom: set() for atom in atoms}
    for edge in bonds:
        if edge <= atoms:
            a, b = tuple(edge)
            neighbors[a].add(b)
            neighbors[b].add(a)
    found = set()

    def visit(path):
        for atom in sorted(neighbors[path[-1]]):
            if atom == path[0] and len(path) >= 3:
                rotations = [tuple(path[i:]+path[:i]) for i in range(len(path))]
                reverse = list(reversed(path))
                rotations += [tuple(reverse[i:]+reverse[:i]) for i in range(len(path))]
                found.add(min(rotations))
            elif atom not in path and len(path) < 8:
                visit([*path, atom])

    for atom in sorted(atoms):
        visit([atom])
    return sorted(found)


def analyze(root, audit_path, output):
    audit = json.loads(audit_path.read_text())
    for name, expected in audit["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("audited repair evidence changed")
    rows, bindings = [], {str(audit_path): digest(audit_path)}
    for row in audit["rows"]:
        case = row["case"]
        work = root / "execution" / f"{case}-seed-11"
        request = DiffusionRequest.from_json((work / "request.json").read_text())
        generated = frozenset(request.generated_atoms)
        paths = [("raw", work / "raw.pdb")]+[(stage["stage"], Path(stage["final_path"])) for stage in row["stages"]]
        for stage, path in paths:
            text = path.read_text()
            points, _ = _coordinates(text)
            structure = _structure(text, "residual.pdb")
            ids = {a.index: AtomIdentity(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name)
                   for a in structure.pdb.topology.atoms()}
            bonds = frozenset(frozenset((ids[a.index], ids[b.index])) for a, b in structure.pdb.topology.bonds())
            residual_bonds = []
            for finding in check_bond_lengths(structure.pdb.topology, structure.pdb.positions):
                if finding.category == "bond_length" and finding.extra:
                    pair = [AtomIdentity(a["chain"], a["resid"], a["icode"], a["atom"]) for a in finding.extra["atoms"]]
                    if any(atom in generated for atom in pair):
                        residual_bonds.append({"atom_identities": [asdict(a) for a in pair], "diagnose_geometry": finding.extra})
            phipsi, rings = [], []
            if case == "1NON":
                for chain in ("A", "D"):
                    residue = ResidueIdentity(chain, "74", "")
                    def xyz(number, name, chain=chain, coordinates=points):
                        return tuple(coordinates[AtomIdentity(chain, str(number), "", name)])
                    phipsi.append({"residue": asdict(residue),
                        "phi_degrees": _dihedral_deg(xyz(73, "C"), xyz(74, "N"), xyz(74, "CA"), xyz(74, "C")),
                        "psi_degrees": _dihedral_deg(xyz(74, "N"), xyz(74, "CA"), xyz(74, "C"), xyz(75, "N"))})
            elif case == "2BE1":
                focus = {a for a in points if a.chain == "B" and a.residue_number == "256" and a.insertion_code == ""}
                for cycle in residue_cycles(focus, bonds):
                    rings.append(asdict(probe_covalent_ring(cycle, points, bonds, generated)))
            ca_cb = None
            if case == "2BE1":
                a, b = AtomIdentity("B", "256", "", "CA"), AtomIdentity("B", "256", "", "CB")
                ca_cb = float(np.linalg.norm(points[a]-points[b]))
            rows.append({"case": case, "stage": stage, "generated_bond_findings": residual_bonds,
                "selected_phi_psi": phipsi, "B256_CA_CB_angstrom": ca_cb, "B256_covalent_ring_probes": rings,
                "topology_basis": "OpenMM endpoint PDB standard-residue graph, not unsaved temporary refiner topology",
                "optimization_trajectory_or_cause_established": False, "new_optimization_calls": 0})
            bindings[str(path)] = digest(path)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "residual-readback.json", {"rows": rows, "bindings": bindings,
        "analyst_sha256": digest(Path(__file__)), "new_sampling_or_repair": False, "fatal_ring_gate": False})
    print(json.dumps({"saved_endpoints_read": len(rows), "ring_piercings": sum(len(probe["piercings"]) for row in rows for probe in row["B256_covalent_ring_probes"])}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.root.resolve(), args.audit.resolve(), args.output.resolve())
