"""Shared force field utilities for dvbfixer.

Provides OpenFF-based parametrization for non-standard residues (glycans,
ligands) that lack templates in standard AMBER force fields. Uses
SMIRNOFFTemplateGenerator from openmmforcefields to auto-parametrize
unknown residues on the fly.
"""

from openmm.app import ForceField

# SMILES for common glycan residues (from PDB Chemical Component Dictionary).
# These are used to create OpenFF Molecule objects for parametrization.
KNOWN_GLYCAN_SMILES = {
    'NAG': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O',   # N-acetyl-D-glucosamine
    'NDG': 'CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@@H]1O',   # N-acetyl-D-glucosamine (alt)
    'BMA': 'OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',              # beta-D-mannose
    'MAN': 'OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',              # alpha-D-mannose
    'FUC': 'C[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',               # L-fucose
    'FUL': 'C[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O',               # L-fucose (alt)
    'GAL': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # D-galactose
    'BGC': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # beta-D-glucose
    'GLC': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',               # D-glucose
    'SIA': 'OC[C@@H](O)[C@@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1NC(C)=O',  # sialic acid (simplified)
}

# Standard protein/water/ion residues that don't need OpenFF parametrization
PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'CYX', 'CYM', 'ASH', 'GLH',
    'LYN', 'MSE', 'ACE', 'NME', 'NHE',
}
SOLVENT_IONS = {
    'HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN',
    'NA+', 'CL-', 'K+', 'MG2+', 'CA2+',
}


def _find_unknown_residue_names(topology):
    """Find residue names in topology that aren't protein/water/ions."""
    known = PROTEIN_RESIDUES | SOLVENT_IONS
    unknown = set()
    for res in topology.residues():
        if res.name not in known:
            unknown.add(res.name)
    return unknown


def create_forcefield_with_openff(ff_xmls, topology, small_mol_ff='openff-2.2.0',
                                  extra_molecules=None, verbose=False):
    """Create OpenMM ForceField with automatic OpenFF parametrization for unknown residues.

    Args:
        ff_xmls: List of force field XML files (e.g., ['amber19/protein.ff19SB.xml', ...])
        topology: OpenMM Topology to scan for unknown residues
        small_mol_ff: OpenFF force field name (default: openff-2.2.0 Sage)
        extra_molecules: Optional list of additional openff.toolkit.Molecule objects
        verbose: Print info about parametrized residues

    Returns:
        ForceField object with SMIRNOFFTemplateGenerator registered for unknown residues
    """
    ff = ForceField(*ff_xmls)

    unknown = _find_unknown_residue_names(topology)
    if not unknown and not extra_molecules:
        return ff

    # Build Molecule objects for known glycans
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    molecules = []
    matched = set()
    for resname in unknown:
        if resname in KNOWN_GLYCAN_SMILES:
            mol = Molecule.from_smiles(KNOWN_GLYCAN_SMILES[resname],
                                       allow_undefined_stereo=True)
            molecules.append(mol)
            matched.add(resname)

    if extra_molecules:
        molecules.extend(extra_molecules)

    unmatched = unknown - matched
    if unmatched and verbose:
        print(f"Warning: no SMILES for residues: {', '.join(sorted(unmatched))}")
        print("  These residues may cause errors. Use --sdf to provide molecule definitions.")

    if molecules:
        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield=small_mol_ff)
        ff.registerTemplateGenerator(smirnoff.generator)
        if verbose:
            print(f"OpenFF parametrization registered for: {', '.join(sorted(matched))}")

    return ff
