"""On-the-fly GAFF2 + AM1-BCC parametrisation of unknown ligands for OpenMM.

Used by `dvbfixer minimize --parametrize-ligands` (and `zbs` when the flag
is forwarded). For each residue in the input PDB whose name has no template
in the resolved base FF (and isn't a standard protein / water / ion), we:

  1. Extract that residue's atoms + coords via OpenBabel → SDF (has bond
     orders and 3D stereochemistry).
  2. Wrap in an `openff.toolkit.Molecule`.
  3. Feed the whole list of molecules to `openmmforcefields.generators.
     GAFFTemplateGenerator` — a proper AMBER GAFF template generator that
     runs `antechamber` + `parmchk2` under the hood and caches to disk.
  4. Return the generator so callers can `ff.registerTemplateGenerator(gen)`
     before `ForceField.createSystem(...)`.

Cached in `~/.cache/dvbfixer/lig_params/` (or `$DVBFIXER_LIG_CACHE`) via
GAFFTemplateGenerator's built-in JSON cache keyed by SMILES.

**Not a replacement for SMIRNOFF** — GAFFTemplateGenerator uses the AMBER
GAFF2 force field with AM1-BCC charges (via antechamber), which is
battle-tested for arbitrary organic ligands. Same limitation as SMIRNOFF /
GLYCAM re. cross-residue bonds: bonds between two ligand residues get no
parameters. Use only for isolated ligands.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path


def _default_cache_dir():
    override = os.environ.get('DVBFIXER_LIG_CACHE')
    if override:
        return Path(override)
    xdg = os.environ.get('XDG_CACHE_HOME')
    base = Path(xdg) if xdg else Path.home() / '.cache'
    return base / 'dvbfixer' / 'lig_params'


def _extract_residue_sdf(pdb_path, chain_id, res_id, out_sdf, verbose=False):
    """Write a single-residue SDF via OpenBabel.

    OpenBabel perceives bond orders + hybridisation from 3D coords, which
    is what antechamber then consumes. Returns True on success.
    """
    try:
        from openbabel import openbabel as ob
    except ImportError:
        if verbose:
            print("  [lig_params] openbabel not available — cannot extract ligand")
        return False

    ob.obErrorLog.SetOutputLevel(0)
    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats('pdb', 'sdf')
    mol = ob.OBMol()
    if not obconv.ReadFile(mol, str(pdb_path)):
        return False

    # Build a sub-OBMol containing only the target residue's atoms.
    sub = ob.OBMol()
    keep = []
    for atom in ob.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res is None:
            continue
        if res.GetChain() != chain_id:
            continue
        if res.GetNum() != res_id:
            continue
        keep.append(atom.GetIdx())

    if not keep:
        return False

    old_to_new = {}
    for old_idx in keep:
        atom = mol.GetAtom(old_idx)
        new_atom = sub.NewAtom()
        new_atom.SetAtomicNum(atom.GetAtomicNum())
        new_atom.SetVector(atom.x(), atom.y(), atom.z())
        old_to_new[old_idx] = new_atom.GetIdx()

    for bond in ob.OBMolBondIter(mol):
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in old_to_new and b in old_to_new:
            sub.AddBond(old_to_new[a], old_to_new[b], bond.GetBondOrder())

    return obconv.WriteFile(sub, str(out_sdf))


def _iter_unique_ligand_residues(topology, base_ff_templates):
    """Yield (chain_id, res_id, resname) for residues that need parametrisation.

    "Need parametrisation" = the resname has no template in the base
    ForceField AND it's not protein / water / ion / already-GLYCAM / already-
    CHARMM (those are handled by the base FF).
    """
    from dvbfixer.ffutils import (
        PROTEIN_RESIDUES, SOLVENT_IONS, GLYCAM_PROTEIN_RESIDUES,
        GLYCAM_CAPS, is_glycam_sugar,
    )

    skip = PROTEIN_RESIDUES | SOLVENT_IONS | GLYCAM_PROTEIN_RESIDUES | GLYCAM_CAPS
    seen = set()
    for res in topology.residues():
        if res.name in skip:
            continue
        if is_glycam_sugar(res.name):
            continue
        if res.name in base_ff_templates:
            continue
        key = res.name
        if key in seen:
            continue
        seen.add(key)
        yield (res.chain.id, res.id, res.name)


def build_ligand_generator(pdb_path, topology, base_ff_templates, *,
                           cache_dir=None, charge_method='bcc',
                           verbose=False):
    """Return an OpenMM template generator for every unknown ligand in `topology`.

    Args:
        pdb_path: Path to the source PDB (used to extract ligand geometry).
        topology: OpenMM Topology scanned for unknown-to-FF residues.
        base_ff_templates: Set of resnames the base FF already knows
            (typically `set(ForceField(*ff_xmls)._templates)`).
        cache_dir: Directory for GAFFTemplateGenerator's on-disk cache
            (default: `~/.cache/dvbfixer/lig_params/`).
        charge_method: `'bcc'` (default; AM1-BCC) or `'gasteiger'` (much
            faster, lower quality). Passed straight to antechamber via
            GAFFTemplateGenerator.
        verbose: Print per-ligand progress.

    Returns:
        Registered `GAFFTemplateGenerator` instance, or None if no unknown
        ligands were found or GAFF isn't available.
    """
    unknown = list(_iter_unique_ligand_residues(topology, base_ff_templates))
    if not unknown:
        if verbose:
            print("  [lig_params] no unknown-to-FF ligand residues")
        return None

    try:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import GAFFTemplateGenerator
    except ImportError as e:
        print(f"WARNING: --parametrize-ligands requires openmmforcefields "
              f"and openff-toolkit ({e}); skipping ligand parametrisation")
        return None

    cache_dir = Path(cache_dir) if cache_dir else _default_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / 'gaff_ligands.json'

    molecules = []
    with tempfile.TemporaryDirectory(prefix='dvbfixer_lig_') as tmp:
        for chain_id, res_id, resname in unknown:
            sdf_path = Path(tmp) / f'{resname}_{chain_id}_{res_id}.sdf'
            if not _extract_residue_sdf(pdb_path, chain_id, res_id,
                                         sdf_path, verbose=verbose):
                print(f"WARNING: could not extract ligand {resname} "
                      f"(chain {chain_id}, res {res_id}) via OpenBabel; "
                      f"skipping")
                continue
            try:
                mol = Molecule.from_file(str(sdf_path),
                                         allow_undefined_stereo=True)
                mol.name = resname  # so GAFF template gets registered by resname
                molecules.append(mol)
                if verbose:
                    print(f"  [lig_params] extracted {resname} "
                          f"({mol.n_atoms} atoms)")
            except Exception as e:
                print(f"WARNING: openff.toolkit could not read "
                      f"{resname}.sdf ({e}); skipping")

        if not molecules:
            return None

        try:
            gen = GAFFTemplateGenerator(molecules=molecules,
                                        forcefield='gaff-2.11',
                                        cache=str(cache_file))
        except Exception as e:
            print(f"WARNING: GAFFTemplateGenerator failed ({e}); "
                  f"skipping ligand parametrisation")
            return None

    if verbose:
        print(f"  [lig_params] built GAFF2 templates for: "
              f"{', '.join(m.name for m in molecules)} "
              f"(cache: {cache_file})")
    return gen
