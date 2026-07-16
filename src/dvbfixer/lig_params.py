"""On-the-fly GAFF2 + AM1-BCC parametrisation of unknown ligands for OpenMM.

Used by `dvbfixer minimize --parametrize-ligands` (and `zbs` when the flag
is forwarded). For each residue in the input topology whose name has no
template in the resolved base FF (and isn't a standard protein / water /
ion), we:

  1. Dump the topology + positions to a temp PDB (guarantees byte-level
     agreement between what we iterate and what OpenBabel sees).
  2. Extract that residue's atoms + coords via OpenBabel → SDF (has bond
     orders and 3D stereochemistry).
  3. Wrap in an `openff.toolkit.Molecule`.
  4. Feed the whole list of molecules to `openmmforcefields.generators.
     GAFFTemplateGenerator` — a proper AMBER GAFF template generator that
     runs `antechamber` + `parmchk2` under the hood and caches to disk.
  5. Return the generator so callers can `ff.registerTemplateGenerator(gen)`
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

import atexit
import os
import tempfile
from pathlib import Path


class LigandParamError(RuntimeError):
    """Raised by `build_ligand_generator(strict=True)` when parametrisation
    of any unknown ligand fails. Message includes the specific ligand and
    a hint. Not caught by the wider minimize createSystem fallback."""


def _default_cache_dir():
    override = os.environ.get('DVBFIXER_LIG_CACHE')
    if override:
        return Path(override)
    xdg = os.environ.get('XDG_CACHE_HOME')
    base = Path(xdg) if xdg else Path.home() / '.cache'
    return base / 'dvbfixer' / 'lig_params'


def _dump_topology_to_pdb(topology, positions):
    """Write `topology` + `positions` to a temp PDB and register cleanup.

    Returns the path. Ensures OpenBabel sees exactly the residues we're
    iterating in `_iter_unique_ligand_residues` — no drift from post-CONECT
    or post-rename processing between the topology and the caller's args.
    """
    from openmm.app import PDBFile

    fd, tmp_path = tempfile.mkstemp(prefix='dvbfixer_ligsrc_', suffix='.pdb')
    os.close(fd)
    with open(tmp_path, 'w') as f:
        PDBFile.writeFile(topology, positions, f, keepIds=True)

    def _cleanup(p=tmp_path):
        try:
            os.unlink(p)
        except OSError:
            pass

    atexit.register(_cleanup)
    return tmp_path


def _resid_to_int(res_id):
    """Best-effort int cast for an OpenMM residue.id string.

    OpenMM stores `Residue.id` as `str`. `OBResidue.GetNum()` returns
    `int`. Strip a trailing insertion code if present; return None on
    total failure so the caller can log + skip cleanly.
    """
    if isinstance(res_id, int):
        return res_id
    s = str(res_id).strip()
    if not s:
        return None
    # Trim insertion-code suffix (e.g. '100A' → '100')
    end = len(s)
    while end > 0 and not s[end - 1].isdigit():
        end -= 1
    if end == 0:
        return None
    try:
        return int(s[:end])
    except ValueError:
        return None


def _extract_residue_sdf(pdb_path, chain_id, res_id, out_sdf, verbose=False):
    """Write a single-residue SDF via OpenBabel.

    OpenBabel perceives bond orders + hybridisation from 3D coords, which
    is what antechamber then consumes. Returns True on success, False
    otherwise (with a printed reason).
    """
    try:
        from openbabel import openbabel as ob
    except ImportError:
        print("  [lig_params] openbabel not available — cannot extract ligand")
        return False

    target_num = _resid_to_int(res_id)
    if target_num is None:
        print(f"  [lig_params] non-integer residue id {res_id!r} — cannot "
              f"extract (insertion codes on ligands are not supported)")
        return False
    target_chain = str(chain_id or ' ')

    ob.obErrorLog.SetOutputLevel(0)
    obconv = ob.OBConversion()
    obconv.SetInAndOutFormats('pdb', 'sdf')
    mol = ob.OBMol()
    if not obconv.ReadFile(mol, str(pdb_path)):
        print(f"  [lig_params] OpenBabel could not read {pdb_path}")
        return False

    # Build a sub-OBMol containing only the target residue's atoms.
    sub = ob.OBMol()
    keep = []
    for atom in ob.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res is None:
            continue
        # OB chain is a single char; empty chain in the PDB shows up as ' '
        obchain = res.GetChain() or ' '
        if obchain != target_chain:
            continue
        if res.GetNum() != target_num:
            continue
        keep.append(atom.GetIdx())

    if not keep:
        print(f"  [lig_params] no atoms found in {pdb_path} for chain "
              f"{target_chain!r} resnum {target_num}")
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

    if not obconv.WriteFile(sub, str(out_sdf)):
        print(f"  [lig_params] OpenBabel could not write {out_sdf}")
        return False
    return True


def _iter_unique_ligand_residues(topology, base_ff_templates):
    """Yield (chain_id, res_id, resname) for residues that need parametrisation.

    "Need parametrisation" = the resname has no template in the base
    ForceField AND it's not protein / water / ion / already-GLYCAM (those
    are handled by the base FF).
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


def build_ligand_generator(topology, positions, base_ff_templates, *,
                           strict=True, cache_dir=None, verbose=False):
    """Return an OpenMM template generator for every unknown ligand in `topology`.

    Args:
        topology: OpenMM Topology scanned for unknown-to-FF residues.
        positions: OpenMM positions matching `topology` (Quantity or list
            of Vec3). Dumped alongside the topology to a temp PDB for
            OpenBabel to read.
        base_ff_templates: Set of resnames the base FF already knows
            (typically `set(ForceField(*ff_xmls)._templates)`).
        strict: When True (default), raise `LigandParamError` if any
            candidate ligand fails to extract or generate. When False,
            fall through to a warning + return None (used e.g. from
            best-effort auto paths).
        cache_dir: Directory for GAFFTemplateGenerator's on-disk cache
            (default: `~/.cache/dvbfixer/lig_params/`).
        verbose: Print per-ligand progress.

    Returns:
        Registered `GAFFTemplateGenerator` instance, or None if there are
        no unknown ligands to parametrise (empty case is NOT an error).

    Raises:
        LigandParamError: (strict=True only) if any ligand fails to
            extract or if GAFFTemplateGenerator fails to build.
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
        msg = (f"--parametrize-ligands requires openmmforcefields and "
               f"openff-toolkit ({e}). Install them or drop the flag.")
        if strict:
            raise LigandParamError(msg)
        print(f"WARNING: {msg}; skipping ligand parametrisation")
        return None

    # AmberTools' sqm/antechamber must be on PATH so
    # AmberToolsToolkitWrapper registers itself in the OpenFF toolkit
    # registry (GAFFTemplateGenerator delegates AM1-BCC charge assignment
    # to it). When dvbfixer is invoked via the env's `bin/dvbfixer`
    # shim, the env's bin dir isn't automatically on PATH — prepend it
    # here so the wrapper's is_available() check succeeds.
    import sys as _sys
    _env_bin = os.path.dirname(_sys.executable)
    if _env_bin and _env_bin not in os.environ.get('PATH', '').split(os.pathsep):
        os.environ['PATH'] = _env_bin + os.pathsep + os.environ.get('PATH', '')

    import shutil as _sh
    for _tool in ('sqm', 'antechamber', 'parmchk2'):
        if _sh.which(_tool) is None:
            msg = (f"--parametrize-ligands: AmberTools binary '{_tool}' "
                   f"not found on PATH. Install ambertools "
                   f"(`micromamba install -n dvbfixer -c conda-forge ambertools`) "
                   f"or drop the flag.")
            if strict:
                raise LigandParamError(msg)
            print(f"WARNING: {msg}; skipping ligand parametrisation")
            return None

    # Force-register the AmberTools toolkit into the global registry so
    # GAFFTemplateGenerator finds it. Some environments load the OpenFF
    # registry before PATH is fixed above; explicit registration ensures
    # AmberToolsToolkitWrapper is present regardless of import order.
    try:
        from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
        from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
        _wrapper = AmberToolsToolkitWrapper()
        _registered = {type(w).__name__ for w in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits}
        if 'AmberToolsToolkitWrapper' not in _registered:
            GLOBAL_TOOLKIT_REGISTRY.register_toolkit(_wrapper)
    except Exception as e:
        msg = (f"--parametrize-ligands: could not register "
               f"AmberToolsToolkitWrapper ({e}). Charge assignment "
               f"will fall back and probably fail.")
        if strict:
            raise LigandParamError(msg)
        print(f"WARNING: {msg}")

    cache_dir = Path(cache_dir) if cache_dir else _default_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / 'gaff_ligands.json'

    src_pdb = _dump_topology_to_pdb(topology, positions)

    molecules = []
    failures = []
    with tempfile.TemporaryDirectory(prefix='dvbfixer_lig_') as tmp:
        for chain_id, res_id, resname in unknown:
            sdf_path = Path(tmp) / f'{resname}_{chain_id}_{res_id}.sdf'
            if not _extract_residue_sdf(src_pdb, chain_id, res_id,
                                         sdf_path, verbose=verbose):
                failures.append(
                    f"{resname} (chain {chain_id!r}, res {res_id}): "
                    f"OpenBabel extraction failed"
                )
                continue
            try:
                mol = Molecule.from_file(str(sdf_path),
                                         allow_undefined_stereo=True)
                mol.name = resname
                molecules.append(mol)
                if verbose:
                    print(f"  [lig_params] extracted {resname} "
                          f"({mol.n_atoms} atoms)")
            except Exception as e:
                failures.append(
                    f"{resname} (chain {chain_id!r}, res {res_id}): "
                    f"openff.toolkit could not read the SDF ({e})"
                )

        if failures:
            hint = ("Check that the ligand has complete heavy atoms and "
                    "a plausible bond graph; drop the flag or pass "
                    "--strip-heterogens if the ligand is not needed.")
            msg = ("--parametrize-ligands failed to parametrise:\n  "
                   + "\n  ".join(failures)
                   + f"\n{hint}")
            if strict:
                raise LigandParamError(msg)
            print(f"WARNING: {msg}")
            if not molecules:
                return None

        if not molecules:
            if strict:
                raise LigandParamError(
                    "--parametrize-ligands: no ligands successfully "
                    "extracted from the topology"
                )
            return None

        try:
            gen = GAFFTemplateGenerator(molecules=molecules,
                                        forcefield='gaff-2.11',
                                        cache=str(cache_file))
        except Exception as e:
            msg = (f"GAFFTemplateGenerator failed to build "
                   f"({e}). Check that AmberTools (antechamber, parmchk2) "
                   f"are installed and on PATH.")
            if strict:
                raise LigandParamError(msg)
            print(f"WARNING: {msg}; skipping ligand parametrisation")
            return None

    if verbose:
        print(f"  [lig_params] built GAFF2 templates for: "
              f"{', '.join(m.name for m in molecules)} "
              f"(cache: {cache_file})")
    return gen
