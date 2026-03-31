"""Multi-template homology modeling with Modeller.

Takes a target FASTA (multi-chain) and one or more template PDB files.
Auto-aligns target to templates, builds a composite model using Modeller's
automodel/LoopModel with multiple knowns, then optionally runs prepare+minimize.

Antibody mode (--antibody): uses ANARCI for numbering/CDR detection.
"""

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from dvbfixer.model import (
    AA3TO1, WATER_RESNAMES, write_target_pir, parse_pir_sequence,
    restore_chain_ids_and_read, build_resnum_mapping, renumber_model_output,
    build_model_dat, remove_water_lines, get_template_mask,
)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta_chains(fasta_path):
    """Parse multi-chain FASTA. Returns list of (chain_id, sequence) tuples.

    Chain ID is extracted from the header: >ChainName or >Protein_X
    The last character of the header ID is used as chain ID if it's a single letter.
    """
    chains = []
    current_id = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    chains.append((current_id, ''.join(current_seq)))
                header = line[1:].split()[0]
                # Extract chain ID: last char after underscore, or last char
                if '_' in header:
                    chain_id = header.split('_')[-1]
                    if len(chain_id) == 1:
                        pass  # single letter chain ID
                    else:
                        chain_id = header[-1]
                else:
                    chain_id = header[-1] if len(header) == 1 else header[0]
                current_id = chain_id
                current_seq = []
            elif current_id is not None:
                current_seq.append(line)
    if current_id is not None:
        chains.append((current_id, ''.join(current_seq)))
    return chains


# ---------------------------------------------------------------------------
# Template analysis
# ---------------------------------------------------------------------------

def get_template_chains(pdb_path):
    """Extract chain sequences from a PDB template file.

    Returns dict: chain_id -> one-letter sequence.
    """
    chains = {}
    current_chain = None
    current_resseq = None
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
            chain = line[21]
            resname = line[17:20].strip()
            resseq = line[22:27].strip()  # includes icode
            if resname in WATER_RESNAMES:
                continue
            aa = AA3TO1.get(resname)
            if aa is None:
                continue
            key = (chain, resseq)
            if key != current_resseq:
                current_resseq = key
                if chain != current_chain:
                    current_chain = chain
                    chains.setdefault(chain, [])
                chains[chain].append(aa)
    return {k: ''.join(v) for k, v in chains.items()}


def sequence_identity(seq1, seq2):
    """Compute approximate sequence identity using k-mer overlap.

    Fast heuristic — counts shared 3-mers normalized by sequence length.
    Sufficient for template selection (not alignment).
    """
    k = 3
    if len(seq1) < k or len(seq2) < k:
        # Fallback for very short sequences
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / max(len(seq1), len(seq2)) if max(len(seq1), len(seq2)) > 0 else 0.0

    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    shared = len(kmers1 & kmers2)
    total = max(len(kmers1), len(kmers2))
    return shared / total if total > 0 else 0.0


def map_target_to_templates(target_chains, template_infos, verbose=False):
    """Map each target chain to the best-matching template chain.

    target_chains: [(chain_id, sequence), ...]
    template_infos: [(template_name, {chain_id: sequence}), ...]

    Returns: dict target_chain_id -> (template_name, template_chain_id, identity)
    """
    mapping = {}
    for tgt_id, tgt_seq in target_chains:
        best_match = None
        best_score = -1
        for tpl_name, tpl_chains in template_infos:
            for tpl_chain_id, tpl_seq in tpl_chains.items():
                score = sequence_identity(tgt_seq, tpl_seq)
                if score > best_score:
                    best_score = score
                    best_match = (tpl_name, tpl_chain_id, score)
        if best_match:
            mapping[tgt_id] = best_match
            if verbose:
                print(f"  Target chain {tgt_id} → {best_match[0]}:{best_match[1]} "
                      f"({best_match[2]*100:.1f}% identity)")
    return mapping


# ---------------------------------------------------------------------------
# Antibody mode (ANARCI)
# ---------------------------------------------------------------------------

def run_antibody_analysis(target_chains, verbose=False):
    """Run ANARCI on target chains for antibody numbering and CDR detection.

    Returns dict with domain info per chain.
    """
    try:
        from anarci import anarci, run_anarci
    except ImportError:
        print("ERROR: ANARCI not installed. Install with: pip install anarci",
              file=sys.stderr)
        sys.exit(1)

    results = {}
    sequences = [(ch_id, seq) for ch_id, seq in target_chains]

    for chain_id, seq in sequences:
        numbered, alignment_details, hit_tables = anarci(
            [("query", seq)], scheme="kabat", output=False)

        if numbered[0] is None:
            if verbose:
                print(f"  Chain {chain_id}: not an antibody chain")
            results[chain_id] = {'type': 'non-antibody'}
            continue

        # Extract domain info
        # numbered[0][0] = (numbering_list, start_index, end_index)
        # chain type is in alignment_details
        numbering, start_idx, end_idx = numbered[0][0]
        domain_type = alignment_details[0][0]['chain_type']  # H, K, L

        # Determine CDR boundaries (Kabat definition)
        cdr_ranges = {}
        if domain_type == 'H':
            cdr_ranges = {
                'CDR-H1': (26, 35),
                'CDR-H2': (50, 65),
                'CDR-H3': (95, 102),
            }
            domain_label = 'VH'
        elif domain_type in ('K', 'L'):
            cdr_ranges = {
                'CDR-L1': (24, 34),
                'CDR-L2': (50, 56),
                'CDR-L3': (89, 97),
            }
            domain_label = 'VL'
        else:
            domain_label = domain_type

        results[chain_id] = {
            'type': domain_label,
            'chain_type': domain_type,
            'numbering': numbering,
            'cdrs': cdr_ranges,
        }

        if verbose:
            print(f"  Chain {chain_id}: {domain_label} ({domain_type})")
            for cdr_name, (start, end) in cdr_ranges.items():
                cdr_seq = ''.join(
                    aa for (pos, icode), aa in numbering
                    if start <= pos <= end and aa != '-'
                )
                print(f"    {cdr_name}: {cdr_seq}")

    return results


# ---------------------------------------------------------------------------
# Modeller multi-template modeling
# ---------------------------------------------------------------------------

def run_homology_modeller(target_chains, template_paths, chain_mapping,
                          args, workdir='.'):
    """Run Modeller multi-template homology modeling.

    Returns path to best model PDB.
    """
    from modeller import Environ, Alignment, Model, log
    from modeller.automodel import automodel, LoopModel
    from modeller import automodel as am

    if args.verbose:
        log.verbose()
    else:
        log.none()

    env = Environ()
    env.io.atom_files_directory = [workdir, '.']
    env.io.hetatm = True

    # Copy template PDBs to workdir with clean names
    template_names = []
    for i, tpl_path in enumerate(template_paths):
        tpl_name = Path(tpl_path).stem
        dst = os.path.join(workdir, f'{tpl_name}.pdb')
        if not os.path.exists(dst):
            shutil.copy2(tpl_path, dst)
        template_names.append(tpl_name)

    # Build target sequence in PIR format (multi-chain with '/')
    target_seq = '/'.join(seq for _, seq in target_chains) + '*'
    target_pir = os.path.join(workdir, 'target.pir')
    write_target_pir(target_seq, target_pir)

    if args.alignment:
        # User-provided alignment
        aln_path = str(Path(args.alignment).resolve())
        shutil.copy2(aln_path, os.path.join(workdir, 'alignment.pir'))
        aln_path = os.path.join(workdir, 'alignment.pir')
    else:
        # Auto-align: build alignment from templates + target
        aln = Alignment(env)
        for tpl_name in template_names:
            mdl = Model(env, file=f'{tpl_name}.pdb')
            aln.append_model(mdl, align_codes=tpl_name,
                             atom_files=f'{tpl_name}.pdb')

        aln.append(file=target_pir, align_codes='target')

        if args.salign:
            aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50),
                       alignment_type='tree', output='ALIGNMENT')
        else:
            aln.align2d(max_gap_length=50)

        aln_path = os.path.join(workdir, 'alignment.pir')
        aln.write(file=aln_path, alignment_format='PIR')

    if args.verbose:
        print(f"Alignment written to {aln_path}")
        with open(aln_path) as f:
            print(f.read()[:500])

    # Select model class
    md_levels = {
        "none": None,
        "fast": am.refine.fast,
        "slow": am.refine.slow,
        "very_slow": am.refine.very_slow,
        "slow_large": am.refine.slow_large,
    }

    if args.no_loop_refine:
        ModelClass = automodel
    else:
        ModelClass = LoopModel

    a = ModelClass(env,
                   alnfile=aln_path,
                   knowns=tuple(template_names),
                   sequence='target')
    a.starting_model = 1
    a.ending_model = args.num_models

    if not args.no_loop_refine:
        a.loop.starting_model = 1
        a.loop.ending_model = args.num_models
        a.loop.md_level = md_levels[args.md_level]

    a.make()

    # Select best model
    if not args.no_loop_refine:
        models = [x for x in a.loop.outputs if x['failure'] is None]
        if not models:
            models = [x for x in a.outputs if x['failure'] is None]
            if models:
                print("WARNING: loop refinement failed, using initial model")
    else:
        models = [x for x in a.outputs if x['failure'] is None]

    if not models:
        print("ERROR: all models failed", file=sys.stderr)
        sys.exit(1)

    best = min(models, key=lambda x: x['molpdf'])
    print(f"Best model: {best['name']} (molpdf={best['molpdf']:.1f})")

    return os.path.join(workdir, best['name']), aln_path


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def restore_chain_ids(model_path, target_chains):
    """Restore chain IDs from Modeller's A,B,C,... to target chain IDs."""
    chain_ids = [ch_id for ch_id, _ in target_chains]
    with open(model_path) as f:
        lines = f.readlines()

    # Modeller assigns chains as A, B, C, ... in order
    modeller_chains = []
    seen = set()
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            ch = line[21]
            if ch not in seen:
                seen.add(ch)
                modeller_chains.append(ch)

    # Build mapping: modeller chain -> target chain
    chain_map = {}
    for i, mod_ch in enumerate(modeller_chains):
        if i < len(chain_ids):
            chain_map[mod_ch] = chain_ids[i]
        else:
            chain_map[mod_ch] = mod_ch

    # Apply mapping
    result = []
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
            old_ch = line[21]
            new_ch = chain_map.get(old_ch, old_ch)
            line = line[:21] + new_ch + line[22:]
        result.append(line)

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog='dvbfixer homology',
        description='Multi-template homology modeling with Modeller. '
                    'Builds a composite model from multiple template structures.',
    )
    p.add_argument('fasta',
                   help='Target sequence FASTA (multi-chain, one >header per chain)')
    p.add_argument('--template', action='append', required=True,
                   help='Template PDB file (repeatable, at least 1)')
    p.add_argument('--alignment', default=None,
                   help='Pre-built PIR alignment file (skip auto-alignment)')
    p.add_argument('--salign', action='store_true',
                   help='Use structure-based alignment instead of align2d')
    p.add_argument('-n', '--num-models', type=int, default=5,
                   help='Number of models to generate (default: 5)')
    p.add_argument('--md-level',
                   choices=['none', 'fast', 'slow', 'very_slow', 'slow_large'],
                   default='fast',
                   help='MD refinement level (default: fast)')
    p.add_argument('--no-loop-refine', action='store_true',
                   help='Use automodel instead of LoopModel (faster, no loop refinement)')
    p.add_argument('--antibody', action='store_true',
                   help='Antibody-aware mode: ANARCI numbering, CDR detection')
    p.add_argument('--prepare', action='store_true',
                   help='Run dvbfixer prepare on output')
    p.add_argument('--minimize', action='store_true',
                   help='Run dvbfixer prepare + minimize on output')
    p.add_argument('--ph', type=float, default=7.0,
                   help='pH for hydrogen addition (default: 7.0)')
    p.add_argument('-o', '--output', default=None,
                   help='Output prefix (default: FASTA stem)')
    p.add_argument('--keep-workdir', action='store_true',
                   help='Keep Modeller working directory')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Verbose output')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    fasta_path = Path(args.fasta).resolve()
    output_prefix = args.output or fasta_path.stem

    if not fasta_path.exists():
        print(f"Error: {fasta_path} not found", file=sys.stderr)
        sys.exit(1)

    for tpl in args.template:
        if not Path(tpl).exists():
            print(f"Error: template {tpl} not found", file=sys.stderr)
            sys.exit(1)

    # 1. Parse target FASTA
    target_chains = parse_fasta_chains(fasta_path)
    print(f"Target: {len(target_chains)} chain(s)")
    for ch_id, seq in target_chains:
        print(f"  Chain {ch_id}: {len(seq)} residues")

    # 2. Analyze templates
    template_infos = []
    template_paths = [str(Path(t).resolve()) for t in args.template]
    for tpl_path in template_paths:
        tpl_name = Path(tpl_path).stem
        tpl_chains = get_template_chains(tpl_path)
        template_infos.append((tpl_name, tpl_chains))
        if args.verbose:
            print(f"Template {tpl_name}: {len(tpl_chains)} chain(s)")
            for ch, seq in tpl_chains.items():
                print(f"  Chain {ch}: {len(seq)} residues")

    # 3. Map target chains to best template chains
    print("\nChain mapping (target → template):")
    chain_mapping = map_target_to_templates(
        target_chains, template_infos, verbose=True)

    # 4. Antibody analysis (optional)
    if args.antibody:
        print("\nAntibody analysis (ANARCI):")
        ab_info = run_antibody_analysis(target_chains, verbose=args.verbose)

    # 5. Run Modeller
    workdir = tempfile.mkdtemp(prefix='dvbfixer_homology_')
    orig_dir = os.getcwd()

    try:
        os.chdir(workdir)
        print(f"\nRunning Modeller ({args.num_models} models, "
              f"md_level={args.md_level})...")

        best_model_path, aln_path = run_homology_modeller(
            target_chains, template_paths, chain_mapping, args, workdir)

        # 6. Post-process: restore chain IDs
        print("Restoring chain IDs...")
        result_lines = restore_chain_ids(best_model_path, target_chains)

        # Remove water if present
        result_lines = remove_water_lines(result_lines)

        # 7. Write output
        os.chdir(orig_dir)
        output_pdb = f"{output_prefix}_homology.pdb"
        with open(output_pdb, 'w') as f:
            f.writelines(result_lines)
        print(f"Wrote {output_pdb}")

        # 8. Write .dat file (record modeled regions)
        # For homology modeling, all atoms are "new" since the entire structure
        # is built. But we mark template-covered regions as "original" and
        # gap regions as "new" for restraint purposes.
        dat_path = f"{output_prefix}_homology.dat"
        dat_info = {
            'description': 'Homology model built by dvbfixer homology',
            'templates': [Path(t).name for t in template_paths],
            'target_chains': {ch: len(seq) for ch, seq in target_chains},
            'total_added': 0,  # all atoms are modeled
            'added_atoms': [],
        }
        with open(dat_path, 'w') as f:
            json.dump(dat_info, f, indent=2)
        print(f"Wrote {dat_path}")

    finally:
        os.chdir(orig_dir)
        if args.keep_workdir:
            print(f"Modeller workdir: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    # 9. Optional pipeline
    if args.prepare or args.minimize:
        print(f"\nRunning prepare...")
        from dvbfixer.prepare import main as prepare_main
        prepare_args = [output_pdb, '--ph', str(args.ph)]
        if args.verbose:
            prepare_args.append('-v')
        prepare_main(prepare_args)

        prepared_pdb = output_pdb.replace('.pdb', '_prepared.pdb')
        if args.minimize and Path(prepared_pdb).exists():
            print(f"\nRunning minimize...")
            from dvbfixer.minimize import main as minimize_main
            minimize_args = [prepared_pdb, '--ph', str(args.ph)]
            if args.verbose:
                minimize_args.append('-v')
            minimize_main(minimize_args)


if __name__ == "__main__":
    main()
