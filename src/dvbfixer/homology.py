"""Multi-template homology modeling with Modeller.

Takes a target FASTA (multi-chain) and one or more template PDB files.
Auto-aligns target to templates, builds a composite model using Modeller's
automodel/LoopModel with multiple knowns, then optionally runs prepare+minimize.

Antibody mode (--antibody): uses ANARCI for numbering/CDR detection.
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

from dvbfixer.model import (
    AA3TO1,
    WATER_RESNAMES,
    parse_pir_sequence,
    remove_water_lines,
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
        from anarci import anarci
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
# Per-chain alignment
# ---------------------------------------------------------------------------

def _extract_chain_pdb(pdb_path, chain_id, output_path):
    """Extract a single chain from a PDB file."""
    with open(pdb_path) as f:
        lines = f.readlines()
    with open(output_path, 'w') as f:
        for line in lines:
            if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
                if line[21] == chain_id:
                    f.write(line)
            elif line.startswith(('HEADER', 'CRYST1', 'REMARK', 'END')):
                f.write(line)


def _build_per_chain_alignment(env, target_chains, template_names,
                               chain_mapping, workdir, use_salign=False,
                               verbose=False):
    """Build multi-chain alignment by aligning each target chain to its best
    template chain independently, then combining into a single PIR file.

    The key insight: Modeller's PIR format for multi-template requires each
    template entry to have the SAME number of chain breaks (/) as the target.
    For templates that don't cover a target chain, that chain position gets
    all-gap characters. The structureX header must specify FIRST:LAST residue
    and chain so Modeller reads the correct atoms from the PDB.

    Returns path to the combined alignment PIR file.
    """
    from modeller import Alignment, Model

    # Step 1: Get template chain sequences as Modeller sees them
    # (important: Modeller may read atoms differently than our parser)
    tpl_modeller_seqs = {}  # tpl_name -> {chain_id: sequence}
    for tpl_name in template_names:
        mdl = Model(env, file=f'{tpl_name}.pdb')
        aln_tmp = Alignment(env)
        aln_tmp.append_model(mdl, align_codes=tpl_name,
                             atom_files=f'{tpl_name}.pdb')
        tmp_pir = os.path.join(workdir, f'_tpl_{tpl_name}.pir')
        aln_tmp.write(file=tmp_pir, alignment_format='PIR')
        full_seq = parse_pir_sequence(tmp_pir, tpl_name)
        # Split by chain breaks
        chain_seqs = full_seq.split('/')
        # Get chain IDs from PDB
        chain_ids = []
        seen = set()
        with open(os.path.join(workdir, f'{tpl_name}.pdb')) as pf:
            for line in pf:
                if line.startswith(('ATOM  ', 'HETATM')):
                    ch = line[21]
                    if ch not in seen:
                        seen.add(ch)
                        chain_ids.append(ch)
        tpl_modeller_seqs[tpl_name] = {}
        for i, ch_id in enumerate(chain_ids):
            if i < len(chain_seqs):
                tpl_modeller_seqs[tpl_name][ch_id] = chain_seqs[i]

    # Step 2: Pairwise alignment for each target chain
    per_chain = {}
    for tgt_chain_id, tgt_seq in target_chains:
        tpl_name, tpl_chain_id, _ = chain_mapping[tgt_chain_id]

        # Get template chain sequence as Modeller sees it
        tpl_chain_seq = tpl_modeller_seqs[tpl_name].get(tpl_chain_id, '')
        if not tpl_chain_seq:
            print(f"WARNING: template {tpl_name} chain {tpl_chain_id} not found",
                  file=sys.stderr)
            # Use target as-is (all modeled)
            per_chain[tgt_chain_id] = {
                'tpl_name': tpl_name,
                'tpl_chain': tpl_chain_id,
                'tpl_aligned': '-' * len(tgt_seq),
                'tgt_aligned': tgt_seq,
            }
            continue

        # Extract single chain from template PDB
        chain_pdb = os.path.join(workdir, f'_chain_{tpl_name}_{tpl_chain_id}.pdb')
        _extract_chain_pdb(
            os.path.join(workdir, f'{tpl_name}.pdb'),
            tpl_chain_id, chain_pdb)

        chain_code = f'{tpl_name}_{tpl_chain_id}'

        # Write single-chain target PIR
        chain_target_pir = os.path.join(workdir, f'_target_{tgt_chain_id}.pir')
        with open(chain_target_pir, 'w') as cf:
            cf.write(f">P1;target_{tgt_chain_id}\n")
            cf.write(f"sequence:target_{tgt_chain_id}::::::::\n")
            seq_with_star = tgt_seq if tgt_seq.endswith('*') else tgt_seq + '*'
            for i in range(0, len(seq_with_star), 75):
                cf.write(seq_with_star[i:i+75] + '\n')

        # Align single chain pair
        aln = Alignment(env)
        mdl = Model(env, file=chain_pdb)
        aln.append_model(mdl, align_codes=chain_code, atom_files=chain_pdb)
        aln.append(file=chain_target_pir, align_codes=f'target_{tgt_chain_id}')

        if use_salign:
            aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50),
                       alignment_type='tree', output='')
        else:
            aln.align2d(max_gap_length=50)

        pair_aln_path = os.path.join(workdir, f'_aln_{tgt_chain_id}.pir')
        aln.write(file=pair_aln_path, alignment_format='PIR')

        tpl_aln_seq = parse_pir_sequence(pair_aln_path, chain_code)
        tgt_aln_seq = parse_pir_sequence(pair_aln_path, f'target_{tgt_chain_id}')

        per_chain[tgt_chain_id] = {
            'tpl_name': tpl_name,
            'tpl_chain': tpl_chain_id,
            'tpl_aligned': tpl_aln_seq,
            'tgt_aligned': tgt_aln_seq,
        }

        if verbose:
            print(f"  Aligned chain {tgt_chain_id} → {tpl_name}:{tpl_chain_id} "
                  f"({len(tpl_aln_seq)} aligned positions)")

    # Step 3: Build combined multi-chain PIR
    # Each template entry must have the same number of '/' as the target.
    # Use Modeller's structureX header with FIRST and LAST residue/chain
    # so Modeller reads the correct PDB atoms.
    aln_path = os.path.join(workdir, 'alignment.pir')

    with open(aln_path, 'w') as f:
        for tpl_name in template_names:
            chain_seqs = []
            for tgt_chain_id, _ in target_chains:
                info = per_chain[tgt_chain_id]
                if info['tpl_name'] == tpl_name:
                    chain_seqs.append(info['tpl_aligned'])
                else:
                    gap_len = len(per_chain[tgt_chain_id]['tgt_aligned'])
                    chain_seqs.append('-' * gap_len)

            combined_seq = '/'.join(chain_seqs)

            # structureX header: FIRST and LAST specify the range Modeller
            # should read from the PDB. Use empty fields to let Modeller
            # auto-detect from the PDB file.
            f.write(f">P1;{tpl_name}\n")
            f.write(f"structureX:{tpl_name}.pdb:FIRST:@:LAST:@::::\n")
            for i in range(0, len(combined_seq), 75):
                f.write(combined_seq[i:i+75] + '\n')
            f.write('*\n\n')

        # Target entry
        chain_seqs = [per_chain[ch]['tgt_aligned'] for ch, _ in target_chains]
        combined_target = '/'.join(chain_seqs)
        f.write(">P1;target\n")
        f.write("sequence:target::::::::\n")
        for i in range(0, len(combined_target), 75):
            f.write(combined_target[i:i+75] + '\n')
        f.write('*\n')

    return aln_path


# ---------------------------------------------------------------------------
# Modeller multi-template modeling
# ---------------------------------------------------------------------------

def run_homology_modeller(target_chains, template_paths, chain_mapping,
                          args, workdir='.'):
    """Run Modeller multi-template homology modeling.

    Returns path to best model PDB.
    """
    from modeller import Alignment, Environ, Model, log
    from modeller import automodel as am
    from modeller.automodel import LoopModel, automodel

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

    md_levels = {
        "none": None,
        "fast": am.refine.fast,
        "slow": am.refine.slow,
        "very_slow": am.refine.very_slow,
        "slow_large": am.refine.slow_large,
    }

    if args.alignment:
        # User-provided alignment — single Modeller run
        aln_path = str(Path(args.alignment).resolve())
        shutil.copy2(aln_path, os.path.join(workdir, 'alignment.pir'))
        aln_path = os.path.join(workdir, 'alignment.pir')

        ModelClass = automodel if args.no_loop_refine else LoopModel
        a = ModelClass(env, alnfile=aln_path,
                       knowns=tuple(template_names), sequence='target')
        a.starting_model = 1
        a.ending_model = args.num_models
        if not args.no_loop_refine:
            a.loop.starting_model = 1
            a.loop.ending_model = args.num_models
            a.loop.md_level = md_levels[args.md_level]
        a.make()

        models = (_get_best_models(a, args.no_loop_refine))
        best = min(models, key=lambda x: x['molpdf'])
        print(f"Best model: {best['name']} (molpdf={best['molpdf']:.1f})")
        return os.path.join(workdir, best['name']), aln_path

    # Per-chain modeling: model each target chain independently with its
    # best template, then assemble into a full multi-chain PDB.
    # This is necessary because templates have different chain structures
    # (e.g. Fab=2 chains vs IgG=4 chains) and Modeller requires matching
    # chain break counts between alignment entries.
    chain_models = {}  # tgt_chain_id -> model PDB path

    for tgt_chain_id, tgt_seq in target_chains:
        tpl_name, tpl_chain_id, _ = chain_mapping[tgt_chain_id]
        print(f"\n  Modeling chain {tgt_chain_id} "
              f"(template {tpl_name}:{tpl_chain_id})...")

        # Extract single chain from template
        chain_pdb = os.path.join(workdir, f'_chain_{tpl_name}_{tpl_chain_id}.pdb')
        _extract_chain_pdb(
            os.path.join(workdir, f'{tpl_name}.pdb'),
            tpl_chain_id, chain_pdb)

        chain_code = f'{tpl_name}_{tpl_chain_id}'

        # Write target PIR for this chain
        chain_target_pir = os.path.join(workdir, f'_target_{tgt_chain_id}.pir')
        with open(chain_target_pir, 'w') as cf:
            cf.write(f">P1;target_{tgt_chain_id}\n")
            cf.write(f"sequence:target_{tgt_chain_id}::::::::\n")
            seq_star = tgt_seq + '*' if not tgt_seq.endswith('*') else tgt_seq
            for i in range(0, len(seq_star), 75):
                cf.write(seq_star[i:i+75] + '\n')

        # Align
        aln = Alignment(env)
        mdl = Model(env, file=chain_pdb)
        aln.append_model(mdl, align_codes=chain_code, atom_files=chain_pdb)
        aln.append(file=chain_target_pir, align_codes=f'target_{tgt_chain_id}')

        if args.salign:
            aln.salign(auto_overhang=True, gap_penalties_1d=(-450, -50),
                       alignment_type='tree', output='')
        else:
            aln.align2d(max_gap_length=50)

        chain_aln_path = os.path.join(workdir, f'_aln_{tgt_chain_id}.pir')
        aln.write(file=chain_aln_path, alignment_format='PIR')

        if args.verbose:
            with open(chain_aln_path) as af:
                print(af.read()[:300])

        # Model this chain
        ModelClass = automodel if args.no_loop_refine else LoopModel
        a = ModelClass(env, alnfile=chain_aln_path,
                       knowns=chain_code,
                       sequence=f'target_{tgt_chain_id}')
        a.starting_model = 1
        a.ending_model = args.num_models
        if not args.no_loop_refine:
            a.loop.starting_model = 1
            a.loop.ending_model = args.num_models
            a.loop.md_level = md_levels[args.md_level]
        a.make()

        models = _get_best_models(a, args.no_loop_refine)
        best = min(models, key=lambda x: x['molpdf'])
        print(f"    Best: {best['name']} (molpdf={best['molpdf']:.1f})")
        chain_models[tgt_chain_id] = os.path.join(workdir, best['name'])

    # Assemble per-chain models into a single multi-chain PDB
    assembled_path = os.path.join(workdir, 'assembled.pdb')
    _assemble_chains(chain_models, target_chains, assembled_path)
    print(f"\nAssembled {len(chain_models)} chains into {assembled_path}")

    return assembled_path, None


def _get_best_models(automodel_obj, no_loop_refine):
    """Extract successful models from automodel output."""
    if not no_loop_refine:
        models = [x for x in automodel_obj.loop.outputs if x['failure'] is None]
        if not models:
            models = [x for x in automodel_obj.outputs if x['failure'] is None]
            if models:
                print("WARNING: loop refinement failed, using initial model")
    else:
        models = [x for x in automodel_obj.outputs if x['failure'] is None]
    if not models:
        print("ERROR: all models failed", file=sys.stderr)
        sys.exit(1)
    return models


def _assemble_chains(chain_models, target_chains, output_path):
    """Assemble per-chain model PDBs into a single multi-chain PDB."""
    serial = 0
    with open(output_path, 'w') as f:
        f.write("REMARK    Assembled by dvbfixer homology\n")
        for tgt_chain_id, _ in target_chains:
            model_path = chain_models[tgt_chain_id]
            with open(model_path) as mf:
                for line in mf:
                    if line.startswith(('ATOM  ', 'HETATM')):
                        serial += 1
                        # Set chain ID and renumber serial
                        line = (f"{line[:6]}{serial % 100000:5d}"
                                f"{line[11:21]}{tgt_chain_id}"
                                f"{line[22:]}")
                        f.write(line)
            # TER between chains
            serial += 1
            f.write(f"TER   {serial % 100000:5d}\n")
        f.write("END\n")


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

        # 6. Post-process
        # Modeller outputs correct chain IDs from the PIR alignment,
        # so restore_chain_ids is only needed if they don't match.
        with open(best_model_path) as f:
            result_lines = f.readlines()

        # Check if chain IDs need fixing
        expected_chains = [ch for ch, _ in target_chains]
        actual_chains = []
        seen = set()
        for line in result_lines:
            if line.startswith(('ATOM  ', 'HETATM')):
                ch = line[21]
                if ch not in seen:
                    seen.add(ch)
                    actual_chains.append(ch)

        if actual_chains != expected_chains:
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
        from dvbfixer.ffutils.dat import DatRecord

        dat_path = f"{output_prefix}_homology.dat"
        DatRecord(
            description="Homology model built by dvbfixer homology",
            templates=[Path(t).name for t in template_paths],
            target_chains={ch: len(seq) for ch, seq in target_chains},
        ).save(dat_path, verbose=False)
        print(f"Wrote {dat_path}")

    finally:
        os.chdir(orig_dir)
        if args.keep_workdir:
            print(f"Modeller workdir: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    # 9. Optional pipeline
    if args.prepare or args.minimize:
        print("\nRunning prepare...")
        from dvbfixer.prepare import main as prepare_main
        prepare_args = [output_pdb, '--ph', str(args.ph)]
        if args.verbose:
            prepare_args.append('-v')
        prepare_main(prepare_args)

        prepared_pdb = output_pdb.replace('.pdb', '_prepared.pdb')
        if args.minimize and Path(prepared_pdb).exists():
            print("\nRunning minimize...")
            from dvbfixer.minimize import main as minimize_main
            minimize_args = [prepared_pdb, '--ph', str(args.ph)]
            if args.verbose:
                minimize_args.append('-v')
            minimize_main(minimize_args)


if __name__ == "__main__":
    main()
