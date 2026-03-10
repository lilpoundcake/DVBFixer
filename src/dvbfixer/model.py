"""Rebuild missing loops/gaps in a PDB structure using Modeller.

Uses the SEQRES section (or a user-provided FASTA) as the complete sequence,
aligns it to the existing ATOM records, and runs Modeller's LoopModel to
rebuild only the missing regions while keeping the rest of the structure fixed.

Non-protein chains (glycans, ligands) are included in the Modeller pipeline
via env.io.hetatm=True. The target sequence includes '.' (BLK residue) entries
for non-protein chains so Modeller preserves them through loop modeling.
The '.' counts are derived from Modeller's own template reading to avoid
any mismatch.
"""

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Nonstandard / protonation variants
    'MSE': 'M', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
    'HID': 'H', 'HIE': 'H', 'HIP': 'H',
    'CYX': 'C', 'CYM': 'C', 'ASH': 'D', 'GLH': 'E', 'LYN': 'K',
}

WATER_RESNAMES = {'HOH', 'WAT', 'TIP3', 'TIP', 'SOL', 'T3P', 'T4P', 'T5P'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="dvbfixer model",
        description="Rebuild missing loops and gaps in a PDB structure using Modeller. "
        "Identifies gaps from SEQRES vs ATOM records (or a provided FASTA), "
        "then uses Modeller's loop modeling to fill them."
    )
    p.add_argument("input", help="Input PDB file (must contain SEQRES or use --fasta)")
    p.add_argument("-o", "--output", help="Output PDB file (default: <input>_model.pdb)")
    p.add_argument(
        "--fasta", help="FASTA file with complete sequence(s). Chain order must match "
        "the PDB. Use instead of SEQRES."
    )
    p.add_argument(
        "-n", "--num-models", type=int, default=1,
        help="Number of initial models to generate (default: 1)"
    )
    p.add_argument(
        "--num-loops", type=int, default=2,
        help="Number of loop refinement models per initial model (default: 2)"
    )
    p.add_argument(
        "--md-level", choices=["none", "fast", "slow", "very_slow", "slow_large"],
        default="fast",
        help="MD refinement level for loop modeling (default: fast)"
    )
    p.add_argument(
        "--no-terminal", action="store_true",
        help="Do not model missing N/C terminal residues (only rebuild internal gaps)"
    )
    p.add_argument(
        "--keep-water", action="store_true",
        help="Keep water molecules (HOH, WAT, TIP3, SOL) in output (default: remove)"
    )
    p.add_argument(
        "--keep-workdir", action="store_true",
        help="Keep the Modeller working directory (for debugging)"
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print Modeller progress"
    )
    return p.parse_args(argv)


def parse_seqres(lines):
    """Return dict: chain_id -> [resname, ...] from SEQRES records."""
    seqres = {}
    for line in lines:
        if not line.startswith("SEQRES"):
            continue
        chain = line[11]
        seqres.setdefault(chain, []).extend(line[19:].split())
    return seqres


def parse_fasta(fasta_path):
    """Parse a FASTA file, return list of (header, sequence) tuples."""
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_path, "fasta"))
    return [(r.id, str(r.seq)) for r in records]


def get_chain_order(lines, seqres_only=None):
    """Get chain IDs from ATOM/HETATM records in order of appearance.

    If seqres_only is a set of chain IDs, only return chains in that set.
    """
    chains = []
    seen = set()
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            c = line[21]
            if c not in seen:
                seen.add(c)
                if seqres_only is None or c in seqres_only:
                    chains.append(c)
    return chains


def write_target_pir(sequence, pir_path):
    """Write target sequence in PIR format."""
    with open(pir_path, 'w') as f:
        f.write(">P1;target\n")
        f.write("sequence:target::::::::\n")
        for i in range(0, len(sequence), 75):
            f.write(sequence[i:i + 75] + '\n')


def parse_pir_sequence(pir_path, code):
    """Parse a PIR file and return the sequence string for a given code.

    Returns the raw sequence (with '/', '*', amino acid letters, and '.').
    """
    with open(pir_path) as f:
        content = f.read()

    blocks = content.split('>P1;')
    for block in blocks:
        if not block.strip():
            continue
        lines = block.strip().split('\n')
        if lines[0].strip() == code:
            seq_lines = lines[2:]  # skip code line and header line
            seq = ''.join(seq_lines).replace('\n', '')
            return seq.rstrip('*')
    return ''


def parse_alignment(aln_path):
    """Parse PIR alignment file to extract the template sequence (with gaps).

    Returns the template sequence string where '-' marks gap positions
    and '/' separates chains.
    """
    with open(aln_path) as f:
        content = f.read()

    blocks = content.split('>P1;')
    if len(blocks) < 2:
        return ''

    template_block = blocks[1]
    seq_lines = template_block.split('\n')[2:]
    seq = ''.join(seq_lines).rstrip('*').replace('\n', '')
    return seq


def count_gaps(aln_path):
    """Count gap residues in the template from a PIR alignment file."""
    return parse_alignment(aln_path).count('-')


def get_atom_sequence(lines, chain_id):
    """Get the one-letter amino acid sequence from ATOM records for a chain."""
    residues = []
    seen = set()
    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        ch = line[21]
        if ch != chain_id:
            continue
        resname = line[17:20].strip()
        resnum = line[22:26].strip()
        icode = line[26]
        key = (resnum, icode)
        if key not in seen:
            seen.add(key)
            aa = AA3TO1.get(resname)
            if aa:
                residues.append(aa)
    return ''.join(residues)


def trim_terminal_gaps(protein_seq_map, pdb_lines, verbose=False):
    """Trim N/C terminal residues from target sequences that aren't in the structure.

    Compares SEQRES-derived target sequence to ATOM-derived sequence for each chain.
    Finds the ATOM sequence as a subsequence of the target, then trims any target
    residues before the first matched position or after the last matched position.
    Internal gaps (missing loops) are preserved.
    """
    trimmed = {}
    for chain, target_seq in protein_seq_map.items():
        atom_seq = get_atom_sequence(pdb_lines, chain)
        if not atom_seq:
            trimmed[chain] = target_seq
            continue

        # Find ATOM residues as subsequence of target — record matched positions
        matched_positions = []
        ti = 0
        for ai in range(len(atom_seq)):
            while ti < len(target_seq) and target_seq[ti] != atom_seq[ai]:
                ti += 1
            if ti >= len(target_seq):
                break
            matched_positions.append(ti)
            ti += 1

        if len(matched_positions) != len(atom_seq):
            # Subsequence match failed, keep original
            trimmed[chain] = target_seq
            continue

        first_pos = matched_positions[0]
        last_pos = matched_positions[-1]
        new_seq = target_seq[first_pos:last_pos + 1]

        n_trimmed_n = first_pos
        n_trimmed_c = len(target_seq) - last_pos - 1
        if verbose and (n_trimmed_n > 0 or n_trimmed_c > 0):
            print(f"  Chain {chain}: trimmed {n_trimmed_n} N-terminal, "
                  f"{n_trimmed_c} C-terminal residues")

        trimmed[chain] = new_seq
    return trimmed


def _fix_terminal_alignment(raw_aln_path, fixed_aln_path, pdb_name,
                             protein_set, all_chains, template_chain_seqs,
                             target_chain_seqs, verbose=False):
    """Fix terminal gap placement in align2d output.

    align2d can misplace terminal gaps (e.g. matching last template residue K
    to last target residue D instead of putting the gap at the C-terminus).

    For each protein chain, builds a correct alignment by finding where the
    template protein residues sit within the target sequence. Non-protein
    chains keep align2d's result. Internal gaps from align2d are preserved.
    """
    with open(raw_aln_path) as f:
        content = f.read()

    blocks = content.split('>P1;')
    if len(blocks) < 3:
        import shutil
        shutil.copy2(raw_aln_path, fixed_aln_path)
        return

    tpl_block = blocks[1]
    tgt_block = blocks[2]

    tpl_lines = tpl_block.strip().split('\n')
    tgt_lines = tgt_block.strip().split('\n')

    tpl_header = tpl_lines[:2]
    tgt_header = tgt_lines[:2]

    raw_tpl_seq = ''.join(tpl_lines[2:]).replace('\n', '').rstrip('*')
    raw_tgt_seq = ''.join(tgt_lines[2:]).replace('\n', '').rstrip('*')

    raw_tpl_chains = raw_tpl_seq.split('/')
    raw_tgt_chains = raw_tgt_seq.split('/')

    fixed_tpl_chains = []
    fixed_tgt_chains = []

    for ci in range(len(raw_tpl_chains)):
        ch = all_chains[ci] if ci < len(all_chains) else '?'
        raw_tpl = raw_tpl_chains[ci]
        raw_tgt = raw_tgt_chains[ci] if ci < len(raw_tgt_chains) else ''

        if ch not in protein_set:
            fixed_tpl_chains.append(raw_tpl)
            fixed_tgt_chains.append(raw_tgt)
            continue

        # Extract pure protein residues from template (skip '.' for non-protein)
        tpl_pure = raw_tpl.replace('-', '')
        tgt_pure = raw_tgt.replace('-', '')

        # Only fix protein residues (letters), skip '.' entries
        tpl_protein = ''.join(c for c in tpl_pure if c != '.')
        tgt_protein = tgt_pure  # target is all protein for protein chains

        if not tpl_protein or not tgt_protein:
            fixed_tpl_chains.append(raw_tpl)
            fixed_tgt_chains.append(raw_tgt)
            continue

        # Find where template protein maps into target via subsequence matching
        match_positions = []
        ti = 0
        for c in tpl_protein:
            while ti < len(tgt_protein) and tgt_protein[ti] != c:
                ti += 1
            if ti >= len(tgt_protein):
                break
            match_positions.append(ti)
            ti += 1

        if len(match_positions) != len(tpl_protein):
            fixed_tpl_chains.append(raw_tpl)
            fixed_tgt_chains.append(raw_tgt)
            continue

        n_prefix = match_positions[0]
        n_suffix = len(tgt_protein) - match_positions[-1] - 1

        # Check: are there internal gaps (non-contiguous matches)?
        has_internal_gaps = False
        for j in range(1, len(match_positions)):
            if match_positions[j] != match_positions[j - 1] + 1:
                has_internal_gaps = True
                break

        # Build new alignment for this chain
        if not has_internal_gaps:
            # Simple case: template is a contiguous block within target
            # N-gap + template + C-gap
            new_tpl = '-' * n_prefix + tpl_pure + '-' * n_suffix
            new_tgt = tgt_protein
        else:
            # Internal gaps exist: use match_positions to build alignment
            new_tpl_chars = []
            new_tgt_chars = []
            tpl_idx = 0
            prev_tgt_pos = -1

            for j, tgt_pos in enumerate(match_positions):
                # Add unmatched target residues before this match as gaps in template
                if j == 0:
                    # N-terminal unmatched
                    for p in range(tgt_pos):
                        new_tpl_chars.append('-')
                        new_tgt_chars.append(tgt_protein[p])
                else:
                    # Internal unmatched between prev and current
                    for p in range(prev_tgt_pos + 1, tgt_pos):
                        new_tpl_chars.append('-')
                        new_tgt_chars.append(tgt_protein[p])

                new_tpl_chars.append(tpl_protein[j])
                new_tgt_chars.append(tgt_protein[tgt_pos])
                prev_tgt_pos = tgt_pos

            # C-terminal unmatched
            for p in range(prev_tgt_pos + 1, len(tgt_protein)):
                new_tpl_chars.append('-')
                new_tgt_chars.append(tgt_protein[p])

            new_tpl = ''.join(new_tpl_chars)
            new_tgt = ''.join(new_tgt_chars)

        if len(new_tpl) != len(new_tgt):
            if verbose:
                print(f"  Chain {ch}: alignment fix length mismatch, keeping original")
            fixed_tpl_chains.append(raw_tpl)
            fixed_tgt_chains.append(raw_tgt)
            continue

        if verbose and (new_tpl != raw_tpl or new_tgt != raw_tgt):
            print(f"  Chain {ch}: fixed terminal alignment "
                  f"(N-gap: {n_prefix}, C-gap: {n_suffix})")

        fixed_tpl_chains.append(new_tpl)
        fixed_tgt_chains.append(new_tgt)

    # Write fixed alignment
    fixed_tpl_seq = '/'.join(fixed_tpl_chains) + '*'
    fixed_tgt_seq = '/'.join(fixed_tgt_chains) + '*'

    with open(fixed_aln_path, 'w') as f:
        f.write('>P1;' + '\n'.join(tpl_header) + '\n')
        for i in range(0, len(fixed_tpl_seq), 75):
            f.write(fixed_tpl_seq[i:i + 75] + '\n')
        f.write('>P1;' + '\n'.join(tgt_header) + '\n')
        for i in range(0, len(fixed_tgt_seq), 75):
            f.write(fixed_tgt_seq[i:i + 75] + '\n')


def run_modeller(input_path, protein_chains, protein_seq_map, all_chains, args):
    """Run Modeller loop modeling. Returns (best_model_path, alignment_path).

    Uses env.io.hetatm=True so non-protein atoms (glycans, ligands) are
    read into the model and preserved through loop modeling.

    The target sequence for non-protein chains is derived from Modeller's own
    template reading to guarantee the '.' counts match exactly.
    """
    from modeller import Environ, Alignment, Model, log
    from modeller.automodel import LoopModel
    from modeller import automodel as am

    if args.verbose:
        log.verbose()
    else:
        log.none()

    env = Environ()
    env.io.atom_files_directory = [str(input_path.parent), '.']
    env.io.hetatm = True

    pdb_name = input_path.stem

    if not Path(f'{pdb_name}.pdb').exists():
        shutil.copy2(input_path, f'{pdb_name}.pdb')

    # Step 1: Load model and get template sequence as Modeller sees it
    mdl = Model(env, file=f'{pdb_name}.pdb')
    aln = Alignment(env)
    aln.append_model(mdl, align_codes=pdb_name, atom_files=f'{pdb_name}.pdb')
    aln.write(file='template_only.pir', alignment_format='PIR')

    # Parse the template sequence Modeller generated
    template_seq = parse_pir_sequence('template_only.pir', pdb_name)
    template_chain_seqs = template_seq.split('/')

    if args.verbose:
        print(f"Modeller template: {len(template_chain_seqs)} chain(s), "
              f"total {len(template_seq) - template_seq.count('/')} residues")

    # Step 2: Build target sequence using SEQRES for protein chains and
    # Modeller's own '.' sequence for non-protein chains (guaranteed match)
    protein_set = set(protein_chains)
    protein_idx = 0
    target_chain_seqs = []
    for ci, ch in enumerate(all_chains):
        if ci < len(template_chain_seqs):
            tpl_chain = template_chain_seqs[ci]
        else:
            tpl_chain = ''

        if ch in protein_set and ch in protein_seq_map:
            target_chain_seqs.append(protein_seq_map[ch])
        else:
            # Non-protein chain: use Modeller's template sequence verbatim
            target_chain_seqs.append(tpl_chain)

    target_seq = '/'.join(target_chain_seqs) + '*'
    write_target_pir(target_seq, 'target.pir')

    # Step 3: Create alignment between template and target
    aln2 = Alignment(env)
    aln2.append_model(mdl, align_codes=pdb_name, atom_files=f'{pdb_name}.pdb')
    aln2.append(file='target.pir', align_codes='target')
    aln2.align2d()
    aln2.write(file='alignment_raw.pir', alignment_format='PIR')

    # Fix terminal gap placement: align2d can misplace terminal gaps
    # (e.g. matching last template residue to last target residue even when
    # they're different). Rewrite alignment with gaps forced to termini.
    _fix_terminal_alignment(
        'alignment_raw.pir', 'alignment.pir',
        pdb_name, protein_set, all_chains, template_chain_seqs,
        target_chain_seqs, args.verbose,
    )

    aln_path = os.path.join(os.getcwd(), 'alignment.pir')
    n_gaps = count_gaps(aln_path)
    if n_gaps == 0:
        print("No gaps found — structure matches SEQRES completely.")
        return str(input_path), aln_path

    print(f"Alignment: {n_gaps} gap residue(s) to rebuild")

    md_levels = {
        "none": None,
        "fast": am.refine.fast,
        "slow": am.refine.slow,
        "very_slow": am.refine.very_slow,
        "slow_large": am.refine.slow_large,
    }

    a = LoopModel(env,
                  alnfile='alignment.pir',
                  knowns=pdb_name,
                  sequence='target')
    a.starting_model = 1
    a.ending_model = args.num_models
    a.loop.starting_model = 1
    a.loop.ending_model = args.num_loops
    a.loop.md_level = md_levels[args.md_level]

    a.make()

    loop_models = [x for x in a.loop.outputs if x['failure'] is None]
    if not loop_models:
        init_models = [x for x in a.outputs if x['failure'] is None]
        if not init_models:
            print("Error: all models failed", file=sys.stderr)
            sys.exit(1)
        best = min(init_models, key=lambda x: x['molpdf'])
        print("Warning: loop refinement failed, using initial model")
    else:
        best = min(loop_models, key=lambda x: x['molpdf'])

    print(f"Best model: {best['name']} (molpdf={best['molpdf']:.1f})")
    return best['name'], aln_path


# ---------------------------------------------------------------------------
# Post-processing: restore chain IDs and residue numbering
# ---------------------------------------------------------------------------

def _set_resid(line, resseq, icode=' '):
    """Set residue sequence number (cols 22-25) and insertion code (col 26)."""
    return line[:22] + f"{resseq:4d}" + icode + line[27:]


def get_template_mask(aln_path):
    """Parse alignment and return per-chain boolean masks.

    True = template position (present in original), False = gap (rebuilt).
    """
    template_seq = parse_alignment(aln_path)
    per_chain = []
    for chain_seq in template_seq.split('/'):
        per_chain.append([ch != '-' for ch in chain_seq])
    flat = []
    for cm in per_chain:
        flat.extend(cm)
    return flat, per_chain


def _get_original_resids_per_chain(original_lines, chain_order):
    """Get ordered list of unique (resSeq, iCode) per chain from original PDB."""
    chain_set = set(chain_order)
    result = {ch: [] for ch in chain_order}
    seen = {ch: set() for ch in chain_order}
    for line in original_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        ch = line[21]
        if ch not in chain_set:
            continue
        key = (int(line[22:26].strip()), line[26])
        if key not in seen[ch]:
            seen[ch].add(key)
            result[ch].append(key)
    return result


def build_resnum_mapping(per_chain_masks, all_chains, protein_chains, original_lines):
    """Build mapping: (model_chain, model_resSeq) -> (original_resSeq, icode).

    Template positions get the original (resSeq, iCode).
    Gap-filled positions get sequential numbers (with blank iCode) between flanking originals.
    If there's not enough room between flanking numbers, uses negative offsets from the
    right neighbor to avoid collisions.
    Only builds mapping for protein chains; non-protein chains are skipped.
    """
    protein_set = set(protein_chains)
    orig_resids = _get_original_resids_per_chain(original_lines, protein_chains)
    mapping = {}

    for ci, chain in enumerate(all_chains):
        if ci >= len(per_chain_masks):
            continue
        mask = per_chain_masks[ci]

        if chain not in protein_set:
            # Non-protein chain — no renumbering needed
            continue

        orig_rids = orig_resids.get(chain, [])

        orig_idx = 0
        full_resids = [None] * len(mask)  # Each entry is (resSeq, iCode)

        for i, is_template in enumerate(mask):
            if is_template and orig_idx < len(orig_rids):
                full_resids[i] = orig_rids[orig_idx]  # (resSeq, iCode)
                orig_idx += 1

        # Fill gap runs: sequential numbers between flanking known positions
        i = 0
        while i < len(mask):
            if mask[i]:
                i += 1
                continue
            gap_start = i
            while i < len(mask) and not mask[i]:
                i += 1
            gap_end = i
            gap_len = gap_end - gap_start

            left = full_resids[gap_start - 1][0] if gap_start > 0 else 0
            right = full_resids[gap_end][0] if gap_end < len(mask) else None

            if right is not None:
                # Internal gap or gap before a matched terminal: fit between left and right
                available = right - left - 1
                if available >= gap_len:
                    # Enough room — number sequentially from left
                    for k in range(gap_len):
                        full_resids[gap_start + k] = (left + k + 1, ' ')
                else:
                    # Not enough room — number backwards from right
                    for k in range(gap_len):
                        full_resids[gap_start + k] = (right - gap_len + k, ' ')
            else:
                # C-terminal gap (nothing to the right): number from left
                for k in range(gap_len):
                    full_resids[gap_start + k] = (left + k + 1, ' ')

        # Modeller numbers all chains continuously: 1..N_total
        offset = sum(len(per_chain_masks[j]) for j in range(ci))
        for pos_in_chain in range(len(mask)):
            model_resnum = offset + pos_in_chain + 1
            mapping[(chain, model_resnum)] = full_resids[pos_in_chain]

    return mapping


def restore_chain_ids_and_read(model_pdb_path, all_chains):
    """Read Modeller output PDB, remap chain IDs back to original.

    Modeller assigns A, B, C, ... to all chains (protein + non-protein).
    """
    modeller_chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    chain_map = {}
    for i, orig_chain in enumerate(all_chains):
        if i < len(modeller_chains):
            chain_map[modeller_chains[i]] = orig_chain

    output_lines = []
    with open(model_pdb_path) as f:
        content = f.read()

    # Modeller can write TER records without trailing newline,
    # causing them to be joined with the next line. Split on known
    # record boundaries to handle this.
    lines = content.split('\n')
    expanded = []
    for line in lines:
        # Check if a line contains a TER followed by ATOM/HETATM (no newline between)
        while True:
            found = False
            for prefix in ("ATOM  ", "HETATM"):
                idx = line.find(prefix, 1)  # skip pos 0 (the line itself might start with it)
                if idx > 0 and line[:idx].startswith("TER"):
                    expanded.append(line[:idx])
                    line = line[idx:]
                    found = True
                    break
            if not found:
                break
        if line:
            expanded.append(line)

    for line in expanded:
        if not line.endswith('\n'):
            line = line + '\n'
        if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
            if len(line) > 21:
                mod_chain = line[21]
                if mod_chain in chain_map:
                    line = line[:21] + chain_map[mod_chain] + line[22:]
            output_lines.append(line)
        elif line.startswith("TER"):
            if len(line) > 21:
                mod_chain = line[21]
                if mod_chain in chain_map:
                    line = line[:21] + chain_map[mod_chain] + line[22:]
            output_lines.append(line)
        elif line.startswith("END"):
            output_lines.append(line)

    return output_lines


def renumber_model_output(result_lines, resnum_mapping):
    """Apply residue renumbering to Modeller output lines.

    resnum_mapping values are (resSeq, iCode) tuples.
    Only renumbers residues that appear in the mapping (protein chains).
    Non-protein chains are left with their numbering.
    """
    output = []
    for line in result_lines:
        if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
            chain = line[21]
            resseq = int(line[22:26].strip())
            key = (chain, resseq)
            if key in resnum_mapping:
                new_resseq, new_icode = resnum_mapping[key]
                line = _set_resid(line, new_resseq, new_icode)
            output.append(line)
        elif line.startswith("TER"):
            if len(line) > 26:
                chain = line[21]
                seq_str = line[22:26].strip()
                if seq_str and seq_str.isdigit():
                    key = (chain, int(seq_str))
                    if key in resnum_mapping:
                        new_resseq, new_icode = resnum_mapping[key]
                        line = _set_resid(line, new_resseq, new_icode)
            output.append(line)
        else:
            output.append(line)
    return output


def patch_missing_hetatm(result_lines, original_lines, all_chains, protein_chains, verbose=False):
    """Patch back HETATM atoms that Modeller dropped (e.g. linkage atoms).

    Modeller's BLK residue handling can drop atoms at protein-glycan linkage
    points (e.g. NAG C1 that forms the N-glycosidic bond to ASN). This function
    compares output vs original non-protein residues and inserts any missing
    atoms from the original PDB.

    Matching is done by (chain, resname, position-in-chain) since Modeller
    renumbers residues.
    """
    protein_set = set(protein_chains)

    def get_nonprotein_residues(lines, chain_set=None):
        """Get ordered list of (chain, resnum, resname) and atoms per residue."""
        residues = {}  # (chain, resnum) -> {resname, atoms: [lines]}
        order = []     # ordered unique (chain, resnum)
        for line in lines:
            if not line.startswith("HETATM"):
                continue
            ch = line[21]
            if chain_set and ch not in chain_set:
                continue
            resnum = line[22:26].strip()
            resname = line[17:20].strip()
            key = (ch, resnum)
            if key not in residues:
                residues[key] = {"resname": resname, "atoms": {}, "lines": {}}
                order.append(key)
            atomname = line[12:16].strip()
            residues[key]["atoms"][atomname] = line
            residues[key]["lines"][atomname] = line
        return residues, order

    nonprotein_set = set(all_chains) - protein_set
    orig_res, orig_order = get_nonprotein_residues(original_lines, nonprotein_set)
    out_res, out_order = get_nonprotein_residues(result_lines, nonprotein_set)

    # Build position-based mapping: (chain, resname, nth_occurrence) -> residue key
    def build_pos_map(order, residues):
        pos_map = {}
        counter = {}
        for key in order:
            ch, _ = key
            resname = residues[key]["resname"]
            ck = (ch, resname)
            counter.setdefault(ck, 0)
            counter[ck] += 1
            pos_map[(ch, resname, counter[ck])] = key
        return pos_map

    orig_pos = build_pos_map(orig_order, orig_res)
    out_pos = build_pos_map(out_order, out_res)

    # Find missing atoms and build patches
    patches = {}  # out_resnum_key -> list of original lines to insert
    n_patched = 0
    for pos_key, orig_key in orig_pos.items():
        if pos_key not in out_pos:
            continue
        out_key = out_pos[pos_key]
        orig_atoms = set(orig_res[orig_key]["atoms"].keys())
        out_atoms = set(out_res[out_key]["atoms"].keys())
        missing = orig_atoms - out_atoms
        if missing:
            out_chain, out_resnum = out_key
            for atomname in sorted(missing):
                orig_line = orig_res[orig_key]["atoms"][atomname]
                # Update chain and residue number to match output
                patched = orig_line[:21] + out_chain + f"{int(out_resnum):4d}" + " " + orig_line[27:]
                patches.setdefault(out_key, []).append(patched)
                n_patched += 1

    if n_patched == 0:
        return result_lines

    if verbose:
        print(f"Patched {n_patched} missing HETATM atom(s) from original")

    # Insert patched atoms after the last atom of each residue
    output = []
    for line in result_lines:
        output.append(line)
        if line.startswith("HETATM"):
            ch = line[21]
            resnum = line[22:26].strip()
            key = (ch, resnum)
            if key in patches:
                # Check if next line is same residue; if not, insert here
                pass  # We'll handle insertion differently

    # Simpler approach: collect all lines per residue, insert patches
    output = []
    prev_key = None
    for line in result_lines:
        if line.startswith("HETATM"):
            ch = line[21]
            resnum = line[22:26].strip()
            cur_key = (ch, resnum)
            # If we moved to a new residue, insert patches for the previous one
            if prev_key and prev_key != cur_key and prev_key in patches:
                output.extend(patches.pop(prev_key))
            prev_key = cur_key
        else:
            # Non-HETATM line — flush patches for previous residue
            if prev_key and prev_key in patches:
                output.extend(patches.pop(prev_key))
            prev_key = None
        output.append(line)

    # Flush any remaining patches
    if prev_key and prev_key in patches:
        # Insert before END if present
        if output and output[-1].startswith("END"):
            end_line = output.pop()
            output.extend(patches.pop(prev_key))
            output.append(end_line)
        else:
            output.extend(patches.pop(prev_key))

    return output


def restore_conect_records(result_lines, original_lines, verbose=False):
    """Restore CONECT records from original PDB with remapped atom serials.

    Modeller strips all CONECT records. This function rebuilds them using:
    - ATOM records (protein): matched by (chain, resname, resnum, atomname)
      directly, since renumber_model_output already restored original numbering
    - HETATM records (non-protein): matched by (chain, resname, nth_occurrence,
      atomname), since Modeller renumbers non-protein residues
    """
    # Build positional map for HETATM residues only (non-protein)
    def _build_hetatm_position_map(lines):
        """Map (chain, resnum) -> (chain, resname, nth_occurrence) for HETATM."""
        counter = {}
        res_to_pos = {}
        seen = set()
        for line in lines:
            if not line.startswith("HETATM"):
                continue
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            key = (chain, resnum)
            if key not in seen:
                seen.add(key)
                ck = (chain, resname)
                counter.setdefault(ck, 0)
                counter[ck] += 1
                res_to_pos[key] = (chain, resname, counter[ck])
        return res_to_pos

    orig_het_pos = _build_hetatm_position_map(original_lines)
    model_het_pos = _build_hetatm_position_map(result_lines)

    # Build original serial -> lookup key
    # ATOM: key = ('atom', chain, resname, resnum, atomname)
    # HETATM: key = ('het', chain, resname, nth, atomname)
    orig_serial_to_key = {}
    for line in original_lines:
        if line.startswith("ATOM  "):
            serial = int(line[6:11].strip())
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            atomname = line[12:16].strip()
            orig_serial_to_key[serial] = ('atom', chain, resname, resnum, atomname)
        elif line.startswith("HETATM"):
            serial = int(line[6:11].strip())
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            atomname = line[12:16].strip()
            pos = orig_het_pos.get((chain, resnum))
            if pos:
                orig_serial_to_key[serial] = ('het', pos[0], pos[1], pos[2], atomname)

    # Build lookup key -> new serial in model output
    key_to_new_serial = {}
    for line in result_lines:
        if line.startswith("ATOM  "):
            serial = int(line[6:11].strip())
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            atomname = line[12:16].strip()
            key = ('atom', chain, resname, resnum, atomname)
            key_to_new_serial[key] = serial
        elif line.startswith("HETATM"):
            serial = int(line[6:11].strip())
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            atomname = line[12:16].strip()
            pos = model_het_pos.get((chain, resnum))
            if pos:
                key = ('het', pos[0], pos[1], pos[2], atomname)
                key_to_new_serial[key] = serial

    # Parse CONECT records from original
    conect_records = []
    for line in original_lines:
        if not line.startswith("CONECT"):
            continue
        serials = []
        text = line.rstrip()
        for i in range(6, len(text), 5):
            field = text[i:i + 5].strip()
            if field and field.isdigit():
                serials.append(int(field))
        if len(serials) >= 2:
            conect_records.append(serials)

    if not conect_records:
        return result_lines

    # Remap serials
    new_conect_lines = []
    n_dropped = 0
    for serials in conect_records:
        new_serials = []
        skip = False
        for s in serials:
            key = orig_serial_to_key.get(s)
            if key is None:
                skip = True
                break
            new_s = key_to_new_serial.get(key)
            if new_s is None:
                skip = True
                break
            new_serials.append(new_s)
        if skip:
            n_dropped += 1
            continue
        line = "CONECT"
        for ns in new_serials:
            line += f"{ns:5d}"
        line = line.ljust(80) + "\n"
        new_conect_lines.append(line)

    if verbose:
        print(f"Restored {len(new_conect_lines)} CONECT records"
              + (f" ({n_dropped} dropped — atoms not found)" if n_dropped else ""))

    if not new_conect_lines:
        return result_lines

    # Insert CONECT before END line
    output = []
    end_line = None
    for line in result_lines:
        if line.startswith("END"):
            end_line = line
        else:
            output.append(line)
    output.extend(new_conect_lines)
    if end_line:
        output.append(end_line)
    else:
        output.append("END\n")

    return output


def build_model_dat(result_lines, per_chain_masks, all_chains, protein_chains, resnum_mapping):
    """Build .dat file recording gap-filled (rebuilt) atoms.

    Uses per_chain_masks (True=template, False=gap) and resnum_mapping
    to identify which residues were rebuilt by Modeller. All atoms in
    rebuilt residues are recorded as "added".
    """
    protein_set = set(protein_chains)

    # Collect (chain, resSeq, iCode) for gap-filled positions
    gap_resids = set()
    for ci, chain in enumerate(all_chains):
        if chain not in protein_set:
            continue
        if ci >= len(per_chain_masks):
            continue
        mask = per_chain_masks[ci]
        offset = sum(len(per_chain_masks[j]) for j in range(ci))
        for pos, is_template in enumerate(mask):
            if not is_template:
                model_resnum = offset + pos + 1
                resid = resnum_mapping.get((chain, model_resnum))
                if resid is not None:
                    gap_resids.add((chain, resid[0], resid[1].strip()))

    # Scan output PDB lines for atoms in gap residues
    added_atoms = []
    residue_summary = {}
    for line in result_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21]
        seq_str = line[22:26].strip()
        if not seq_str or not seq_str.lstrip('-').isdigit():
            continue
        resnum = int(seq_str)
        icode = line[26].strip() if len(line) > 26 else ''
        if (chain, resnum, icode) in gap_resids:
            resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 77 else ''
            resid_str = str(resnum)
            added_atoms.append({
                "chain": chain,
                "resid": resid_str,
                "icode": icode,
                "resname": resname,
                "atom": atom_name,
                "element": element,
            })
            rkey = f"{chain}/{resname}{resid_str}"
            if rkey not in residue_summary:
                residue_summary[rkey] = {"heavy": 0, "hydrogen": 0}
            if element == 'H':
                residue_summary[rkey]["hydrogen"] += 1
            else:
                residue_summary[rkey]["heavy"] += 1

    return {
        "description": "Modeller gap-fill data. Rebuilt residue atoms get weak/no restraints "
                       "during minimization.",
        "total_added": len(added_atoms),
        "residue_summary": residue_summary,
        "added_atoms": added_atoms,
    }


def remove_water_lines(lines):
    """Remove water molecules from PDB lines."""
    output = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM", "ANISOU")):
            resname = line[17:20].strip()
            if resname in WATER_RESNAMES:
                continue
        elif line.startswith("TER") and len(line) > 20:
            resname = line[17:20].strip()
            if resname in WATER_RESNAMES:
                continue
        output.append(line)
    return output


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        print(f"File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        output_path = Path(args.output).resolve()
    else:
        output_path = input_path.with_stem(input_path.stem + "_model")

    with open(input_path) as f:
        lines = f.readlines()

    # Determine protein chains (those with SEQRES)
    seqres = parse_seqres(lines)
    if args.fasta:
        fasta_seqs = parse_fasta(args.fasta)
        protein_chains = get_chain_order(lines, seqres_only=set(seqres.keys())) if seqres else get_chain_order(lines)
        if len(fasta_seqs) != len(protein_chains):
            print(f"Error: FASTA has {len(fasta_seqs)} sequences but PDB has "
                  f"{len(protein_chains)} protein chains ({', '.join(protein_chains)})",
                  file=sys.stderr)
            sys.exit(1)
        protein_seq_map = {ch: seq for ch, (_, seq) in zip(protein_chains, fasta_seqs)}
    else:
        if not seqres:
            print("No SEQRES records found — no gaps to model, copying input")
            import shutil as _sh
            _sh.copy2(str(input_path), str(output_path))
            print(f"Wrote {output_path}")
            return
        protein_chains = get_chain_order(lines, seqres_only=set(seqres.keys()))
        protein_seq_map = {}
        for ch in protein_chains:
            protein_seq_map[ch] = ''.join(AA3TO1.get(r, 'X') for r in seqres[ch])

    if not protein_chains:
        print("No modelable chains found.", file=sys.stderr)
        sys.exit(1)

    # ALL chains in order of appearance (protein + non-protein)
    all_chains = get_chain_order(lines)
    nonprotein = [c for c in all_chains if c not in set(protein_chains)]

    if args.verbose:
        print(f"Protein chains: {', '.join(protein_chains)}")
        if nonprotein:
            print(f"Non-protein chains (included via hetatm): {', '.join(nonprotein)}")

    if args.no_terminal:
        print("Trimming N/C terminal residues (--no-terminal)")
        protein_seq_map = trim_terminal_gaps(protein_seq_map, lines, args.verbose)

    n_protein_res = sum(len(protein_seq_map[ch]) for ch in protein_chains)
    print(f"Target sequence: {n_protein_res} protein residues across {len(protein_chains)} chain(s)")

    # Work in temp directory
    workdir = tempfile.mkdtemp(prefix='dvbfixer_model_')
    orig_dir = os.getcwd()

    try:
        os.chdir(workdir)

        # Write FULL PDB (all chains) for Modeller
        input_pdb = workdir + '/' + input_path.name
        with open(input_pdb, 'w') as f:
            f.writelines(lines)

        # Run Modeller — target is built inside using Modeller's own template
        # reading, so non-protein '.' counts are guaranteed to match
        best_model, aln_path = run_modeller(
            Path(input_pdb), protein_chains, protein_seq_map, all_chains, args
        )

        # Restore original chain IDs (Modeller reassigns A,B,C,... to all chains)
        result_lines = restore_chain_ids_and_read(best_model, all_chains)

        # Restore residue numbering for protein chains only
        _, per_chain_masks = get_template_mask(aln_path)
        resnum_mapping = build_resnum_mapping(
            per_chain_masks, all_chains, protein_chains, lines
        )
        result_lines = renumber_model_output(result_lines, resnum_mapping)

        # Patch back any HETATM atoms Modeller dropped (e.g. NAG C1 linkage atom)
        result_lines = patch_missing_hetatm(
            result_lines, lines, all_chains, protein_chains, args.verbose
        )

        # Restore CONECT records from original PDB with remapped atom serials
        result_lines = restore_conect_records(result_lines, lines, args.verbose)

        if args.verbose:
            n_gaps_filled = sum(
                not m for ci, cm in enumerate(per_chain_masks)
                for m in cm
                if ci < len(all_chains) and all_chains[ci] in set(protein_chains)
            )
            print(f"Restored numbering ({n_gaps_filled} gap positions assigned)")

        if not args.keep_water:
            result_lines = remove_water_lines(result_lines)

        # Final pass: split any remaining concatenated TER+ATOM/HETATM lines
        # that may have survived post-processing, and ensure all lines end with \n
        sanitized = []
        for line in result_lines:
            while True:
                found = False
                for prefix in ("ATOM  ", "HETATM"):
                    idx = line.find(prefix, 1)
                    if idx > 0 and line[:idx].rstrip().startswith("TER"):
                        ter = line[:idx].rstrip() + "\n"
                        sanitized.append(ter)
                        line = line[idx:]
                        found = True
                        break
                if not found:
                    break
            if not line.endswith("\n"):
                line = line + "\n"
            sanitized.append(line)
        result_lines = sanitized

        with open(str(output_path), 'w') as f:
            f.writelines(result_lines)
        print(f"Wrote {output_path}")

        # Write .dat file only if there were actual gaps
        n_gaps = sum(
            not m for ci, chain in enumerate(all_chains)
            if chain in set(protein_chains) and ci < len(per_chain_masks)
            for m in per_chain_masks[ci]
        )
        if n_gaps > 0:
            dat = build_model_dat(
                result_lines, per_chain_masks, all_chains, protein_chains, resnum_mapping
            )
            dat_path = output_path.with_suffix('.dat')
            with open(dat_path, 'w') as f:
                json.dump(dat, f, indent=2)
            print(f"Saved restraint data: {dat_path} ({dat['total_added']} atoms in rebuilt regions)")

    finally:
        os.chdir(orig_dir)
        if args.keep_workdir:
            print(f"Workdir kept: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)
