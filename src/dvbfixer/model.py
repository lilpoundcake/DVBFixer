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
import re
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


def _explain_modeller_error(exc, protein_chains, protein_seq_map):
    """Translate a Modeller exception into a useful diagnostic.

    Modeller error messages typically include either a residue index
    (1-based position in the target sequence) or a (chain, resseq)
    reference. They look like:

      "Number of residues in the alignment and pdb files are different"
      "Residue type 'X' too long: 'BLK'"
      "Heavy atom at index N has zero coordinates"
      "Atom 'XX' not found in residue NNN:C"

    This helper looks for `(\\d+):([A-Z])` or `residue (\\d+)` style
    references in the message and resolves them to the actual
    (chain, resseq, resname) using `protein_seq_map`. Returns a
    multi-line string with the resolved residue + neighbour context, or
    None if no recognizable reference was found.

    `protein_chains`: list of chain IDs in the target order.
    `protein_seq_map`: dict {chain_id: [(resseq, icode, resname), ...]}.
    """
    import re

    msg = str(exc)
    lines = []

    # Recognise well-known Modeller error classes and add a plain-language
    # explanation that points at the likely cause + remediation. The raw
    # Modeller message is still printed unchanged before this diagnostic.
    if "No aligned template residues for BLK residue" in msg:
        lines.append(
            "Cause: Modeller's PIR alignment has at least one '.' (BLK / "
            "non-protein residue placeholder) in the target sequence that "
            "is not paired with a matching '.' in the template. This "
            "usually means a chain ID appears in two separate file-order "
            "blocks (e.g. protein ATOMs early, glycan HETATMs later, "
            "separated by other chains), which makes Modeller emit one "
            "extra template chain segment that the target sequence builder "
            "does not produce. dvbfixer reorders chains automatically — if "
            "you still hit this, check that no chain ID appears in two "
            "disjoint segments of the input PDB."
        )
    elif "Sequence difference between alignment and  pdb" in msg or \
         "alignment sequence must match that from the atom file" in msg.lower():
        lines.append(
            "Cause: The target sequence built from SEQRES (or --fasta) does "
            "not line up with the template residues Modeller read from the "
            "ATOM records. Common reasons: a chain ID split across multiple "
            "file-order blocks (protein + late HETATM); a SEQRES that omits "
            "or doubles a non-protein chain; or non-standard residue names "
            "Modeller does not recognise."
        )
    elif "Residue type" in msg and "too long" in msg and "BLK" in msg:
        lines.append(
            "Cause: A non-standard 3-letter residue name (HETATM) is in the "
            "input but Modeller cannot map it to a known type. Check for "
            "typos in HETATM resname columns 18-20 or convert exotic ligand "
            "names to a known het code before running model."
        )

    # protein_seq_map values are 1-letter AA sequence STRINGS (one letter
    # per protein residue in order). The exact PDB resseq/resname for each
    # position is in the original PDB lines, not in protein_seq_map.
    # Best-effort: report position-in-sequence + chain ID.

    # Try `<resseq>:<chain>` Modeller convention first — this is already
    # the most useful format because <resseq>:<chain> directly identifies
    # the residue in the user's PDB.
    for m in re.finditer(r'(\d+)\s*:\s*([A-Za-z0-9])', msg):
        resseq, chain = int(m.group(1)), m.group(2)
        if chain in protein_seq_map:
            lines.append(
                f"Modeller error refers to residue {resseq}:{chain} "
                f"(chain {chain}, resseq {resseq} — check this position "
                f"in your input PDB)"
            )

    # Try `residue <N>` where N is a 1-based index into target sequence
    for m in re.finditer(r'residue\s+(\d+)', msg, re.IGNORECASE):
        idx = int(m.group(1))
        # Walk the per-chain protein sequences to find which chain idx
        # falls into, and the position within that chain.
        cursor = 0
        for ch in protein_chains:
            seq_len = len(protein_seq_map.get(ch, ''))
            if cursor < idx <= cursor + seq_len:
                pos_in_chain = idx - cursor   # 1-based
                letter = protein_seq_map[ch][pos_in_chain - 1]
                lines.append(
                    f"Modeller error refers to residue index {idx} → "
                    f"chain {ch}, position {pos_in_chain} of {seq_len} "
                    f"(AA letter '{letter}' in target sequence)"
                )
                break
            cursor += seq_len

    if not lines:
        # Generic fallback: list per-chain residue counts so the user
        # can at least see what chains were being modeled.
        for ch in protein_chains:
            seq = protein_seq_map.get(ch, '')
            if not seq:
                continue
            lines.append(f"chain {ch}: {len(seq)} protein residue(s) in target")
        if lines:
            lines.insert(0, "No specific residue identified in error message. "
                            "Chain sizes being modeled:")

    return '\n'.join(lines) if lines else None


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
        "--fasta", help="FASTA file with complete sequence(s). Headers must encode "
        "chain IDs: '>chain_X', '>PDBID_X', or '>X'. Mapping is by chain ID, "
        "not file order. Use instead of SEQRES."
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


_CHAIN_PATTERNS = [
    re.compile(r'^chain[_\s]*([A-Za-z0-9])$', re.IGNORECASE),
    re.compile(r'_([A-Za-z0-9])$'),
    re.compile(r'^([A-Za-z0-9])$'),
]


def _chain_id_from_header(header):
    """Extract chain ID from a FASTA header. Returns None on failure."""
    h = header.strip()
    for pat in _CHAIN_PATTERNS:
        m = pat.search(h)
        if m:
            return m.group(1).upper()
    return None


def parse_fasta(fasta_path):
    """Parse a FASTA file, return dict {chain_id: sequence}.

    Headers are parsed for a chain ID. Supported header formats:
      >chain_X / >chainX  (explicit prefix, case-insensitive)
      >XXXX_X             (PDB-style, e.g. '1abc_A')
      >X                  (single-character header)

    Raises ValueError if a chain ID cannot be determined or if duplicate
    chain IDs are found.
    """
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_path, "fasta"))
    result = {}
    for r in records:
        ch = _chain_id_from_header(r.id)
        if ch is None:
            raise ValueError(
                f"Could not extract chain ID from FASTA header '{r.id}'. "
                "Use '>chain_X', '>PDBID_X', or '>X' format."
            )
        if ch in result:
            raise ValueError(f"Duplicate chain ID '{ch}' in FASTA file.")
        result[ch] = str(r.seq)
    return result


def _reorder_chains_for_modeller(lines):
    """Group ATOM/HETATM lines so each chain ID appears as one contiguous block.

    Modeller starts a new PIR chain segment on every chain-ID change (and on
    every TER). When a chain ID is split into multiple file-order blocks by
    interleaved other chains (e.g. chain A protein → chains B/C/D/E glycans →
    chain A HETATM glycan), Modeller emits multiple template segments for the
    same ID. That breaks downstream code that pairs Modeller's segments
    one-for-one with the unique chain IDs returned by `get_chain_order`, and
    can manifest as 'No aligned template residues for BLK residue' from
    align2d.

    This reorderer:
    - preserves all non-ATOM/HETATM/TER/CONECT lines (REMARK, SEQRES, LINK …)
    - groups ATOM/HETATM by chain ID in first-appearance order
    - drops original TER and CONECT records (Modeller infers HETATM bonds via
      `env.io.hetatm=True`; CONECT is restored downstream from atom identity)
    - emits a single TER after each chain block
    - preserves atom serials (CONECT is restored by identity, not serial)
    """
    chain_lines = {}
    chain_order = []
    headers = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            ch = line[21]
            if ch not in chain_lines:
                chain_lines[ch] = []
                chain_order.append(ch)
            chain_lines[ch].append(line)
        elif line.startswith(("TER", "CONECT")):
            continue
        elif line.startswith("END"):
            continue
        else:
            headers.append(line)

    out = list(headers)
    for ch in chain_order:
        out.extend(chain_lines[ch])
        last = chain_lines[ch][-1]
        last_resseq = last[22:26]
        last_resname = last[17:20]
        last_icode = last[26] if len(last) > 26 else ' '
        ter = f"TER   {0:>5d}      {last_resname} {ch}{last_resseq}{last_icode}\n"
        out.append(ter)
    out.append("END\n")
    return out


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

        # Strip leading/trailing dots — these are HETATM (BLK) slots that
        # must be aligned 1:1 (target & template both got the same `.` count
        # when target was built from Modeller's template). Treating them as
        # protein letters or as suffix/prefix gaps in the matcher below
        # collapses them to `-`, which then drops the HETATM from the
        # template and crashes Modeller with 'No aligned template residues
        # for BLK residue'.
        def _strip_outer_dots(s):
            n_lead = 0
            while n_lead < len(s) and s[n_lead] == '.':
                n_lead += 1
            n_trail = 0
            while n_trail < len(s) - n_lead and s[-(n_trail + 1)] == '.':
                n_trail += 1
            middle = s[n_lead:len(s) - n_trail] if n_trail else s[n_lead:]
            return n_lead, middle, n_trail

        tpl_lead, tpl_pure, tpl_trail = _strip_outer_dots(tpl_pure)
        tgt_lead, tgt_pure, tgt_trail = _strip_outer_dots(tgt_pure)

        # Only fix protein residues (letters); any `.` left inside `tpl_pure`
        # would be unexpected — bail out conservatively in that case.
        tpl_protein = ''.join(c for c in tpl_pure if c != '.')
        if tpl_protein != tpl_pure:
            fixed_tpl_chains.append(raw_tpl)
            fixed_tgt_chains.append(raw_tgt)
            continue
        tgt_protein = tgt_pure

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
            new_tpl = '-' * n_prefix + tpl_protein + '-' * n_suffix
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

        # Re-attach leading/trailing dots: target keeps its dot count; for
        # template we pair as many as it had, padding any shortfall with `-`
        # (extra target slots) or absorbing extras as gaps in target (extra
        # template slots — shouldn't normally happen).
        def _pair_dots(n_tpl, n_tgt):
            if n_tpl == n_tgt:
                return '.' * n_tpl, '.' * n_tgt
            if n_tpl < n_tgt:
                return '.' * n_tpl + '-' * (n_tgt - n_tpl), '.' * n_tgt
            return '.' * n_tpl, '.' * n_tgt + '-' * (n_tpl - n_tgt)

        lead_tpl, lead_tgt = _pair_dots(tpl_lead, tgt_lead)
        trail_tpl, trail_tgt = _pair_dots(tpl_trail, tgt_trail)
        new_tpl = lead_tpl + new_tpl + trail_tpl
        new_tgt = lead_tgt + new_tgt + trail_tgt

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
            # Preserve HETATM residues attached to a protein chain by
            # counting '.' (BLK) entries in Modeller's template sequence
            # for this chain and re-appending the same count to the SEQRES-
            # derived target. Without this, glycans / ligands / ions
            # attached to a protein chain (e.g. N-linked NAG covalently
            # bonded to ASN.ND2 on the same chain) get dropped by Modeller
            # during loop modeling because the target has no slot for them.
            n_dots = tpl_chain.count('.')
            seq = protein_seq_map[ch]
            if n_dots:
                seq = seq + '.' * n_dots
            target_chain_seqs.append(seq)
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

    try:
        a.make()
    except Exception as e:
        diag = _explain_modeller_error(e, protein_chains, protein_seq_map)
        print(f"\nERROR: Modeller failed during model building:\n  {e}\n",
              file=sys.stderr)
        if diag:
            for line in diag.split('\n'):
                print(f"  {line}", file=sys.stderr)
        raise

    loop_models = [x for x in a.loop.outputs if x['failure'] is None]
    if not loop_models:
        init_models = [x for x in a.outputs if x['failure'] is None]
        if not init_models:
            # Some/all initial models failed too — Modeller stores the per-
            # model failure string. Print the failures with context so the
            # user can see which residue/chain triggered the problem.
            failed_init = [x for x in a.outputs if x.get('failure') is not None]
            failed_loop = [x for x in a.loop.outputs if x.get('failure') is not None]
            print("\nERROR: All Modeller models failed.\n", file=sys.stderr)
            for i, x in enumerate(failed_init):
                print(f"  initial model {i+1} failure: {x['failure']}",
                      file=sys.stderr)
            for i, x in enumerate(failed_loop):
                print(f"  loop model {i+1} failure: {x['failure']}",
                      file=sys.stderr)
            # Best-effort: parse the first failure for residue context.
            first_fail = (failed_init or failed_loop)[0] if (failed_init or failed_loop) else None
            if first_fail:
                diag = _explain_modeller_error(
                    Exception(str(first_fail['failure'])),
                    protein_chains, protein_seq_map)
                if diag:
                    print("", file=sys.stderr)
                    for line in diag.split('\n'):
                        print(f"  {line}", file=sys.stderr)
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


def _get_original_resletters_per_chain(original_lines, chain_order):
    """Per protein chain: ordered list of (resseq, icode, AA-letter) for ATOM records."""
    chain_set = set(chain_order)
    result = {ch: [] for ch in chain_order}
    seen = {ch: set() for ch in chain_order}
    for line in original_lines:
        if not line.startswith("ATOM  "):
            continue
        ch = line[21]
        if ch not in chain_set:
            continue
        resseq = int(line[22:26].strip())
        icode = line[26]
        key = (resseq, icode)
        if key not in seen[ch]:
            resname = line[17:20].strip()
            letter = AA3TO1.get(resname, 'X')
            seen[ch].add(key)
            result[ch].append((resseq, icode, letter))
    return result


def _find_seqres_offset_by_resseq(orig_residues, seqres_seq, tolerance=0.1):
    """Find offset K such that SEQRES[K + (resseq - first_resseq)] == atom letter
    for as many atoms as possible.

    K is the number of N-terminal SEQRES residues absent from ATOM (signal
    peptide etc.). This formula assumes input resseq is monotonic and
    contiguous within each stretch — i.e. the resseq jumps in the input
    correspond exactly to missing SEQRES positions. That makes it the
    "right" alignment for the user's perspective (gap-fill resseqs match
    the gap in the input).

    Picks the K with the highest match count. Accepts K when
    matches >= (1 - tolerance) * len(atoms), allowing a few mutations.

    Returns K (int >= 0) or None if no usable offset exists (icodes
    present, resseq exceeds SEQRES, or too many mismatches).
    """
    if not orig_residues:
        return 0
    if any(icode != ' ' for _, icode, _ in orig_residues):
        return None
    first_resseq = orig_residues[0][0]
    L = len(seqres_seq)
    max_offset = orig_residues[-1][0] - first_resseq
    if max_offset >= L:
        return None
    best_K = None
    best_matches = -1
    threshold = (1.0 - tolerance) * len(orig_residues)
    for K in range(L - max_offset):
        matches = 0
        for resseq, _, letter in orig_residues:
            pos = K + (resseq - first_resseq)
            if seqres_seq[pos] == letter or letter == 'X':
                matches += 1
        if matches > best_matches:
            best_matches = matches
            best_K = K
        if matches == len(orig_residues):
            return K  # perfect match — short-circuit
    if best_matches >= threshold:
        return best_K
    return None


def _align_atoms_to_seqres(orig_residues, seqres_seq):
    """Align ATOM letters to SEQRES letters via semi-global Needleman-Wunsch.

    The result is a list of length `len(orig_residues)` giving the SEQRES
    position (0-indexed) each ATOM residue maps to, or `None` for that entry
    if the residue could not be aligned. Returns `None` for the whole call
    when no usable alignment exists (e.g. ATOMs longer than SEQRES, no
    letter information).

    Semi-global = free end gaps on the SEQRES side, no skipping of ATOM
    letters. Affine gaps with a heavy gap-open penalty so consecutive
    gaps stay clumped (align2d's blosum scoring can split them, which is
    exactly the FcgRI failure mode this is bypassing).

    Robustness vs the previous `_find_seqres_offset`:
    - tolerates point mutations (mismatch score absorbed, alignment still
      placed in the only sensible spot)
    - tolerates `X` (unknown) residues with zero cost — they don't anchor
      or punish
    - handles N-terminal AND C-terminal SEQRES extras automatically via
      free end gaps
    - falls back gracefully when SEQRES is shorter than ATOMs
    """
    Q = len(orig_residues)
    R = len(seqres_seq)
    if Q == 0:
        return []
    if R == 0 or Q > R:
        return None

    query = [r[2] for r in orig_residues]

    MATCH, MISMATCH, X_NEUTRAL = 2, -1, 0
    GAP_OPEN, GAP_EXTEND = -10, -1
    NEG = float('-inf')

    # Three matrices for affine gaps:
    #   M[i][j]  = best score ending with query[i-1] aligned to seqres[j-1]
    #   IX[i][j] = best score ending with a SEQRES gap (skip seqres[j-1]; not used since gap-in-seqres = skip-query, disallowed)
    #   IY[i][j] = best score ending with a query gap (skip seqres[j-1])
    # Because we forbid skipping query letters, IX is dropped. Only IY is needed.
    M = [[NEG] * (R + 1) for _ in range(Q + 1)]
    IY = [[NEG] * (R + 1) for _ in range(Q + 1)]
    M[0][0] = 0
    # free end gaps on SEQRES side at the START: M[0][j] = 0 (any j)
    for j in range(R + 1):
        M[0][j] = 0
        IY[0][j] = NEG

    bt_M = [[None] * (R + 1) for _ in range(Q + 1)]
    bt_IY = [[None] * (R + 1) for _ in range(Q + 1)]

    for i in range(1, Q + 1):
        for j in range(1, R + 1):
            q, r = query[i - 1], seqres_seq[j - 1]
            if q == 'X' or r == 'X':
                sub = X_NEUTRAL
            elif q == r:
                sub = MATCH
            else:
                sub = MISMATCH
            # M[i][j]: match/mismatch from diag
            from_M = M[i - 1][j - 1] + sub
            from_IY = IY[i - 1][j - 1] + sub
            if from_M >= from_IY:
                M[i][j] = from_M
                bt_M[i][j] = 'M'
            else:
                M[i][j] = from_IY
                bt_M[i][j] = 'IY'
            # IY[i][j]: gap in query (skip seqres[j-1])
            open_iy = M[i][j - 1] + GAP_OPEN
            ext_iy = IY[i][j - 1] + GAP_EXTEND
            if open_iy >= ext_iy:
                IY[i][j] = open_iy
                bt_IY[i][j] = 'M'
            else:
                IY[i][j] = ext_iy
                bt_IY[i][j] = 'IY'

    # Free end gaps on SEQRES side at the END: best ending = max of M[Q][j] for j in 1..R
    best_j = 0
    best_score = NEG
    for j in range(1, R + 1):
        if M[Q][j] > best_score:
            best_score = M[Q][j]
            best_j = j
    if best_score == NEG:
        return None

    result = [None] * Q
    i, j = Q, best_j
    state = 'M'
    while i > 0 and j > 0:
        if state == 'M':
            result[i - 1] = j - 1
            prev = bt_M[i][j]
            i -= 1
            j -= 1
            state = prev
        else:  # IY
            prev = bt_IY[i][j]
            j -= 1
            state = prev
    if any(r is None for r in result):
        return None
    return result


def _interpolate_gaps(full_resids, region_len):
    """Fill `None` entries in `full_resids[:region_len]` by interpolating
    between flanking non-None entries.

    Handles:
    - internal gap (left and right templates present): place sequentially
      from `left + 1` if room, else number backward from `right`.
    - N-terminal gap (no left): number backward from the first template.
    - C-terminal gap (no right within region): number forward from the
      last template.
    - all-None region: number sequentially from 1.

    Operates in-place on `full_resids`. Entries beyond `region_len` are
    left alone (HETATM `.` slots are filled by the trailing sequential
    pass in `build_resnum_mapping`).
    """
    region_len = min(region_len, len(full_resids))
    i = 0
    while i < region_len:
        if full_resids[i] is not None:
            i += 1
            continue
        gap_start = i
        while i < region_len and full_resids[i] is None:
            i += 1
        gap_end = i
        gap_len = gap_end - gap_start

        left = None
        for j in range(gap_start - 1, -1, -1):
            if full_resids[j] is not None:
                left = full_resids[j][0]
                break
        right = None
        for j in range(gap_end, region_len):
            if full_resids[j] is not None:
                right = full_resids[j][0]
                break

        if left is not None and right is not None:
            available = right - left - 1
            if available >= gap_len:
                for k in range(gap_len):
                    full_resids[gap_start + k] = (left + k + 1, ' ')
            else:
                for k in range(gap_len):
                    full_resids[gap_start + k] = (right - gap_len + k, ' ')
        elif left is None and right is not None:
            for k in range(gap_len):
                full_resids[gap_start + k] = (right - gap_len + k, ' ')
        elif left is not None and right is None:
            for k in range(gap_len):
                full_resids[gap_start + k] = (left + k + 1, ' ')
        else:
            for k in range(gap_len):
                full_resids[gap_start + k] = (k + 1, ' ')


def build_resnum_mapping(per_chain_masks, all_chains, protein_chains, original_lines,
                          protein_seq_map=None):
    """Build mapping: (model_chain, model_resSeq) -> (original_resSeq, icode).

    Strategy (preferred): when input chain numbering is contiguous within the
    target sequence length and free of insertion codes, use a DETERMINISTIC
    resseq-based mapping: target-sequence position N (0-indexed) →
    resseq = first_resseq + N. Modeller produces one residue per target
    position, so this gives the correct resseq for both template and
    gap-filled positions without depending on align2d's gap placement
    (which can put a gap one position too early/late and then create cascade
    collisions in the dedup step, pushing modeled residues to the end of
    the chain).

    Fallback (when input has insertion codes, N-terminal SEQRES extras, or
    other non-contiguous numbering): walk align2d's mask one True position
    at a time, consuming `orig_rids` in order, and interpolate gap-filled
    positions between flanking originals (with N-term/C-term branching).

    Trailing positions beyond the protein region (HETATM `.` slots) are
    filled with sequential resseqs after the last protein residue.

    Only builds mapping for protein chains; non-protein chains are skipped.
    """
    protein_set = set(protein_chains)
    orig_resids = _get_original_resids_per_chain(original_lines, protein_chains)
    orig_letters = _get_original_resletters_per_chain(original_lines, protein_chains)
    mapping = {}
    protein_seq_map = protein_seq_map or {}

    for ci, chain in enumerate(all_chains):
        if ci >= len(per_chain_masks):
            continue
        mask = per_chain_masks[ci]

        if chain not in protein_set:
            # Non-protein chain — no renumbering needed
            continue

        orig_rids = orig_resids.get(chain, [])
        full_resids = [None] * len(mask)  # Each entry is (resSeq, iCode)

        seq_len = len(protein_seq_map.get(chain, ''))
        seqres_seq = protein_seq_map.get(chain, '')

        # Mapping strategy (tried in order):
        #
        # 1. resseq-based K-finder (preferred). Picks K = N-terminal-extras
        #    count and maps position N → resseq (first_resseq + N - K).
        #    Trusts the resseq jumps in the input PDB — when there's a
        #    jump from resseq 218 to 224 in chain A, that gap of 5 maps
        #    to SEQRES positions 198..202 (i.e. resseqs 219..223). This is
        #    what the user wants: gap-fill residues land inside the gap,
        #    not at the C-terminus.
        #
        # 2. Letter-only NW alignment (mutation-tolerant fallback). Used
        #    when K-finder bails (too many letter mismatches). NW handles
        #    point mutations but doesn't use resseq jumps, so it may place
        #    the gap one position off when the surrounding letters are
        #    ambiguous (e.g. repetitive sequence).
        #
        # 3. align2d's mask (last resort). Uses Modeller's alignment +
        #    flank interpolation. Susceptible to align2d's gap placement
        #    choices.
        seqres_index = None
        used = None
        if seq_len > 0 and orig_rids:
            K = _find_seqres_offset_by_resseq(orig_letters.get(chain, []), seqres_seq)
            if K is not None:
                used = 'K-finder'
                first_resseq = orig_rids[0][0]
                # Build seqres_index from K + (resseq - first_resseq)
                seqres_index = []
                for resseq, icode, _letter in orig_letters.get(chain, []):
                    if icode == ' ':
                        seqres_index.append(K + (resseq - first_resseq))
                    else:
                        seqres_index = None
                        break
            if seqres_index is None:
                # Try NW as fallback
                nw_result = _align_atoms_to_seqres(orig_letters.get(chain, []), seqres_seq)
                if nw_result is not None:
                    used = 'NW'
                    seqres_index = nw_result

        atom_only = orig_letters.get(chain, [])
        used_deterministic = (
            seqres_index is not None and len(seqres_index) == len(atom_only)
        )
        if used_deterministic:
            # Place each ATOM (protein-only) residue at its aligned SEQRES position.
            for k, (resseq, icode, _letter) in enumerate(atom_only):
                pos = seqres_index[k]
                if pos is None or pos < 0 or pos >= len(mask):
                    continue
                full_resids[pos] = (resseq, icode)
            _interpolate_gaps(full_resids, seq_len)

            # Place HETATM residues (attached glycans, ligands) at the
            # trailing `.` slots beyond the protein region, preserving
            # their ORIGINAL resseqs (so an N-linked NAG keeps the resseq
            # it had in the input PDB rather than getting renumbered to
            # the next sequential integer after the last protein residue).
            atom_keys = {(r, i) for r, i, _ in atom_only}
            hetatm_rids = [rid for rid in orig_rids if rid not in atom_keys]
            het_pos = seq_len
            for het_rid in hetatm_rids:
                while het_pos < len(mask) and full_resids[het_pos] is not None:
                    het_pos += 1
                if het_pos >= len(mask):
                    break
                full_resids[het_pos] = het_rid
                het_pos += 1
        else:
            # Fallback: align2d mask-based consumption + interpolation
            orig_idx = 0
            for i, is_template in enumerate(mask):
                if is_template and orig_idx < len(orig_rids):
                    full_resids[i] = orig_rids[orig_idx]
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

                # Find nearest non-None positions flanking the gap
                left = None
                for j in range(gap_start - 1, -1, -1):
                    if full_resids[j] is not None:
                        left = full_resids[j][0]
                        break
                right = None
                for j in range(gap_end, len(mask)):
                    if full_resids[j] is not None:
                        right = full_resids[j][0]
                        break

                if left is not None and right is not None:
                    # Internal gap: fit between left and right
                    available = right - left - 1
                    if available >= gap_len:
                        for k in range(gap_len):
                            full_resids[gap_start + k] = (left + k + 1, ' ')
                    else:
                        # Not enough room — number backwards from right
                        for k in range(gap_len):
                            full_resids[gap_start + k] = (right - gap_len + k, ' ')
                elif left is None and right is not None:
                    # N-terminal gap: extend backward from the first template residue
                    for k in range(gap_len):
                        full_resids[gap_start + k] = (right - gap_len + k, ' ')
                elif left is not None and right is None:
                    # C-terminal gap: extend forward from the last template residue
                    for k in range(gap_len):
                        full_resids[gap_start + k] = (left + k + 1, ' ')
                else:
                    # Whole chain is a gap — start at 1
                    for k in range(gap_len):
                        full_resids[gap_start + k] = (k + 1, ' ')

        # Fill any remaining None entries (mostly HETATM `.` slots beyond
        # the protein region) with sequential resseqs after the last
        # assigned residue.
        last_num = 0
        for i in range(len(full_resids)):
            if full_resids[i] is not None:
                last_num = full_resids[i][0]
            else:
                last_num += 1
                full_resids[i] = (last_num, ' ')

        # Deduplicate: only on the mask-based fallback path. The
        # deterministic path places residues at their natural SEQRES
        # positions with the original resseqs; running dedup would shift
        # HETATM resseqs (which can be numerically below protein resseqs)
        # to the wrong values.
        if not used_deterministic:
            for i in range(1, len(full_resids)):
                if full_resids[i][0] <= full_resids[i-1][0]:
                    full_resids[i] = (full_resids[i-1][0] + 1, ' ')

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
    # record boundaries to handle this. Also split on CONECT so that a
    # joined "TER...CONECT..." line gets disentangled and the CONECT
    # part can be dropped explicitly below.
    lines = content.split('\n')
    expanded = []
    for line in lines:
        while True:
            found = False
            for prefix in ("ATOM  ", "HETATM", "CONECT"):
                idx = line.find(prefix, 1)  # skip pos 0
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
        elif line.startswith("CONECT"):
            # Explicitly drop Modeller-emitted CONECT records. Bonds are
            # restored separately by restore_conect_records using the
            # ORIGINAL PDB's CONECT entries, remapped by atom identity.
            continue
        # Any other line (REMARK, HEADER, SEQRES, SSBOND, etc.) is dropped.

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
            if key in resnum_mapping and resnum_mapping[key] is not None:
                new_resseq, new_icode = resnum_mapping[key]
                line = _set_resid(line, new_resseq, new_icode)
            output.append(line)
        elif line.startswith("TER"):
            if len(line) > 26:
                chain = line[21]
                seq_str = line[22:26].strip()
                if seq_str and seq_str.isdigit():
                    key = (chain, int(seq_str))
                    if key in resnum_mapping and resnum_mapping[key] is not None:
                        new_resseq, new_icode = resnum_mapping[key]
                        line = _set_resid(line, new_resseq, new_icode)
            output.append(line)
        else:
            output.append(line)
    return output


def _renumber_atom_serials(lines):
    """Assign sequential unique serials (1..N) to every ATOM/HETATM line.

    Used after `patch_missing_hetatm`, which inserts atoms from the input
    PDB keeping their original serials — those collide with Modeller's
    renumbered serials in the result. Without re-serialization, the final
    PDB has duplicate serials, and CONECT records resolve to whichever
    atom the viewer happens to pick first → spurious bonds.

    Returns a new list of lines with cols 7-11 (atom serial) rewritten.
    Non-atom lines pass through unchanged.
    """
    out = []
    serial = 0
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            serial += 1
            new_serial = serial if serial < 100000 else (serial % 100000)
            out.append(line[:6] + f"{new_serial:5d}" + line[11:])
        else:
            out.append(line)
    return out


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


def restore_conect_records(result_lines, original_lines, verbose=False,
                            per_chain_masks=None, all_chains=None):
    """Restore CONECT records from original PDB with remapped atom serials.

    Modeller strips all CONECT records. This function rebuilds them by:
    1. Walking residues in chain order in BOTH the input and the model
       output to build a position-based residue index per chain.
    2. For each input residue, map its INPUT position N to its MODEL
       output position via `per_chain_masks` (True=template, False=gap).
       The Nth template-position in the mask is the Nth input residue;
       its model output position is the mask index where that True sits.
    3. Within each aligned residue, atom-name matching links input
       atoms to model output atoms regardless of serial renumbering.

    Without `per_chain_masks` (e.g. when called without alignment info),
    falls back to naive position-by-position alignment which is correct
    only when no gaps were filled.
    """
    def _build_residue_index(lines):
        """Walk all ATOM/HETATM lines and return:
          - by_chain_pos: dict {chain: [(record_type, resname, resnum, icode), ...]}
            (entries in file order, one per unique residue per chain)
          - serial_to_pos: dict serial -> (chain, position_index, atomname)
        """
        by_chain_pos = {}
        seen_per_chain = {}
        serial_to_pos = {}
        for line in lines:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            chain = line[21]
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            icode = line[26] if len(line) > 26 else ' '
            atomname = line[12:16].strip()
            try:
                serial = int(line[6:11].strip())
            except ValueError:
                continue
            res_key = (resnum, icode)
            seen = seen_per_chain.setdefault(chain, {})
            if res_key not in seen:
                idx = len(by_chain_pos.setdefault(chain, []))
                seen[res_key] = idx
                by_chain_pos[chain].append(
                    ('het' if line.startswith('HETATM') else 'atom',
                     resname, resnum, icode))
            pos = seen[res_key]
            serial_to_pos[serial] = (chain, pos, atomname)
        return by_chain_pos, serial_to_pos

    orig_by_chain, orig_serial_to_pos = _build_residue_index(original_lines)
    model_by_chain, _ = _build_residue_index(result_lines)

    # Build (chain, model_pos) -> serial map by re-walking result_lines
    model_atom_serial = {}
    seen_per_chain = {}
    for line in result_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain = line[21]
        resnum = line[22:26].strip()
        icode = line[26] if len(line) > 26 else ' '
        atomname = line[12:16].strip()
        try:
            serial = int(line[6:11].strip())
        except ValueError:
            continue
        res_key = (resnum, icode)
        seen = seen_per_chain.setdefault(chain, {})
        if res_key not in seen:
            seen[res_key] = len(seen)
        pos = seen[res_key]
        model_atom_serial[(chain, pos, atomname)] = serial

    # Build per-chain map: input_position → model_output_position.
    # When per_chain_masks is provided, use it: the Nth template-True
    # in the mask is the Nth input residue, located at the model
    # position corresponding to that mask index. When per_chain_masks
    # is missing, fall back to naive identity (works only when no gaps).
    chain_input_to_model_pos = {}  # {chain: {input_pos: model_pos}}
    if per_chain_masks is not None and all_chains is not None:
        for ci, chain in enumerate(all_chains):
            if ci >= len(per_chain_masks):
                continue
            mask = per_chain_masks[ci]
            input_pos = 0
            mapping = {}
            for model_pos, is_template in enumerate(mask):
                if is_template:
                    mapping[input_pos] = model_pos
                    input_pos += 1
            chain_input_to_model_pos[chain] = mapping

    # Sanity: warn if input vs model chain lengths don't match the mask
    # expectation. Non-protein chains aren't in the mask but should still
    # align position-by-position (Modeller doesn't reorder them).
    if verbose:
        for chain, orig_residues in orig_by_chain.items():
            model_residues = model_by_chain.get(chain, [])
            if chain in chain_input_to_model_pos:
                expected_input_n = sum(1 for v in chain_input_to_model_pos[chain].values())
                if expected_input_n != len(orig_residues):
                    print(f"  WARN: chain {chain}: mask says {expected_input_n} "
                          f"input residues but PDB has {len(orig_residues)} "
                          f"— bond restoration may miss some atoms")
            elif len(orig_residues) != len(model_residues):
                print(f"  WARN: chain {chain}: input has {len(orig_residues)} "
                      f"residues, model has {len(model_residues)} — bonds may "
                      f"be partially restored")

    # Helper: resolve an input serial to the corresponding model output serial
    def _remap(serial):
        pos_info = orig_serial_to_pos.get(serial)
        if pos_info is None:
            return None
        chain, orig_pos, atomname = pos_info
        # Translate input_pos → model_pos via the mask (or identity).
        pos_map = chain_input_to_model_pos.get(chain)
        if pos_map is not None:
            model_pos = pos_map.get(orig_pos)
        else:
            model_pos = orig_pos  # non-protein chain: position-by-position
        if model_pos is None:
            return None
        # Check resname consistency to avoid misrouting if mask is wrong
        model_residues = model_by_chain.get(chain, [])
        orig_residues = orig_by_chain.get(chain, [])
        if model_pos >= len(model_residues) or orig_pos >= len(orig_residues):
            return None
        if orig_residues[orig_pos][1] != model_residues[model_pos][1]:
            # resname mismatch — drop the bond rather than misroute it.
            return None
        return model_atom_serial.get((chain, model_pos, atomname))

    # Build a per-input-serial → per-model-serial map used below
    orig_serial_to_model_serial = {}
    for serial in orig_serial_to_pos:
        new_serial = _remap(serial)
        if new_serial is not None:
            orig_serial_to_model_serial[serial] = new_serial

    # Direct input-serial → model-serial mapping used for CONECT remap.
    orig_to_model_serial = orig_serial_to_model_serial

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
            new_s = orig_to_model_serial.get(s)
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
        try:
            fasta_map = parse_fasta(args.fasta)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)
        protein_chains = (
            get_chain_order(lines, seqres_only=set(seqres.keys()))
            if seqres else get_chain_order(lines)
        )
        missing = [c for c in protein_chains if c not in fasta_map]
        if missing:
            print(
                f"Error: FASTA missing sequences for chain(s): {', '.join(missing)}. "
                f"FASTA has: {', '.join(sorted(fasta_map.keys()))}. "
                f"PDB has: {', '.join(protein_chains)}.",
                file=sys.stderr,
            )
            sys.exit(1)
        extra = [c for c in fasta_map if c not in protein_chains]
        if extra and args.verbose:
            print(f"Note: FASTA has extra chain(s) not in PDB (ignored): "
                  f"{', '.join(extra)}")
        protein_seq_map = {ch: fasta_map[ch] for ch in protein_chains}
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

        # Write FULL PDB (all chains) for Modeller.
        # Rename non-standard protonation variants to standard names so Modeller
        # can read them (Modeller only knows HIS, not HIE/HID/HIP/ASH/GLH etc.)
        # Save the original names to restore in the output.
        _VARIANT_TO_STD = {
            'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
            'HSE': 'HIS', 'HSD': 'HIS', 'HSP': 'HIS',
            'ASH': 'ASP', 'ASPP': 'ASP',
            'GLH': 'GLU', 'GLUP': 'GLU',
            'CYX': 'CYS', 'CYM': 'CYS',
            'LYN': 'LYS',
            # GLYCAM glycoprotein residues — Modeller doesn't know NLN/
            # OLS/OLT; pre-rename to the standard parent residue so the
            # protein chain is recognized correctly. Restored on output
            # via _input_variants downstream.
            'NLN': 'ASN', 'OLS': 'SER', 'OLT': 'THR',
        }
        # Capture: (chain, resseq) -> original variant name
        _input_variants = {}
        for line in lines:
            if line.startswith(('ATOM  ', 'HETATM')):
                resname = line[17:20].strip()
                if resname in _VARIANT_TO_STD:
                    ch = line[21]
                    rs = line[22:26].strip()
                    _input_variants[(ch, rs)] = resname

        # Apply variant rename, then reorder chains so each chain ID forms
        # one contiguous block in the temp file Modeller reads. Without
        # reordering, Modeller emits one PIR segment per file-order chain
        # block — if a chain ID appears in two file positions (e.g. protein
        # ATOMs early, glycan HETATMs later), the template ends up with
        # one more segment than `all_chains` knows about, and align2d either
        # silently drops a segment (3ry6: chain C's NAG is lost) or aborts
        # with 'No aligned template residues for BLK residue' (FcgRI).
        renamed_lines = []
        for line in lines:
            if line.startswith(('ATOM  ', 'HETATM')):
                resname = line[17:20].strip()
                std = _VARIANT_TO_STD.get(resname)
                if std:
                    line = line[:17] + f"{std:>3s}" + line[20:]
            renamed_lines.append(line)
        reordered_lines = _reorder_chains_for_modeller(renamed_lines)

        input_pdb = workdir + '/' + input_path.name
        with open(input_pdb, 'w') as f:
            f.writelines(reordered_lines)

        # Run Modeller — target is built inside using Modeller's own template
        # reading, so non-protein '.' counts are guaranteed to match
        best_model, aln_path = run_modeller(
            Path(input_pdb), protein_chains, protein_seq_map, all_chains, args
        )

        # If no gaps were found, run_modeller returns the input path unchanged.
        # Just copy input to output with variant restoration and skip post-processing.
        if best_model == str(Path(input_pdb)):
            result_lines = lines[:]
            # Restore variant names
            if _input_variants:
                restored = []
                for line in result_lines:
                    if line.startswith(('ATOM  ', 'HETATM')):
                        ch = line[21]
                        rs = line[22:26].strip()
                        orig = _input_variants.get((ch, rs))
                        if orig:
                            line = line[:17] + f"{orig:>3s}" + line[20:]
                    restored.append(line)
                result_lines = restored
            if not args.keep_water:
                result_lines = remove_water_lines(result_lines)
            with open(output_path, 'w') as f:
                f.writelines(result_lines)
            print(f"Wrote {output_path}")
            # Write empty .dat (no gaps filled)
            dat = {"description": "No gaps — input copied unchanged",
                   "total_added": 0, "residue_summary": {}, "added_atoms": []}
            if _input_variants:
                dat['variant_overrides'] = {
                    f"{ch}:{rs}": var for (ch, rs), var in _input_variants.items()
                }
            dat_path = output_path.with_suffix('.dat')
            with open(dat_path, 'w') as f:
                json.dump(dat, f, indent=2)
            print(f"Wrote {dat_path}")
            return

        # Restore original chain IDs (Modeller reassigns A,B,C,... to all chains)
        result_lines = restore_chain_ids_and_read(best_model, all_chains)

        # Restore residue numbering for protein chains only
        _, per_chain_masks = get_template_mask(aln_path)
        resnum_mapping = build_resnum_mapping(
            per_chain_masks, all_chains, protein_chains, lines,
            protein_seq_map=protein_seq_map,
        )
        result_lines = renumber_model_output(result_lines, resnum_mapping)

        # Patch back any HETATM atoms Modeller dropped (e.g. NAG C1 linkage atom)
        result_lines = patch_missing_hetatm(
            result_lines, lines, all_chains, protein_chains, args.verbose
        )

        # Re-serialize ALL atoms after patching so serials stay unique.
        # patch_missing_hetatm copies orig_line[:21] + new_chain + new_resnum
        # + orig_line[27:] — which keeps the ORIGINAL input serial in cols
        # 7-11. That serial collides with Modeller's renumbered atoms,
        # producing two atoms with the same serial in the output. CONECT
        # records (added next by restore_conect_records) then resolve to
        # whichever duplicate the viewer picks first, creating spurious
        # bonds between unrelated atoms ("strange bonds inside glycans").
        result_lines = _renumber_atom_serials(result_lines)

        # Restore CONECT records from original PDB with remapped atom serials
        result_lines = restore_conect_records(
            result_lines, lines, args.verbose,
            per_chain_masks=per_chain_masks, all_chains=all_chains,
        )

        if args.verbose:
            n_gaps_filled = sum(
                not m for ci, cm in enumerate(per_chain_masks)
                for m in cm
                if ci < len(all_chains) and all_chains[ci] in set(protein_chains)
            )
            print(f"Restored numbering ({n_gaps_filled} gap positions assigned)")

        if not args.keep_water:
            result_lines = remove_water_lines(result_lines)

        # Restore original protonation variant names (HIE, ASH, etc.)
        if _input_variants:
            restored = []
            for line in result_lines:
                if line.startswith(('ATOM  ', 'HETATM')):
                    ch = line[21]
                    rs = line[22:26].strip()
                    orig = _input_variants.get((ch, rs))
                    if orig:
                        line = line[:17] + f"{orig:>3s}" + line[20:]
                restored.append(line)
            result_lines = restored

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
