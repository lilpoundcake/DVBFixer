"""Modeller invocation for ``dvbfixer model``.

Split out of the flat ``model.py`` in the Phase 2.2 follow-up work.
Owns everything on the Modeller side of the pipeline: PIR generation
and parsing, the LoopModel subclass that pins the input during MD,
the terminal-alignment fixer, the chain-block reorderer that keeps
Modeller from emitting duplicate template segments for disjoint chain
IDs, and the ``run_modeller`` orchestrator itself. Also the
:func:`_explain_modeller_error` translator that surfaces useful causes
for the three common Modeller failure classes.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
from pathlib import Path


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

def _add_chirality_restraints(model, verbose: bool = False) -> tuple[int, int]:
    """Add an improper Gaussian restraint on CA-N-C-CB to every non-GLY
    residue of a Modeller ``LoopModel``/``AutoModel`` instance.

    Target mean +34° (L configuration), stdev 5°. The dihedral atom
    order CA-N-C-CB matches the AMBER canonical Cβ chirality improper —
    L geometry gives θ ≈ +34°, D gives θ ≈ -34°. A previous revision
    used (N, CA, C, CB) which measures a DIFFERENT dihedral (L ≈ -126°,
    D ≈ +145°) — with mean=+34° the restraint was actually pulling
    residues TOWARDS D. Fixed 2026-07-24 after the 1FR2 zbs output
    surfaced 40+ D residues per chain.

    Returns ``(n_added, n_skipped)``. Missing atoms (e.g. BLK residues
    or targets without CB in the alignment) are skipped silently.
    """
    import math

    from modeller import features, forms, physical

    n_added = 0
    n_skipped = 0
    for res in model.residues:
        try:
            resname = getattr(res, 'pdb_name', None) or getattr(res, 'name', '')
        except Exception:
            resname = ''
        if resname == 'GLY':
            n_skipped += 1
            continue
        try:
            n_atom = res.atoms['N']
            ca = res.atoms['CA']
            c = res.atoms['C']
            cb = res.atoms['CB']
        except (KeyError, IndexError, AttributeError):
            n_skipped += 1
            continue
        try:
            # Atom order CA-N-C-CB — AMBER canonical Cβ chirality
            # improper convention (L ≈ +34°, D ≈ -34°).
            # stdev=2° is tight: energy penalty at ±34° from mean
            # (i.e. at planar 0°) is (34/2)²/2 = 144 pseudo-kcal/mol,
            # enough to overwhelm any local packing preference for D.
            # stdev=5° was insufficient — GLN85 on 2VLQ still flipped
            # to D in ~1/3 of Modeller runs during 2026-07-24 testing.
            model.restraints.add(forms.Gaussian(
                group=physical.improper,
                feature=features.Dihedral(ca, n_atom, c, cb),
                mean=math.radians(34.0),
                stdev=math.radians(2.0),
            ))
            n_added += 1
        except Exception:
            n_skipped += 1
    if verbose:
        print(f"[modeller] chirality restraints: {n_added} added, "
              f"{n_skipped} skipped")
    return n_added, n_skipped


def run_modeller(input_path, protein_chains, protein_seq_map, all_chains, args):
    """Run Modeller loop modeling. Returns (best_model_path, alignment_path).

    Uses env.io.hetatm=True so non-protein atoms (glycans, ligands) are
    read into the model and preserved through loop modeling.

    The target sequence for non-protein chains is derived from Modeller's own
    template reading to guarantee the '.' counts match exactly.
    """
    from modeller import Alignment, Environ, Model, Selection, log
    from modeller import automodel as am
    from modeller.automodel import LoopModel

    _verbose_restraints = bool(getattr(args, 'verbose', False))

    class _LoopModelWithChirality(LoopModel):
        """LoopModel with a Gaussian improper restraint on N-CA-C-CB for
        every non-GLY residue (target +34°, stdev 5°). Prevents loop MD
        refinement from settling into the D basin — the root cause of the
        occasional D-Cα that survived to the ``model`` output before.
        """
        def special_restraints(self, aln):
            super().special_restraints(aln)
            _add_chirality_restraints(self, verbose=_verbose_restraints)

    class _PinnedLoopModel(_LoopModelWithChirality):
        """LoopModel variant that lets ONLY gap residues move during
        loop refinement MD.

        Overrides only `select_loop_atoms` — the loop refinement stage
        that produces the visible ±~3 residue flank drift. Deliberately
        does NOT override `select_atoms` (initial automodel CG). Doing so
        (earlier version of this class did) held every non-gap atom at
        its raw input coords, which preserved any pre-existing input
        close contacts. Downstream `dvbfixer prepare` runs OpenBabel
        `ConnectTheDots` on the model output, which then infers a
        spurious bond across a close contact and trips OpenMM template
        matching (e.g. `1 N atom too many externally bonded` on a GLN
        far away from any gap). Letting the initial automodel CG run on
        all atoms lets those contacts relax normally.

        Uses `self.loops(..., insertion_ext=0, deletion_ext=0)` to strip
        the default extension margin that stock LoopModel adds.
        """
        def select_loop_atoms(self):
            aln = self.read_alignment()
            loops = self.loops(aln, minlength=1, maxlength=9999,
                               insertion_ext=0, deletion_ext=0,
                               include_termini=True)
            sel = Selection(loops).only_std_residues()
            return sel if len(sel) > 0 else Selection(self)

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
        # Return the same shape as the normal path: list[(path, molpdf)].
        # Fake molpdf = 0.0 since we never actually scored anything. The
        # caller detects the shortcut by comparing this first path to
        # str(Path(input_pdb)) — see main()'s `no_gap_shortcut` flag.
        return [(str(input_path), 0.0)], aln_path

    print(f"Alignment: {n_gaps} gap residue(s) to rebuild")

    md_levels = {
        "none": None,
        "fast": am.refine.fast,
        "slow": am.refine.slow,
        "very_slow": am.refine.very_slow,
        "slow_large": am.refine.slow_large,
    }

    LoopModelCls = (_PinnedLoopModel if getattr(args, 'pin_input', True)
                    else _LoopModelWithChirality)
    a = LoopModelCls(env,
                     alnfile='alignment.pir',
                     knowns=pdb_name,
                     sequence='target')
    if args.verbose and LoopModelCls is _PinnedLoopModel:
        print("Pinning input residues: only gap residues will move during MD "
              "refinement (--no-pin-input to disable)")
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
        candidates_sorted = sorted(init_models, key=lambda x: x['molpdf'])
        print("Warning: loop refinement failed, using initial model(s)")
    else:
        candidates_sorted = sorted(loop_models, key=lambda x: x['molpdf'])

    # Post-make() Cα chirality sweep: even with the special_restraints
    # improper term, a candidate can carry a residual D geometry if the
    # residue was outside the loop-refinement selection (initial CG
    # can leave marginal D at low weight). Reflect and rewrite each
    # candidate PDB. Log any repair — with the restraint active this
    # should be zero.
    import numpy as _np
    from openmm.app import PDBFile
    from openmm.unit import Quantity, nanometer

    from dvbfixer.ffutils.geometry import fix_ca_chirality
    for c in candidates_sorted:
        cand_path = c['name']
        try:
            _cand = PDBFile(cand_path)
        except Exception as e:
            # This candidate skips the chirality sweep entirely — it can
            # still reach the output with a residual D-Cα this pass would
            # have caught. The pipeline's own later `assert_all_l` call is
            # a second safety net, but a silent skip here leaves no
            # breadcrumb pointing back at THIS candidate as the cause.
            print(f"  WARNING: [modeller] could not load {cand_path} for "
                  f"the post-refinement chirality sweep ({type(e).__name__}: "
                  f"{e}) — skipping it; any residual D-Cα here will only "
                  f"surface later, if at all.")
            continue
        # Numpy-backed Quantity: fix_ca_chirality mutates via item
        # assignment (supported), and PDBFile.writeFile's np.isnan
        # guard passes regardless of the OpenMM version (older
        # builds don't strip units before np.isnan).
        _pos_arr = _np.array(
            [list(_cand.positions[i].value_in_unit(nanometer))
             for i in range(len(_cand.positions))],
            dtype=float,
        )
        _pos = Quantity(_pos_arr, nanometer)
        _n = fix_ca_chirality(_cand.topology, _pos, verbose=args.verbose)
        if _n:
            with open(cand_path, 'w') as _f:
                PDBFile.writeFile(_cand.topology, _pos, _f, keepIds=True)
            print(f"  [modeller] repaired {_n} D-Cα in {cand_path}")

    # Return a list of (path, molpdf) tuples sorted by molpdf ascending
    # (best first). Caller picks the top-N via --num-output.
    candidates = [(c['name'], c['molpdf']) for c in candidates_sorted]
    best_name, best_molpdf = candidates[0]
    print(f"Best model: {best_name} (molpdf={best_molpdf:.1f})"
          + (f" of {len(candidates)} candidates" if len(candidates) > 1 else ""))
    return candidates, aln_path
