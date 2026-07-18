# dvbfixer model — Loop/Gap Rebuilding with Modeller

[← command index](index.md) · [← README](../../README.md)

Rebuilds missing loops and gaps using Modeller's LoopModel. Identifies missing regions by aligning ATOM records to the SEQRES sequence (or a user-provided FASTA), then runs Modeller's loop modeling with MD refinement to fill them.

Non-protein chains (glycans, ligands) are included in the Modeller pipeline via `env.io.hetatm=True` with `'.'` (BLK residue) entries in the target sequence, so Modeller preserves them through loop modeling. Original chain IDs and residue numbering are restored automatically.

Writes a `.dat` file recording all atoms in rebuilt gap residues. This is merged by `dvbfixer prepare` with its own additions, so that `dvbfixer minimize` applies appropriate restraints to all rebuilt regions. Water molecules are removed by default (`--keep-water` to preserve).

**Key behaviour:**
- **Robust multi-chain handling** — chains are reordered before Modeller so disjoint file-order blocks sharing a chain ID (ATOM protein + HETATM glycans split by other chains) are grouped into one contiguous segment. Fixes silent glycan drops (e.g. 3ry6 chain C, FcgRI chain A).
- **Deterministic residue numbering** — gap-fill residues are numbered from input resseq jumps (not from `align2d`'s score-based placement). N-terminal extras extend backward, internal gaps fill between input neighbours, C-terminal extras extend forward — all via `first_resseq + N - K`. Mutation-tolerant via Needleman-Wunsch fallback. HETATMs attached to a protein chain keep their original resseqs.
- **`--fasta` headers must encode chain IDs**: `>chain_X`, `>PDBID_X` (e.g. `>1abc_A`), or `>X`. Matched by ID, not file order. Clear error if unparseable.
- **Plain-language Modeller diagnostics** — common errors (BLK alignment, sequence difference, unknown residue type) get a clear cause + remediation alongside the raw Modeller message.
- **Input residues pinned by default** — a `_PinnedLoopModel` subclass overrides `select_loop_atoms()` so only the newly-modeled gap residues move during Modeller's loop refinement MD. Prevents the ±~3 residue flank drift produced by stock LoopModel. The initial automodel CG is left untouched (all atoms) so any pre-existing input close contacts get relaxed as usual — otherwise they'd survive into the model output and trip OpenBabel's CONECT inference in the downstream `prepare` step (spurious external bonds → OpenMM template match failures on residues far from any gap). Downstream `minimize` still refines the whole system under a proper force field. Pass `--no-pin-input` to restore the legacy behaviour.
- **Fast no-gap pre-check** (since v0.3) — before invoking Modeller, model.py compares each protein chain's ATOM sequence to its SEQRES via string equality. If every chain matches, no gaps → copy input verbatim and skip Modeller entirely. Runtime drops from ~2 min (align2d on a 5-chain × 829-residue complex is O(N²)) to ~0.25 s. Falls through to the full Modeller path when any chain has a mismatch (real gaps).

## Usage

```bash
# Basic usage — writes input_model.pdb
dvbfixer model input.pdb -v

# Higher quality (more sampling, slower)
dvbfixer model input.pdb --num-models 2 --num-loops 4 --md-level slow -v

# Save the top 5 candidates instead of just the best one
# (useful for ensemble analysis or visual inspection)
dvbfixer model input.pdb --num-models 2 --num-loops 5 --num-output 5 -v
# → input_model_1.pdb, input_model_2.pdb, ..., input_model_5.pdb
#   (sorted ascending by Modeller's molpdf — best first; matching .dat files)

# Use FASTA instead of SEQRES for complete sequence
# (FASTA headers must encode chain IDs: >chain_A, >1abc_A, or >A)
dvbfixer model input.pdb --fasta sequence.fasta -v

# Keep Modeller working directory for debugging
dvbfixer model input.pdb --keep-workdir -v
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o`, `--output` | `<input>_model.pdb` | Output file path |
| `--fasta` | none | FASTA file with complete sequence(s) (alternative to SEQRES) |
| `-n`, `--num-models` | 1 | Number of initial models to generate |
| `--num-loops` | 2 | Number of loop refinement models per initial model |
| `--num-output` | 1 | Number of top-ranked candidates to save (ceiling: `num_models × num_loops`). Sorted ascending by Modeller's `molpdf` (best first). With `--num-output > 1`, output filenames get a `_N` suffix |
| `--md-level` | fast | MD refinement level: none, fast, slow, very_slow, slow_large |
| `--pin-input` / `--no-pin-input` | **on** | Freeze every input-structure residue during MD refinement — only gap residues move. Prevents flanking-residue drift; the downstream `minimize` step refines the whole system properly. Pass `--no-pin-input` for the legacy LoopModel behaviour (gap ±~3 residue flank mobile). |
| `--no-terminal` | off | Do not model missing N/C terminal residues (only rebuild internal gaps) |
| `--keep-water` | off | Keep water molecules (HOH, WAT, TIP3, SOL) — removed by default |
| `--keep-workdir` | off | Keep Modeller temp directory |
| `-v`, `--verbose` | off | Print Modeller progress |

## How It Works

1. Parses SEQRES (or FASTA) for the complete target sequence
2. Builds a target PIR with protein chains from SEQRES and non-protein chains as `'.'` (BLK residues)
3. Reorders chains so disjoint file-order blocks sharing a chain ID (e.g. protein + HETATM glycans on the same chain split by other chains) are grouped into one contiguous segment before Modeller runs
4. Reads the full PDB with `env.io.hetatm=True` so non-protein atoms are included
5. Creates a PIR alignment between the structure (with gaps as `-`) and the target sequence using Modeller's `align2d`
6. Runs `LoopModel` to generate initial model(s) filling the gaps
7. Refines loop regions with configurable MD level
8. Selects the best model by lowest `molpdf` score
9. Restores original chain IDs (Modeller A,B,C,... -> original letters)
10. Restores original residue numbering: template positions keep their original `(resSeq, iCode)`; gap-filled positions are numbered deterministically from input resseq jumps (N-terminal extras extend backward, internal gaps fill between input neighbours, C-terminal extras extend forward — all via `first_resseq + N - K`). HETATM residues attached to a protein chain keep their original resseqs.

## See also

- [`renumber`](renumber.md) — run before `model` when input has insertion codes
- [`prepare`](prepare.md) — next step after `model`, merges the `.dat` file
- [`homology`](homology.md) — multi-template alternative
- [ARCHITECTURE.md](../../ARCHITECTURE.md) — `build_resnum_mapping`, `_reorder_chains_for_modeller`, `parse_fasta` internals

## How it works
Rebuilds missing loops/gaps using Modeller's LoopModel. Takes SEQRES (or --fasta) as complete sequence, aligns to ATOM records via `align2d`, runs loop modeling with configurable MD refinement. Non-protein chains (glycans, ligands) are included in the Modeller pipeline via `env.io.hetatm=True` with `'.'` (BLK) entries in the target PIR sequence so Modeller keeps them through loop modeling. Post-processing restores: (1) original chain IDs from Modeller's A,B,C,..., (2) original residue numbering with insertion codes using the alignment. `--no-terminal` trims N/C terminal missing residues from the target sequence. Terminal alignment is auto-fixed after `align2d` to prevent misplaced terminal gaps. Writes a `.dat` file recording all atoms in rebuilt gap residues — prepare merges this with its own additions. Water removed by default (`--keep-water` to preserve). Protonation variant names (HIE/HID/HIP/ASH/GLH etc.) are renamed to standard (HIS/ASP/GLU) before Modeller reads the PDB (Modeller only knows standard names), then restored in the output.

**Fast Python-only no-gap pre-check** (runs before Modeller loads the PDB): compares each protein chain's ATOM-derived sequence to its SEQRES via string equality. If every chain matches, no gaps → skip Modeller entirely, copy input verbatim + write empty `.dat`. Runtime drops from ~2 min (align2d on a 5-chain × 829-residue TCR/MHC complex is O(N²) DP) to ~0.25 s. The old `run_modeller` no-gap shortcut still exists as defensive fallback and now returns `[(str(input_path), 0.0)]` (list-of-one-tuple, same shape as the normal-path return) so the caller's `first_path, _ = candidates[0]` unpacking works uniformly. Prior bug: the shortcut returned raw `str(input_path)` → caller tried to unpack the first char of the path → `ValueError: not enough values to unpack (expected 2, got 1)`.

**`--pin-input` / `--no-pin-input`** (default ON): during Modeller's LOOP refinement MD, allow only the newly-modeled gap residues to move — no ±flank margin. Implemented via `_PinnedLoopModel(LoopModel)` subclass inside `run_modeller()` that overrides only `select_loop_atoms()` (loop refinement MD). Body: `self.loops(aln, minlength=1, maxlength=9999, insertion_ext=0, deletion_ext=0, include_termini=True)` + `Selection(loops).only_std_residues()` — stripping the ±~3 residue flank margin that stock `LoopModel.select_loop_atoms` adds via `insertion_ext=2, deletion_ext=1`. Conditional instantiation: `LoopModelCls = _PinnedLoopModel if getattr(args, 'pin_input', True) else LoopModel`. **Deliberately does NOT override `select_atoms`** (the initial automodel CG). An earlier version of this class did — the effect was that every non-gap atom (both flankers and entire other chains without gaps) was held at its raw input coordinates through the whole pipeline. Any pre-existing close contact in the input (crystal packing artefact, upstream tool output, etc.) then survived into the model output. `dvbfixer prepare` runs OpenBabel `ConnectTheDots` on that output, which infers a spurious bond across the close contact — most commonly a bogus external N bond on a random GLN — and OpenMM's addHydrogens then fails with `1 N atom too many externally bonded`, even for residues far from any gap. Letting the initial automodel CG run on all atoms lets those contacts relax normally, matching stock behaviour outside the loop refinement stage. The loop refinement MD is the dominant source of user-visible flank drift, so restricting only that stage is sufficient for the pinning goal. Downstream `minimize` still refines the entire system under a proper AMBER/CHARMM force field, so `--no-pin-input` remains available for the legacy full-flank-mobile behaviour. Threaded through `zbs` (same `BooleanOptionalAction` default True; forwards `--no-pin-input` when the user opts out). No-op on the no-gap shortcut path (Modeller isn't called there).

**`--num-output N`** (default 1): save the top-N candidate models instead of just the best. `run_modeller()` now returns a list of `(path, molpdf)` tuples sorted ascending (best first). `main()` wraps the post-processing block in a `for idx, (candidate_path, molpdf) in enumerate(candidates[:n_save])` loop. With `--num-output 1` the output filename is unchanged (`<stem>_model.pdb` + `<stem>_model.dat`); with `--num-output > 1` each output gets a `_N` suffix (`<stem>_model_1.pdb`, `_model_2.pdb`, ...) and matching `.dat`. The same chain-ID / residue-renumbering / variant restoration / CONECT / water-filter helpers are applied to each candidate. Validation: `--num-output < 1` errors fast (before Modeller); `--num-output > num_models × num_loops` warns but doesn't fail. The no-gap shortcut (input copied verbatim, no Modeller) emits an INFO line when `--num-output > 1` is set and falls back to a single output (no alternatives exist without gaps).

**Chain block reordering** (`_reorder_chains_for_modeller`, runs before writing the temp PDB Modeller reads): when a chain ID appears in two disjoint file-order blocks (e.g. chain A protein ATOMs early, chain A HETATM glycan late, separated by chains B/C/D/E), Modeller emits an extra template segment for the duplicate appearance. The reorderer groups each chain's records into one contiguous block (preserves chain first-appearance order, drops original TER/CONECT, emits one TER after each block). Fixes silent NAG drops on 3ry6 chain C (the glycosylation was being lost) and 'No aligned template residues for BLK residue' crashes on FcgRI.

**Terminal-alignment dot pairing** (in `_fix_terminal_alignment`): the terminal-fix subseq matcher previously treated target's trailing `.` (HETATM slot for an attached glycan) as a suffix gap, dropping the template's matching `.` and triggering BLK errors during `a.make()`. It now strips leading/trailing dots from both target and template before protein-letter matching, then re-attaches them 1:1 (pads any imbalance with `-`).

**Error diagnostics** (`_explain_modeller_error`, wired into both the `a.make()` exception handler and the all-models-failed exit path): recognises three Modeller error classes and prints a plain-language cause + remediation alongside the raw Modeller message — 'No aligned template residues for BLK residue' (chain ID in disjoint file blocks; dvbfixer reorders automatically but still warns), 'Sequence difference between alignment and pdb' (target/template misalignment, usually the same root cause), and 'Residue type ... too long: BLK' (unknown HETATM resname).

**FASTA chain-ID matching** (`parse_fasta` returns dict keyed by chain ID; `main()` maps PDB chains to FASTA sequences by chain ID, not file order): headers are parsed in priority order `>chain_X` / `>chainX`, then `>PDBID_X` (e.g. `1abc_A`), then `>X`. `parse_fasta` raises `ValueError` on missing, duplicate, or unparseable headers; `main()` reports any PDB chains absent from the FASTA dict (with both lists) so the user can correct headers. The `--fasta` help text documents the accepted header forms.

**Residue numbering — three-tier deterministic strategy** (`build_resnum_mapping`, `_find_seqres_offset_by_resseq`, `_align_atoms_to_seqres`, with the gap-fill logic extracted into reusable `_interpolate_gaps`): replaces the old brittle align2d-mask-based renumbering, which depended on BLOSUM gap placement and could shift modeled residues to the C-terminus when input resseq jumps disagreed with align2d's gap positions (FcgRI chain A regression: 5 missing residues at 219-223 ended up numbered 283-287). The new strategy is tried in order — (1) **K-finder (preferred)**: picks K = number of N-terminal SEQRES residues absent from ATOM, then assigns SEQRES position N → `resseq = first_resseq + N - K`. Trusts input resseq jumps as authoritative gap positions, tolerates up to 10% letter mismatches (mutations), and handles N-term extras, internal gaps, and C-term extras uniformly via the same formula. Only requirement: no insertion codes in input. (2) **NW (mutation-tolerant fallback)**: semi-global Needleman-Wunsch with affine gaps (open=-10, extend=-1, X-neutral, free end gaps on the SEQRES side). Used when K-finder fails due to excessive letter mismatches. (3) **align2d mask (last resort)**: walks Modeller's alignment mask, consumes input residues in order, and interpolates gaps between flanking originals via `_interpolate_gaps`. Kept for chains that defeat both deterministic paths (icode-bearing antibody Kabat numbering, etc.). The dedup loop that resolves duplicate (resSeq, iCode) tuples only runs on the fallback (mask-based) path — running it on the deterministic paths would shift HETATM resseqs (which sit numerically below protein) to wrong values.

**N-terminal gap numbering** (in `build_resnum_mapping`): N-terminal gaps now extend backward from the first template residue. Example — SEQRES has 2 extra N-term residues, input starts at resseq 235 → new residues are numbered 233, 234, then 235... Previously the gap-fill default `left=0` made N-term gaps numbered 1, 2, ..., producing a huge numbering jump at the N-term junction.

**HETATM resseq preservation** (in the deterministic paths of `build_resnum_mapping`): HETATMs attached to a protein chain (N-linked NAG and friends) keep their ORIGINAL input resseq instead of being renumbered to the next sequential integer after the last protein residue. Verified on FcgRI chain A (NAG stays at resseq 4) and 3ry6 chain C (NAG stays at resseq 206).
