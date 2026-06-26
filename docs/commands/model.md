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
