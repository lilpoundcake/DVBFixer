# dvbfixer homology — Multi-Template Homology Modeling

[← command index](index.md) · [← README](../../README.md)

Multi-template homology modeling with Modeller. Takes a target FASTA (multi-chain) and one or more template PDB files. Auto-aligns target to templates via pairwise `align2d` per chain (or `--salign` for structure-based). Each target chain is modeled independently against its best template chain, then assembled into a multi-chain PDB. Point mutations are handled naturally by the differing target sequence. Antibody mode (`--antibody`): uses ANARCI for Kabat/IMGT numbering, CDR detection, and auto-mapping of Fv/constant domains to different templates.

## Usage

```bash
# Basic multi-template
dvbfixer homology target.fasta --template fab.pdb --template fullsize.pdb -v

# With pipeline (prepare + minimize)
dvbfixer homology target.fasta --template fab.pdb --template fullsize.pdb --minimize -v

# Antibody mode
dvbfixer homology target.fasta --template fab.pdb --template igg.pdb --antibody -v
```

## See also

- [`model`](model.md) — single-template loop rebuilding
- [`prepare`](prepare.md) — post-modeling structure preparation
- [BEST_PRACTICES.md](../../BEST_PRACTICES.md) — antibody workflow recipe

## How it works
Multi-template homology modeling with Modeller. Takes a target FASTA (multi-chain) and one or more template PDB files. Auto-aligns target to templates via Modeller's `align2d` (or `salign` with `--salign`). Builds model with `automodel` or `LoopModel` using multiple `knowns`. Point mutations handled naturally by differing target sequence from templates. Post-processing restores chain IDs and residue numbering. Writes `.dat` file for downstream `prepare`/`minimize` restraints. `--prepare` and `--minimize` flags run the full pipeline automatically. Antibody mode (`--antibody`): uses ANARCI for Kabat/IMGT numbering, CDR detection, VH/VL/CH/CL domain classification, and auto-mapping of Fv from one template + constant domains from another. Dependencies: Modeller (required), ANARCI (for `--antibody` mode).
