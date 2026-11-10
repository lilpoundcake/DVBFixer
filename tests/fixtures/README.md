# Tracked test fixtures

These are the input structures and companion sequences read by the pytest
suite. They were copied from the historical, mostly untracked `test/` working
tree so CI and fresh clones exercise the same structural regressions. Generated
outputs (`*_zbs`, prepared/minimized files, GRO files, and topology trees) are
deliberately excluded.

Run `sha256sum -c tests/fixtures/MANIFEST.sha256` from the repository root to
verify that none of the source fixtures changed accidentally.

## Fixture inventory

| Tracked fixture | Historical source | Purpose |
|---|---|---|
| `ASN.pdb` | `test/ASN.pdb` | Minimal two-residue ASN input |
| `protein_context/5j13_selected.pdb` | User-supplied `test/errors_diffusion/selected.pdb`, deposited structure 5J13 after CIF normalization and alternate-state selection | Multi-gap diffusion regression: three case-sensitive chains, 10 missing regions, partial observed sidechains. This is a supplied preparation input, not a pristine RCSB download; no generated diffusion coordinates are included. |
| `protein_context/5j13_selected.fasta` | User-supplied `test/errors_diffusion/selected.fasta`, copied byte-for-byte | Authoritative target for the supplied 5J13 derivative: chain lengths A/B/C 134/213/226, seven missing regions / 42 residues. Distinct from deposited SEQRES; used by the successful seed-7 terminal-context control. No predicted sequence or coordinates are included. |
| `protein_context/9cqd_split.{pdb,fasta}` | Source PDB and companion FASTA copied byte-for-byte from user-supplied `test/errors_diffusion/diagnostics.tar.gz`, verified against its embedded checksum manifest | Diffusion metadata regression: source SEQRES declares M/W/X 49/220/261, FASTA selects 27/220/221, and `--no-terminal` targets 27/214/221. Source contains foreign-chain annotations; valid disulfides and observed coordinates must survive reconciliation. No generated coordinates are tracked. |
| `assemblies/8XJ0.pdb`, `assemblies/8XJ0_chain_A_observed.fasta` | `test/homology_modelling/8XJ0.pdb`; FASTA derived from its observed chain A | Four Fab biological assemblies declared by REMARK 350/BIOMT; post-2023 temporal-holdout diffusion pilot |
| `hinge_CH3_glycosylated.pdb` | `test/default.pdb` | GLYCAM hinge/CH3, CYX, rename/convert/CONECT tests |
| `multistate.pdb`, `multistate.fasta` | `test/multistate/test_multistate.pdb`, `test/multistate/test.fasta` | Eleven-model split input and companion sequence |
| `8cz8/*` | `test/8cz8/*` | Pure-protein renumbering and truncated-LYS deterministic rebuild; companion FASTA |
| `lipid/*` | `test/lipid/*` | 7X35 protein/PLM residue-number collision; companion FASTA |
| `regressions/1DQJ.pdb` | `test/shit/1DQJ_original.pdb` | Disulfide-rich antibody |
| `regressions/{1EMV,1FR2,2VLN,2VLQ}.pdb` | Corresponding `test/shit/*_original.pdb` | Historical coincident-H/chirality/NaN failures |
| `regressions/4tut_1uc.pdb` | Unmodified [`Amber-MD/AmberClassic` 4TUT example](https://github.com/Amber-MD/AmberClassic/blob/76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002/examples/4tut/4tut_1uc.pdb) at commit `76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002` | Same-chain protein fragments separated by `TER`; force-field atom-name conversion at both C termini |
| `trastuzumab.pdb` | `test/shit/trastuzumab.pdb` | PDB-named glycoprotein stress case |
| `1VCU.pdb` | `test/protein_ligand/1VCU.pdb` | Protein with DAN and two EPE ligands |
| `3ry6.pdb` | `test/3ry6/3ry6.pdb` | Under-annotated four-site glycoprotein |
| `glycosilated_mAb_CHARMM.pdb` | `test/glycosilated_mAb_Charmm-GUI/conf.pdb` | CHARMM-GUI glycoprotein topology input |
| `numbering/*` | `test/numbering_problem/*` | 8B01 bound/unbound structures and reference FASTA |
| `c_glh/*` | `test/C_GLH/{8cde_t_u.pdb,8cde_renamed.fasta}` | Terminal GLH/capping preparation regression input and companion FASTA |
| `warnings/*` | `test/warnings/{8ct6_t_b.pdb,8ct6_renamed.fasta}` | Real addHydrogens/connectivity-warning input and companion FASTA |
| `overlap/8dis_t_u.pdb` | `test/overlap/8dis_t_u.pdb` | Coordinate-identical chains `d`/`D` caused by missing MODEL/ENDMDL separators |
| `insertion_codes/{7K8S.pdb,7K8S.fasta}` | RCSB PDB entry 7K8S coordinate and FASTA downloads | Antibody heavy-chain insertion-code diffusion benchmark (`H/82A`-`H/82C`) |
| `protein_context/{1UBQ,1TIM,1A3N}.pdb` | Unmodified RCSB PDB downloads on 2026-10-14 from `https://files.rcsb.org/download/<ID>.pdb` | Frozen protein-context expansion: ubiquitin nonlocal-contact strata, triosephosphate-isomerase dimer, and hemoglobin protein interface |
| `protein_context/scoring_followup/{1CRN,1AKI,1BTA}.pdb` | Unmodified RCSB PDB downloads on 2026-11-03 from `https://files.rcsb.org/download/<ID>.pdb`, after predeclaring the three-system scoring follow-up | Crambin/lysozyme disulfide-preserving projection and CONECT-order regression; barstar's 714 deposited H coordinates intentionally remain unsupported by the private heavy-only completion contract. Raw deposited sources only; no generated candidates or scoring snapshots. |
| `protein_context/independent_followup/{2CI2,1R69,1CTF}.pdb` | Unmodified RCSB PDB downloads on 2026-11-03, 14:27 UTC, from `https://files.rcsb.org/download/<ID>.pdb`, after pre-download enrollment of all three sources | New-system three-chi/template scoring controls: first-model/first-ATOM-chain selection preserves 2CI2 chain I, and 1CTF has only one observed THR. Source-only regressions retain raw native gate deviations and template-stereo control inputs, not generated proposals, logs or public completion outputs. |
| `protein_context/ncac_followup/{1PGB,1BPI,1CSP}.pdb` | Unmodified RCSB downloads on 2026-11-03 from `https://files.rcsb.org/download/<ID>.pdb`, after frozen three-system NCAC/relative-rotamer enrollment | Prospective negative controls: 1BPI retains altloc/occupancy refusals and a single-VAL shortfall; 1CSP retains incomplete context and scorer movement/import blockers. Original sources only, no repairs or generated outputs. |
| `cc89_apple_handoff/{9WW8,8EK6,9Z0U}.{pdb,fasta}`, `manifest.json` | Whole observed protein references and unchanged FASTAs from the independently audited prospective cc89 cohorts; per-source RCSB CIF URLs and hashes are bound by the manifest and linked cohort manifests | Fresh-clone Apple handoff: four source-only masks, seeds 7/19; no Linux archives or checkpoint needed to materialize inputs |

The protein-context materializer derives protein-only, heavy-atom inputs from
these deposited sources. Heterogens (including 1A3N hemes) are excluded explicitly;
this is a protein-context mechanism study, not acceptance for complete hemoglobin
chemistry. Low/high contact counts are structural proxies, not measured solvent
accessibility. Derived masks, controls, and outcomes live in ignored evidence
archives; only deposited source inputs and frozen manifests are tracked.

The scoring-follow-up sources retain all original headers, bond annotations,
waters/other deposited records and (for 1BTA) explicit hydrogens. The research
projection is separate and preserves chain-A ATOM bytes and relevant linkage
metadata; it never overwrites these sources. Follow-up SHA-256 values:

- 1CRN: `42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811`
- 1AKI: `c75d7a689617248cdd92dc6633531d2506fb9bef1e6e21e26c8f579ae6955abb`
- 1BTA: `d6cf4b679fa9d9dc738d0e1ba588d58bc5d0f1758c7b186b27a5d858926673c9`

See [enrollment, explicit protocol amendment and bounded scoring results](../../docs/research/rosetta-scoring-followup.md).

The independently enrolled follow-up sources preserve every original PDB byte;
first-model/first-chain projection and deletion masks live only in private research
artifacts. Source SHA-256 values:

- 2CI2: `8a42f06a0a59074eddf9ac6274962057fc965a1256dce268d96b9d903ba4ae4d`
- 1R69: `4ceed57b7b51edfd6e873fba0aff4edf6c87d330eea318b67eabd5520d525a1b`
- 1CTF: `d0d6bbbe48bcb09c737d06d65cd0e22f86c220b3c68fa7166449e42b43294266`

See [pre-download rules, complete denominators and retained stereo failures](../../docs/research/rosetta-independent-controls.md).

The NCAC follow-up also preserves every original download byte. SHA-256 values:

- 1PGB: `1f7a38180b264722853820aa636f9ff66b458f914061f5a12c1687fd951c2f6d`
- 1BPI: `2289dd009b6eae05ac142836277e23098e178185df9441e2a3b9c1bed7dd75df`
- 1CSP: `0b3335beb32172e4e8fab1cf6533db795f88de71f46192904aff67a23a6f5b9e`

The private enrollment and all generated evidence remain below
`tmp/pdbfixer-research/ncac-new-system-*`. The tracked sources are not selected
successes: unsupported observations and incomplete context are intentionally kept.

The Apple handoff sources are the already-normalized experimental input frames:
protein heavy atoms only, water/H excluded, original observed chain sequences and
coordinates retained, residue IDs densely mapped by the audited mosaic builder.
They are not generated model outputs. `scripts/prepare_apple_cc89_handoff.py`
recreates masked inputs and separates evaluation references. See
[`docs/plans/apple-cc89-handoff.md`](../../docs/plans/apple-cc89-handoff.md).

## 4TUT TER-boundary naming

`regressions/4tut_1uc.pdb` is the complete upstream AmberClassic example,
copied byte-for-byte from the commit-pinned source above (SHA-256
`ad65166d9f69cff807ddc915730cd7670ecfc8dcb98072b0ef58eec7f43cc030`).
It contains two protein fragments that reuse chain `A` but have distinct residue
numbers and are separated by `TER` records. The regression verifies only that
`atom-names --target-ff amber --profile gromacs` recognizes both fragment
termini while preserving every non-name field. The source file has no
hydrogens, so it does not exercise N-terminal hydrogen renaming. It also makes
no claim about force-field suitability, topology generation, or molecular
dynamics behavior.

## 8UCD molecule retention and connectivity

`rename_mols/structure.pdb` is the user-provided 8UCD preparation input copied
unchanged from `test/rename_mols/structure.pdb`. Split regression tests verify
that all 807 heterogen atoms (nine LBN/FAD/HEM residues) survive retention and
unique-chain assignment with unchanged coordinates and serials.

## 1VCU ligand chemistry

The optional-SMILES regression uses RCSB Chemical Component Dictionary
connectivity and DAN stereochemistry. The selected forms represent the
dominant DAN carboxylate and a representative HEPES/EPE zwitterion near
physiological pH, so the expected acidic oxygens remain unprotonated:

- DAN: `CC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)[O-])O`
- EPE: `C1CN(CC[NH+]1CCO)CCS(=O)(=O)[O-]`

Sources: [RCSB DAN](https://www.rcsb.org/ligand/DAN),
[RCSB EPE](https://www.rcsb.org/ligand/EPE), and
[PubChem HEPES (CID 23831)](https://pubchem.ncbi.nlm.nih.gov/compound/Hepes).
At pH 7.4, HEPES is close to its approximately 7.55 buffer pKa, so protonated
zwitterionic and deprotonated anionic microspecies coexist; the test chooses
the slightly more populated zwitterion. A protein binding site can shift this
balance, so users should supply the microspecies appropriate to their system.
