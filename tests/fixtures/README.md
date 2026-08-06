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
| `hinge_CH3_glycosylated.pdb` | `test/default.pdb` | GLYCAM hinge/CH3, CYX, rename/convert/CONECT tests |
| `multistate.pdb`, `multistate.fasta` | `test/multistate/test_multistate.pdb`, `test/multistate/test.fasta` | Eleven-model split input and companion sequence |
| `8cz8/*` | `test/8cz8/*` | Pure-protein renumbering and truncated-LYS deterministic rebuild; companion FASTA |
| `lipid/*` | `test/lipid/*` | 7X35 protein/PLM residue-number collision; companion FASTA |
| `regressions/1DQJ.pdb` | `test/shit/1DQJ_original.pdb` | Disulfide-rich antibody |
| `regressions/{1EMV,1FR2,2VLN,2VLQ}.pdb` | Corresponding `test/shit/*_original.pdb` | Historical coincident-H/chirality/NaN failures |
| `trastuzumab.pdb` | `test/shit/trastuzumab.pdb` | PDB-named glycoprotein stress case |
| `1VCU.pdb` | `test/protein_ligand/1VCU.pdb` | Protein with DAN and two EPE ligands |
| `3ry6.pdb` | `test/3ry6/3ry6.pdb` | Under-annotated four-site glycoprotein |
| `glycosilated_mAb_CHARMM.pdb` | `test/glycosilated_mAb_Charmm-GUI/conf.pdb` | CHARMM-GUI glycoprotein topology input |
| `numbering/*` | `test/numbering_problem/*` | 8B01 bound/unbound structures and reference FASTA |

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
