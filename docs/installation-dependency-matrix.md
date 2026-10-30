# Installation dependency matrix

Generated from `src/dvbfixer/data/installation-profiles.toml`.

Desktop inherits every Standard scientific dependency. `science` means both profiles;
`gui` means Desktop only. Other groups are explicit optional add-ons.

| Package / constraint | Source | Group | Consumers | Runtime / manual check |
|---|---|---|---|---|
| `python >=3.11,<3.14` | conda-forge | science | CLI interpreter; PROPKA annotation compatibility | python-version |
| `numpy <2.5` | conda-forge | science | geometry, alignment, RESP fit; Python 3.11 typing baseline | numpy |
| `python <3.13` | conda-forge | ambertools-compat | NumPy 1.26 wheels require Python <=3.12; intersect Standard's existing Python bounds |  |
| `numpy >=1.26,<2` | conda-forge | ambertools-compat | AmberTools proprep, ndfes, fetkutils and edgembar declared Python dependencies |  |
| `biopython >=1.83,<1.86` | conda-forge | ambertools-compat | AmberTools proprep declared Python dependency |  |
| `pdb2pqr` | PyPI | ambertools-compat | AmberTools proprep and packmol-memgen declared Python dependency | pdb2pqr |
| `libblas=*=*openblas` | conda-forge | science | NumPy/SciPy native linear algebra | linalg |
| `libcblas=*=*openblas` | conda-forge | science | NumPy/SciPy native linear algebra | linalg |
| `liblapack=*=*openblas` | conda-forge | science | NumPy/SciPy native linear algebra | linalg |
| `scipy` | conda-forge | science | diagnose sterics, RESP and alignment | scipy |
| `threadpoolctl` | conda-forge | science | in-process native BLAS runtime verification; MDAnalysis | threadpoolctl |
| `openmm >=8.1` | conda-forge | science | prepare, protonate, minimize, topology | openmm |
| `pdbfixer >=1.9` | conda-forge | science | legacy heavy-atom repair | pdbfixer |
| `openmmforcefields ==0.15.1` | PyPI | science | lig_params.GAFFTemplateGenerator | openmmforcefields |
| `tinydb` | conda-forge | science | openmmforcefields GAFF template cache (undeclared by PyPI metadata) | tinydb |
| `lxml` | conda-forge | science | openmmforcefields force-field XML processing | lxml |
| `pyyaml` | conda-forge | science | openmmforcefields force-field configuration | yaml |
| `validators` | conda-forge | science | openmmforcefields force-field source validation | validators |
| `openff-toolkit-base` | conda-forge | science | lig_params OpenFF Molecule and AmberTools wrapper | openff.toolkit |
| `rdkit` | conda-forge | science | prepare glycan/SMILES ligand chemistry | rdkit |
| `biopython` | conda-forge | science | sequence_alignment, align, salign, homology | Bio |
| `gemmi >=0.7.5` | conda-forge | science | structure_input loss-aware PDBx/mmCIF normalization | gemmi |
| `mdanalysis` | conda-forge | science | cluster, align, GRO and trajectory inputs | MDAnalysis |
| `modeller` | salilab | science | model, homology, explicit SALIGN engine | modeller; Academic license: https://salilab.org/modeller/registration.html |
| `parmed` | conda-forge | science | acpype_export and AMBER/GROMACS topology conversion | parmed |
| `plotly` | conda-forge | science | cluster --plot HTML reports | plotly |
| `hmmer` | bioconda | science | ANARCI runtime dependency (hmmscan) | hmmscan |
| `anarci` | bioconda | science | antibody and homology numbering/CDR selection | anarci; anarci-data |
| `mafft` | conda-forge | science | msa auto/default, salign and Homology GUI | mafft |
| `openbabel-wheel >=3.1.1.23` | PyPI | science | small-molecule CIF, SDF extraction, heterogen H, obminimize | openbabel.openbabel; obabel; obminimize |
| `ambertools >=23` | conda-forge | science | GAFF, RESP and tleap-reduce; independent of Reduce | antechamber; parmchk2; sqm; tleap |
| `reduce` | bioconda | science | MolProbity tautomer/proton assignment and tleap-reduce | reduce |
| `acpype` | PyPI | science | top --acpype and AMBER/GROMACS export | acpype; acpype |
| `pip` | conda-forge | science | managed-environment Python package installation | pip |
| `propka >=3.5` | PyPI | science | prepare and protonate pKa decisions | propka |
| `pyscf` | PyPI | science | parametrize --charge-method resp --qm-engine pyscf | pyscf; resp |
| `nodejs >=22` | conda-forge | gui | GUI lockfile build and standalone host | node; npm |
| `muscle >=5,<6` | bioconda | msa-alternatives | msa explicit MUSCLE 5 engine | muscle |
| `clustalo` | bioconda | msa-alternatives | msa explicit Clustal Omega engine | clustalo |
| `probe =2.18` | bioconda | probe | optional diagnose engine; Python fallback when unavailable | probe; probe-version |

Channels, in priority order: conda-forge, bioconda, salilab.

## Managed add-on platforms

These restrictions apply to provisioning, not to doctor checks for an operator-provided tool.

| Add-on | Managed platforms | Manual alternative |
|---|---|---|
| `msa-alternatives` | linux-64, osx-arm64 | — |
| `ambertools-compat` | linux-64, osx-arm64 | [Installation guide](installation.md#ambertools-auxiliary-metadata-compatibility) |
| `probe` | linux-64 | [Installation guide](installation.md#native-apple-probe-alternative) |

Version/build and execution evidence: [profile acceptance](installation-profile-acceptance.md).
Licenses/notices remain governed by upstream distributions and `THIRD_PARTY_NOTICES/`.
