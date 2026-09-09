# Whole-complex relaxation: evidence and proposed work

Status: research recorded for 0.8.5, 2026-09-09. No new minimization backend is
implemented in this release. Repository observations below were checked against
the release source and tracked fixture. External engine capabilities are drawn
from the linked documentation; comparative performance and proposed designs
have not been benchmarked here.

## Problem and verified current behavior

A useful preparation result must preserve component identity, connectivity,
observed coordinates, and chemically meaningful geometry across protein,
ligands, glycans, cofactors, solvent, and ions. Keeping every atom in a PDB and
successfully minimizing every retained component are separate requirements.

The Split regression input [8UCD](../../tests/fixtures/rename_mols/structure.pdb)
contains 807 HETATM records in nine LBN/FAD/HEM residues. Its provenance and
checksum are recorded in the [fixture catalog](../../tests/fixtures/README.md).
The empirical Split writer previously changed atom serials while copying CONECT
unchanged. Version 0.8.5 preserves those serials; regression tests compare bond
endpoint identities through Split and Model's restoration routine. Unique-chain
naming uses PDB residue boundaries, not general molecular graph partitioning.
Viewer loading now selects the complete first coordinate model so deposited
assembly chain lists do not hide newly named ligand chains.

`minimize` attempts OpenMM whole-system construction. Native templates and
validated additional XML can provide parameter coverage. Unknown eligible
organic ligands trigger automatic GAFF2/AM1-BCC parameterization; the explicit
flag makes ordinary failures strict. Unknown HEM/FAD/FMN/NAD/NAP cofactors are
rejected by the domain guard even in automatic mode. Template-name screening
precedes OpenMM's atom/bond matching, so a same-name template can still fail.
These observations come from `lig_params.py`, `domain/parameterization.py`, and
`minimize/pipeline.py`, not a claim of complete chemical validation.

Some other setup failures still fall back to protein-only strip-and-splice.
Noncovalent ligands can then remain at their original positions while the pocket
moves; covalent glycan trees have a separate rigid-tracking path. Heterogen heavy
atoms use weak positional restraints in whole-system optimization. Hydrogens
are free. Chirality guards and post-repair checks remain mandatory. A zero-D-Cα
result alone does not demonstrate good ligand geometry or an unclashed interface.

## Incomplete FAD: input evidence

The [RCSB FAD chemical component](https://www.rcsb.org/ligand/FAD) has formula
C27 H33 N9 O15 P2, hence 53 heavy atoms. Comparing unique non-hydrogen atom names
in the tracked input against the [CCD atom table](https://files.rcsb.org/ligands/download/FAD.cif)
gives the following. These are counts of observed atoms, not minimized outputs.

| Input residue | Heavy atoms present | Missing from the 53-atom CCD set |
|---|---:|---:|
| A/FAD402 | 43 | 10 |
| B/FAD402 | 37 | 16 |
| C/FAD401 | 35 | 18 |

All three lack adenine atoms C2A, C4A, C5A, C6A, C8A, N1A, N3A, N6A, N7A,
and N9A. B/FAD402 additionally lacks C1B, C2B, C3B, O2B, O3B, and O4B;
C/FAD401 also lacks C4B and C5B. There are no unexpected heavy-atom names after
unquoting CIF values. This establishes incomplete FAD relative to the CCD; it
does not establish why atoms are absent, the intended redox/protonation state,
or an experimentally justified conformation for the missing groups.

Reproduce by reading `_chem_comp_atom.atom_id` and `type_symbol` from FAD.cif,
unquoting CIF strings, excluding H, and comparing to columns 13–16 of HETATM
records grouped by chain/resSeq/insertion code. Preserve the input checksum and
record the CCD retrieval date when repeating this comparison.

**Inference:** neither changing AMBER to CHARMM nor supplying an XML parameter
file solves the missing-coordinate problem. A repair must explicitly choose
component chemistry and rebuild/validate absent atoms, or intentionally retain
an incomplete model with stated limits, or exclude the component. Renaming the
fragment as a complete FAD template would conceal this decision. No cofactor
reconstruction or full-complex relaxation benchmark was performed for this note.

## Backend comparison

| Route | Available in dvbfixer 0.8.5? | Capability and limitation |
|---|---|---|
| OpenMM with native/additional templates | Yes | Joint energy minimization when all retained atom/bond graphs are parameterized; template coverage is required. |
| GAFF2 template generation | Yes, for eligible organic ligands | Extends OpenMM coverage; does not supply missing heavy atoms or validated metal/cofactor chemistry. |
| Protein-only strip-and-splice | Yes | Preserves otherwise unsupported material in output but does not optimize its interface jointly. |
| xtb GFN-FF / OpenBabel UFF or MMFF post-pass | Yes, optional | Runs after OpenMM; cannot rescue an earlier rejection. Heterogen-only mode freezes anchors/protein and limits interface accommodation. |
| Phenix/cctbx geometry regularization | No | Candidate dictionary/restraint-based route; requires valid component/link restraints and integration work. |
| A custom OpenMM geometry-restraint objective | No | Possible execution engine for the same proposed domain contract; chemistry dictionaries, weighting and validation would have to be implemented. |

OpenMM residue templates specify atoms and bonds and support template generators;
this is the mechanism behind parameter coverage, not automatic recognition of
all chemical states. See [OpenMM force-field construction](https://docs.openmm.org/7.7.0/userguide/application/05_creating_ffs.html).

Phenix documents geometry minimization with standard geometry restraints plus
optional secondary-structure, rotamer, Ramachandran, reference-position and
custom restraints. Its interface also exposes linkage and unknown-component
handling. This makes it a candidate for evaluation, not evidence that an
arbitrary incomplete cofactor will be repaired automatically. See
[Phenix geometry minimization](https://phenix-online.org/version_docs/dev-2222/reference/geometry_minimization.html).

GFN-FF is a force-field method exposed by xtb; it should not be described as a
quantum electronic-structure calculation merely because xtb also offers such
methods. Its applicability must be tested for the intended chemistry. See the
[GFN-FF documentation](https://xtb-docs.readthedocs.io/en/latest/gfnff.html).
The preparation choices `legacy` and `tleap-reduce` are separate from all these
minimization choices; the latter remains an opt-in pure-protein prep route.

## Proposed architecture (unimplemented)

Extend the existing [domain model](../domain-model.md) incrementally. A component
inventory would carry stable atom/residue identity, observed versus rebuilt
atoms, chemical graph, links, completeness, and provenance. A separate relaxation
request would express movable selections, coordinate tolerances, required
geometry restraints, and whether full MD parameterization is needed. Adapters
would report unsupported chemistry before optimization and return coordinates
through an explicit identity map with validation results.

Geometry regularization and MD parameterization should have distinct success
criteria. The former could optimize bonds, angles, planarity, chirality, steric
contacts, and reference restraints across the entire retained complex without
claiming MD-ready charges or energies. The latter still requires compatible
force-field parameters. This boundary is a design recommendation, not an
existing service or flag. The present DDD integration covers chain allocation
and ligand-routing policy; parsers and scientific pipelines retain substantial
responsibility.

Risks include wrong ligand bond orders, missing atoms, alternate conformers,
metal coordination mistaken for ordinary covalent bonds, unsupported redox
states, loss of insertion codes, and component graph changes across adapters.
Overstrong reference restraints can preserve clashes; weak ones can permit
drift. Geometry-only objectives do not establish binding energy or correct
protonation. A Phenix adapter would also require deployment/license review and
versioned dictionary provenance. Universal-force-field coverage is not proof
of accuracy for every linkage.

## Proposed benchmarks and acceptance criteria (not yet run)

1. Establish a parameterized OpenMM baseline and compare candidate geometry
   routes using identical observed coordinates and explicit chemistry choices.
   Use `--no-solvent` for development iterations. Separate any later explicit
   solvent evaluation from this baseline.
2. Cover a complete protein/organic ligand, an N-linked glycan tree, complete
   validated FAD/HEM models, the incomplete 8UCD fixture, antibody insertion
   codes, and multi-MODEL inputs. Unsupported cases should fail explicitly;
   do not count silent component removal as success.
3. Check atom inventory, stable identity, CONECT endpoints, covalent link lengths,
   zero D-Cα, ligand stereocentres, bond/angle outliers, planarity, severe clashes,
   and protein–ligand interface contacts before and after. Measure protein,
   ligand and glycan RMSD separately against the same reference frame, alongside
   maximum observed-atom displacement.
4. Record runtime, peak memory, tool/dictionary versions, parameters, seeds,
   failures and repeated-run variability. Do not compare raw energies across
   different objectives as if they shared a scale.
5. Require no unexplained atom loss or graph change and no new stereochemical
   inversion. Set numerical geometry/displacement thresholds before running,
   with chemistry-specific exceptions reviewed explicitly. Choose a default
   only after coverage and quality evidence, not merely successful execution.

Version 0.8.5 adds documentation and Split/viewer corrections. It makes no claim
that the proposed regularization architecture or these benchmarks are complete.
