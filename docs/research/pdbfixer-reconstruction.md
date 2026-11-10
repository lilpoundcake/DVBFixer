# PDBFixer reconstruction: causal findings and a constrained alternative

Date: 2026-11-02. **Private research, not a production backend change.**
Baseline: `main` at `4a328a6`; branch `research/pdbfixer-reconstruction`.
The primary question is loss of amino-acid stereochemistry during reconstruction.
Hydrogen geometry is a separate secondary investigation.

**Implementation priority after user review:** consult the
[actual-usage, chemistry and replacement specification](pdbfixer-usage-and-replacement-scope.md)
before extending the prototype. Coordinate generation, discovery and complete
PDBFixer dependency removal have different functional/test requirements.

## Results in brief

1. **An initialization defect is reproduced:** PDBFixer's `_overlayPoints`
   allows determinant −1 transforms. In the declared 1UBQ gap, THR12 and ILE13
   are D **before any minimization or dynamics**. Complete sidechain deletions
   at LYS27 and ALA28 independently reproduce D initialization and final output.
2. **A separate optimization defect is reproduced:** forcing a proper rotation
   makes every gap centre initially L, but ILE13 becomes D during the first
   minimization. No Langevin fallback runs in these controls. Adding the existing
   signed-volume wall prevents that inversion in this case.
3. **A determinant fix is not sufficient:** in partial residues, template atoms
   are fitted but observed atoms are kept. The resulting hybrid is not a rigid
   transform of the template. Proper all-anchor fitting still initializes D THR7
   and ILE13. Proper fitting restricted to N/CA/C plus the wall avoids final D in
   all six declared sidechain controls. This is not broad qualification.
4. **An independent Cβ primitive is feasible:** solving bond-length/angle
   constraints with an explicitly L signed-volume branch gives positive volumes
   in six withheld-Cβ controls, with 0.011–0.237 Å deposited-Cβ errors. This
   constructs one atom, not complete sidechains or loops.
5. **The SER hydrogen collapse has a concrete force-level cause in the tested
   fixture:** a close fixed OG/OXT environment, attractive HG/OXT Coulomb term,
   and zero HG Lennard-Jones epsilon. Removing only that pair's electrostatics
   recovers ~0.964 Å OG–HG in all five declared seeds. Increasing 50 iterations
   to 500 alone does not repair it. This does not establish the cause of every
   hydrogen defect or excuse the source's already abnormal heavy geometry.

Neither seeds nor post-hoc reflection are explanations of these mechanisms.
None of these findings authorizes removing production safeguards.

## Pinned evidence and scope

The [portable evidence](pdbfixer-reconstruction-evidence.json) contains source and
resource SHA-256s, stage volumes, determinant logs, first negative minimizer
callbacks, final readback audits and all declared outcomes. The
[live handoff](../plans/pdbfixer-reconstruction-handoff.md) owns continuation.

Execution: Python 3.12.16, native Apple arm64, **OpenMM Reference**,
PDBFixer 1.12.0 and OpenMM 8.5.2 (`8.5.2.dev-36a30cb`). Installed source bytes
were compared with pinned upstream files and **match exactly**:

- PDBFixer commit `94cfa4c0ca551cdc5f13320f9a658efd59f2b881` (`v1.12`),
  `pdbfixer/pdbfixer.py`, SHA `a52e39125e3e0e2d1b58b83841fb05b8669e0bcdd73461e1b0ed1d5a8e8261b9`.
- OpenMM commit `36a30cbca54e727b216b606f3c011b67201eb8b4` (`8.5.2`),
  `wrappers/python/openmm/app/modeller.py`, SHA
  `376d8c8c06d1cfed4e40da20d77a2094099e872bda783b9c9bfc03a67305e813`.
- `soft.xml`, available PDBFixer PDB/XML assets, OpenMM hydrogen definitions and
  Amber14/ff14SB definitions are separately digest-bound in the evidence.

The inspected general Standard prefix currently lacks `soft.xml`; no files there
were changed. Runs instead used the intact isolated cc94 v3 CPU refinement
runtime. This is not evidence about why package resources disappeared.

All source inputs come from tracked fixtures or an existing in-test fixture:

- `tests/fixtures/protein_context/1UBQ.pdb`, explicitly projected to protein chain
  A residues 1–72; original fixture bytes and occupancies are not changed.
- Sidechain-only control: remove LYS27 CG/CD/CE/NZ, retaining observed CB.
- Gap: delete A10–A13 and explicitly declare GLY/LYS/THR/ILE at that location.
  There is no invented SEQRES or second alignment algorithm. This is bare
  PDBFixer insertion, **not** the production MODELLER gap-reconstruction baseline.
- Six-site sidechain cohort, declared before execution: THR7, ILE13, SER20,
  LYS27, ALA28 and ALA46. Remove all sidechain heavy atoms, retaining N/CA/C/O.
- Hydrogen fixture: `_SER_INPUT_HG_MISSING_OXT_PRESENT` extracted by AST from
  `tests/test_prepare_broken_geom.py`, rather than a copied undocumented input.

Only seed 1 is needed for the heavy controls: no fallback dynamics occurs.
Hydrogen seeds are Python `random` seeds 1, 7, 19, 42 and 101. Cross-platform
bit-identical results are not inferred. Crystal context, hydrogens and secondary
THR/ILE stereocentres are not comprehensively validated by this heavy-atom study.

## How PDBFixer actually reconstructs atoms

Line ranges refer to the pinned `pdbfixer.py`, not moving upstream `main`.

| Stage | Symbols / lines | Meaning and limitations |
|---|---|---|
| Whole-residue detection | `findMissingResidues`, 874–944 | Uses sequence records, residue numbers/insertion codes and matching offsets; absence of SEQRES is not proof of a complete chain. DVBfixer's shared sequence placement remains authoritative. |
| Atom detection | `findMissingAtoms`, 1108–1205 | Compares residue atom names with templates and separately records terminal atoms. Its residue-object keys do not survive topology replacement. |
| Partial residue construction | `_addAtomsToTopology`, 579–738 | Overlays one residue template on matched observed atoms; inserts only missing coordinates, leaving observations untouched. This is not environment-aware rotamer packing. |
| Whole-gap construction | `_addMissingResiduesToChain`, 744–806 | Places template copies along an arc between residue centres; attempts trans orientation. It does not initially solve exact peptide closure. |
| Template superposition | `_overlayPoints`, 214–257 | SVD with no determinant correction; returns an orthogonal transform, potentially a reflection, rather than necessarily a proper rotation. |
| Temporary objective | `_createForceField`, 1582–1674 | Uses heavy-atom `soft.xml`, not a fully protonated production Amber system. |
| Optimization and transfer | `addMissingAtoms`, 1207–1331 | Mass-zero observations; minimize new atoms; possibly clash-escape MD; rebuild full topology and copy optimized new coordinates into it. |

The soft objective contains bonds, angles, torsions, planar impropers and bounded
repulsion `C/((r/0.2)^4+1)`. It has **no explicit L-Cα signed-volume guard**.
Repulsion is finite at zero distance; small energy or successful minimization
cannot certify geometry. If a new atom remains within 0.13 nm of a non-excluded
neighbour, the fallback uses 300 K Langevin dynamics, 5 fs, up to 10×200 steps,
ramps the repulsive coefficient and minimizes again. Same-residue pairs and
direct bonds are excluded from that trigger. Seeds control the stochastic
fallback, not the correctness of initial stereochemistry or the objective.

This agrees with the historical upstream explanation in
[PDBFixer issue 145](https://github.com/openmm/pdbfixer/issues/145), but our
specific conclusions below come from current pinned, locally executed controls.

## Three distinct stereochemical failure mechanisms

Define `V = ((N−CA) × (C−CA)) · (CB−CA)` in nm³. Here positive means L for the
supported canonical residue naming; glycine has no Cβ centre. Undefined,
degenerate and almost planar centres require separate rejection, not an L label.

### 1. Reflection during template placement

The original SVD computes `R = (U @ Vt).T`, without enforcing `det(R)=+1`.
Replacing it privately with `R = (U @ diag(1,1,det(U@Vt)) @ Vt).T` removes the
reflections in the gap experiment. Both THR and ILE templates themselves have
positive signed volume; the D configuration is not inherent to those assets.

| Gap A10–A13 control | D at construction | D after minimization | D after PDB readback |
|---|---|---|---|
| Original | THR12, ILE13 | LYS11, THR12, ILE13 | same three |
| Proper rotation only | none | ILE13 | ILE13 |
| Proper rotation + signed-volume wall | none | none | none |

For the original gap, initial THR12/ILE13 volumes are approximately
−0.002457/−0.002487 nm³. Final LYS11/THR12/ILE13 volumes are
−0.002073/−0.002848/−0.002706 nm³. Observed coordinates move **0.0 nm**.

### 2. Inversion during optimization, independently of initialization

In the proper-only gap arm, ILE13 is initially L and becomes negative by
minimizer callback **iteration 0**; final volume is −0.002878 nm³. A callback is
an optimizer report, not a continuous trajectory or every internal line-search
evaluation. We therefore locate the first failing reported optimizer stage,
not an exact continuous-time planar crossing.

The same initial proposal plus the existing
`ffutils.geometry.build_ca_chirality_force` wall produces no final D centres.
This isolates optimization as another failure source. No hydrogen placement,
MD fallback, source-Cα movement or downstream full-system minimization is needed
to produce this example.

**Important negative result:** the wall-protected gap still has **nine findings**
from independent bond-length/backbone-angle/peptide-omega readback checks,
versus five for the unguarded arms. Preventing D is necessary, not sufficient for
a useful loop. These are rejected research structures, not accepted repairs.

### 3. Hybrid fixed-anchor geometry can defeat a proper rotation

For partially observed residues, the original N/CA/C/O coordinates are retained;
only missing sidechain positions come from the fitted template. A proper rotation
preserves chirality of an **entire transformed template**, but not necessarily
of this hybrid when fitted anchors and retained anchors disagree.

In the six-site cohort, original fitting ends with D at **LYS27 and ALA28**.
Proper all-anchor fitting removes those two failures but ends with D at **THR7
and ILE13**, already D at construction. The wall does not repair these because
it intentionally guards only centres that start L. Guarding a wrong starting
sign or silently calling it L would be a scientific error.

The declared additional intervention fits the proper template frame to **N/CA/C
only**, excluding O from the least-squares fit, then applies the same wall.
All six final centres are L, observations remain unchanged, and the three
limited readback geometry checks produce zero findings. This is one small
same-protein causal cohort, **not a demonstrated universal solution**, rotamer
accuracy claim or full stereochemical/steric qualification.

## Why the tested hydrogens collapse

OpenMM `Modeller.addHydrogens`, pinned `modeller.py` lines 753–1107, chooses
variants and expected H names, builds provisional positions from a sum of
directions away from bonded atoms plus random positive-component perturbations,
then minimizes newly added H only. Existing heavy atoms and existing H are
immobile. Automatic variant selection retains existing hydrogens; explicit
variant selection can remove extra H. This path is distinct from PDBFixer
heavy-atom rebuilding and from Salilab MODELLER.

With a supplied force field, OpenMM uses its real bonded/nonbonded objective,
with a nominal 50-iteration minimization. Without one, it uses a different
repulsive/distance objective; that is not equivalent force-field preparation.

The tested source has **OG–OXT = 1.67804 Å before adding HG**, already an abnormal
nonbonded heavy-atom separation. Amber14 parameters in the actual constructed
system are `q(HG)=+0.4474 e`, `q(OXT)=−0.8132 e`, and `epsilon(HG)=0`. The pair
is not a topology exception. Its LJ mixed epsilon is zero; attractive Coulomb
energy has no opposing pairwise LJ core. Its singular attraction can overwhelm
the finite cost of stretching OG–HG while fixed heavy atoms cannot separate.
For seed 19 the final HG–OXT separation is about **0.000013 Å** and the total
reported energy about **−39.6 million kJ/mol**: a lower energy is not a sane H.

| Intervention (five declared seeds) | OG–HG result | Interpretation |
|---|---|---|
| Original Amber objective, 50 iterations | 1.568–2.472 Å | All five detected as misplaced |
| Same, 500 iterations | identical to 50 | Iteration budget alone does not explain/repair this fixture |
| Strip all H first, then original placement | 1.646–1.679 Å | Stripping alone does not solve it |
| No force field | 1.186–1.190 Å | Avoids gross detector failures, but stretched; not a validated replacement |
| Same initial coordinates, remove only HG/OXT Coulomb | 0.96427–0.96436 Å | Pair-specific causal intervention, all five |
| Same initial coordinates, add pair repulsive core | 0.96424–0.96429 Å | Prevents singular collapse; still very short HG–OXT ~1.128 Å |

These interventions are diagnostics, **not recommended force-field edits**.
The result is more specific than “CSER template bug” and narrower than “all
OpenMM hydrogens are wrong.” Other failures, including sibling-H angles, frozen
bad existing H, variant/bond mismatches and emergency reflection strain, need
their own controls. Existing hydrogen repair remains enabled.

## Can we build our own geometrically disciplined atom rebuilder?

**Yes for a bounded canonical-protein completion service; no justification for
a universal reconstruction replacement yet.** Missing distal sidechains are
conformationally ambiguous. Valid geometry is achievable without knowing the
experimental conformation, but no atom generator can recover absent experimental
information merely by minimizing an energy.

### Executed minimal independent primitive

`scripts/research_chiral_completion.py::construct_l_cb` does not overlay a
PDBFixer template. Let `u` and `v` be normalized CA→N and CA→C. For declared
CA–CB length `l` and angles `alpha`, `beta`, solve the two projection constraints
for an in-plane vector `p`. The solutions are

```text
CB = CA + p ± sqrt(l² − |p|²) * normalize(u × v)
```

Choose **only the positive-volume branch** for this canonical L-protein scope.
Reject coincident/collinear/nonfinite anchors and incompatible angle constraints.
No observation is moved. The research defaults are `l=1.53 Å`, angles
110.5°/110.1°: declared illustrative geometry, not a complete calibrated
residue/backbone-dependent parameter library.

At the same six withheld-Cβ sites: all volumes positive, CA–CB exactly 1.53 Å,
deposited-Cβ errors **0.054, 0.237, 0.052, 0.011, 0.040, 0.082 Å** respectively.
Reference CB is used only for scoring, never the constructor. Tests verify
lengths, angles, L sign, proper-motion equivariance over 100 rotations and
explicit refusal of impossible anchors. This is a feasibility control, **not**
six successful full sidechain predictions or a comparison to a full legacy run.

### Proposed service architecture

1. **Resolve chemistry and exact identity first.** Canonical residues initially;
   preserve chain case, insertion codes, explicit links and atom names. Unknown
   chemistry, ambiguous variants/altlocs and incompatible fixed atoms refuse.
2. **Classify the missing subgraph.** One terminal atom, branched sidechain,
   aromatic ring, missing backbone and whole gap need different solvers.
3. **Build in stereochemically constrained internal coordinates.** Explicit L
   Cα branch; correct secondary THR/ILE stereochemistry; dictionary bond lengths,
   angles, planar rings and closures. Preserve every observed atom, not just CA.
4. **Sample actual conformational freedom.** Enumerate backbone-dependent
   rotamers and local chi perturbations consistent with all retained sidechain
   atoms. Partial branches with multiple fixed anchors require constrained closure,
   not rigid template overwrite. PRO and rings require dedicated closure.
5. **Rank in retained context.** Reject covalent/stereo/clash violations first;
   then rank rotamer likelihood, sterics and chemistry-aware interactions, keeping
   scores separate and retaining ambiguous alternatives. Coupled neighbouring
   missing sidechains require packing, not independent greedy placement.
6. **Refine without escaping the feasible stereochemical set.** Prefer torsion/
   internal-coordinate refinement. If Cartesian optimization is needed, retain
   signed-volume barriers and hard postchecks. A finite penalty is never itself
   a mathematical guarantee; reject any final forbidden/nearplanar centre.
7. **Independently audit and publish.** Exact generated mask; observations byte-
   identical; expected graph, all stereo centres, junctions, clashes, planarity,
   rotamers and completeness checked. Use canonical `DatRecord` provenance and
   atomic publication. No output on unresolved geometry, even if energy is low.

Whole missing residues remain a **separate loop-closure/modeling task**. An
internal-coordinate chain must satisfy both anchors, peptide length/angles,
omega and Ramachandran constraints, then pack sidechains. Fragment/CCD/KIC
methods are candidates, not an implemented feature. MODELLER remains the proper
production gap comparator; the failed arc insertion is not its substitute.

### Relevant precedents and implementation decision

- [SCWRL4, Krivov et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2885146/)
  uses rotamer/subrotamer enumeration, environmental scoring and coupled packing.
  It illustrates the algorithm class absent from PDBFixer's single-template fit;
  its published accuracy is not evidence about our prototype.
- [Dunbrack 2010 library](https://dunbrack.fccc.edu/lab/bbdep2010) supplies
  backbone-dependent probabilities, chi means and variances. The provider states
  CC BY 4.0 for the library. Pin the exact artifact, attribution and terms before
  use; this does not license SCWRL binaries/source or other third-party code.
- Existing [backend research](reconstruction-and-modeling-backends.md) already
  separates chemistry authority, coordinate materialization and MD parameters.
  Retain that boundary rather than creating a second discovery/chemistry stack.

**Recommendation:** build a private canonical sidechain/internal-coordinate
prototype first, not another unconstrained Cartesian “sprinkle and minimize”
engine. Use the L-Cβ primitive as a tested building block, then add residue-aware
geometry, rotamer data and observation-constrained packing. Benchmark separately
against seeded PDBFixer, `tleap` where applicable, and full legacy for gaps.
Do not change defaults or promise conformation recovery from these results.

## Verification and limitations

- 30 executed heavy rebuilding arms (6 initial controls + 24 six-site controls),
  20 hydrogen placement arms and 15 pair-force re-minimizations; six independent
  Cβ constructions. Re-minimizations reuse saved proposals, not new samplings.
- Independent PDB readback: **0.0 nm observation movement in all 30 heavy arms**;
  failures and the guard-protected gap's remaining geometry defects are retained.
- Focused new/existing seeded-rebuild/chirality tests: **20 passed**; targeted
  Ruff check passed. No full-suite or CPU/GPU/production acceptance is claimed.
- `run-01` contains a research-instrumentation failure (forcefield argument
  introspection), not executed heavy science; `run-02` corrected it; `run-03`
  added optimizer callbacks. Earlier attempts remain intact and are not counted
  as independent confirmations. `sidechains-02` includes the predeclared
  backbone-frame intervention after the first cohort exposed its need.
- Only Cα stereochemistry is comprehensively traced here. Whole-complex sterics,
  THR/ILE secondary chirality, broad source cohorts, production platforms,
  clinical/noncanonical chemistry and complete own-builder accuracy remain open.
- Full coordinates, source snapshots, input derivatives, systems, callbacks and
  logs are retained locally under `tmp/pdbfixer-research/`. They are ignored by
  Git; the portable record and reproduction scripts carry the transferable
  findings, not a claim that raw files transfer to another host automatically.
