# Rosetta repack: mechanism, ownership and local probes

Status: **source analysis plus bounded mechanism experiments**, not a replacement
backend or scientific accuracy qualification. This follows the
[48-mask comparison and native-contact control](protein-completion-unfiltered.md).
Continuation: a separate [atom-fixed scoring pilot](rosetta-fixed-scoring.md)
now evaluates retained generated-only candidates without repacking. It preserves
physical observations and improves three retrospective mask selections, but does
not make ordinary residue-level repack atom-preserving or change a steric gate.
The user authorized use of the existing local PyRosetta environment. No package
installation, environment repair, production routing or steric-policy change.

## Revisions and evidence

Primary source revision matches the installed **PyRosetta 2026.06**, Python 3.11.14:
[`1a56185c2592611dec4c9c75ddc9468cd2227c1f`](https://github.com/RosettaCommons/rosetta/tree/1a56185c2592611dec4c9c75ddc9468cd2227c1f).
The build identifies itself as native Apple Silicon, cxx11thread/serialization,
release dated 2026-01-30. Current upstream `main` was separately inspected at
`eddc15bbac60cf91302d3faf76895bd6d215cc50`; do not attribute its revision to local
execution. Primary evidence retains **32 source files** with URLs and SHA-256;
two unsuccessful guessed source filenames remain recorded as discovery errors.

[Portable source/probe evidence](rosetta-repack-mechanism-evidence.json) binds the
plan, harness, source manifest, snapshots, diagnostic poses, score weights and
contact controls. Local artifacts:

- `tmp/pdbfixer-research/rosetta-installed-source-01/`: matching source snapshots.
- `tmp/pdbfixer-research/rosetta-mechanism-plan-02/`: frozen inputs and 16-cell plan.
- `tmp/pdbfixer-research/rosetta-mechanism-run-02/`: all executed cells and logs.
- `tmp/pdbfixer-research/rosetta-contact-02/`: native and contracted pair controls.

Upstream sources remain local research material, **not vendored or redistributed**
in this package. PyRosetta uses its bundled empirical library internally; no
rotamer tables were exported or supplied to our own reader. This does not resolve
the separate authentic-data acquisition/licensing gate of our private generator.

## 1. Import already reconstructs atoms, and may repack

The observed source-to-pose path is:

```text
PDB -> StructFileRep -> PoseFromSFRBuilder
    -> map coordinates and missing-atom identities
    -> fill_missing_atoms from residue internal coordinates
    -> optional pack_missing_sidechains
    -> optional hydrogen optimization
```

[`PoseFromSFRBuilder.cc:1341–1387`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/io/pose_from_sfr/PoseFromSFRBuilder.cc#L1341-L1387)
initializes the missing mask and calls `Conformation::fill_missing_atoms`.
[`import_pose.cc:891–900`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/import_pose/import_pose.cc#L891-L900)
then conditionally packs missing sidechains and optimizes hydrogens. Consequently,
**disabling import-time packing does not leave absent heavy atoms absent**.

Important documentation trap:
[`pack_missing_sidechains.cc:47–83`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/pack_missing_sidechains.cc#L47-L83)
still describes rotamer trials, including in the Python docstring, but its actual
call is **`pack_rotamers`**, not `rotamer_trials`. The atom mask identifies which
**residues** become repackable; it is not an atom-level immobility mask.

## 2. PackerTask controls residue freedom, not missing-atom ownership

`PackerTask` specifies packable/designable/fixed residues and allowed residue
types. Repacking requires explicit sequence restriction: generic packing APIs can
otherwise allow design. Our probe uses `restrict_to_repacking()`, freezes every
non-target residue and enables `IncludeCurrent` for the selected target.

`IncludeCurrent` adds the current conformation as **one candidate**; it does not
force its selection or freeze its atom coordinates. For incomplete inputs, the
"current" conformation already contains import-generated atoms. It must never be
mistaken for a fully observed experimental sidechain.

## 3. Rotamer generation is backbone-conditioned and typed

[`RotamerSet_.cc:346–425`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/rotamer_set/RotamerSet_.cc#L346-L425)
selects a residue-specific rotamer library, requests extra samples, generates
candidates, bump-filters them and optionally includes the current rotamer.

The canonical Dunbrack path:

1. Read backbone-dependent variables and periodic interpolation bins.
2. Interpolate rotamer probabilities, means and widths; accumulate a configured
   probability mass while building candidates.
3. Expand requested extra-χ and proton-χ samples.
4. Construct a typed residue in the existing backbone frame and set its χ angles.

See
[`RotamericSingleResidueDunbrackLibrary.tmpl.hh:1085–1129`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh#L1085-L1129)
and
[`create_rotamers_from_chisets:1321–1341`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh#L1321-L1341).
This is not our three fixed staggered angles or a Cartesian atom minimizer.

The `fa_dun` term also contributes to scoring: roughly a backbone-dependent
negative log rotamer probability plus periodic χ-deviation penalties, with
implementation-specific interpolation/normalization corrections. It is not just
a prior used to enumerate candidates.

## 4. Backbone preservation does not imply sidechain preservation

[`Residue::place:979–1122`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/conformation/Residue.cc#L979-L1122)
orients the new typed residue and copies matching backbone coordinates. Other
sidechain coordinates come from the new rotamer, not the original atom records.

The installed flag is **`-packing:preserve_c_beta`**. The older application-doc
spelling `preserve_input_cb` is not accepted by this build. Its implementation
aligns the CA–CB bond direction with the source direction by rotating sidechain
coordinates; it does **not copy the observed CB position or bond length**. Distal
observations are not pinned. With missing CB, its source direction is itself an
import reconstruction, not a measured anchor.

Residue topology supplies stereochemical internal coordinates (e.g. VAL is typed
`L_AA`); Rosetta also explicitly supports mirrored/D residue types. Do not infer
a universal zero-D guarantee from "repack". Our independent Cα/THR-Cβ and
post-rounding checks must remain mandatory for any future adapter.

## 5. Interaction graph and stochastic joint selection

[`pack_rotamers_setup:192–245`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/pack_rotamers.cc#L192-L245)
builds the packer-neighbor graph and rotamer sets and initializes the interaction
graph. For pairwise-decomposable terms, the conceptual objective is:

```text
E(states) = constant + sum_i E_i(state_i) + sum_(i,j) E_ij(state_i, state_j)
```

One-body contributions include internal rotamer terms and fixed-context effects;
pair contributions encode interactions between variable residues. Factories can
use dense/precomputed, lazy or linear-memory graphs; specialized terms can require
more than this simple pairwise model. Do not assume every modern score term is
pairwise-decomposable or every graph eagerly stores every pair energy.

[`FixbbSimAnnealer.cc:174–243`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/annealer/FixbbSimAnnealer.cc#L174-L243)
proposes state substitutions, queries energy deltas, applies temperature-dependent
Metropolis acceptance, and retains the best visited assignment. Finally,
[`pack_rotamers_run:281–288`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/pack/pack_rotamers.cc#L281-L288)
**replaces each selected residue**, not just its originally missing atoms.
This is stochastic search, not MD, Cartesian minimization or a global-optimum
guarantee. Our small ≤81-combination controls can retain exact enumeration.

Bump filtering checks complete fixed neighbors but only the backbone of other
packable residues. Their mutual sidechain conflicts belong in joint optimization,
not premature rejection against their old sidechain placements.

## 6. Count-pair: three bonds is 1–4, and cp_half is 0.2

Pinned definitions:
[`CountPairCrossover3.hh:40–51`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/scoring/etable/count_pair/CountPairCrossover3.hh#L40-L51),
[`CountPairCrossover4.hh:40–51`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/scoring/etable/count_pair/CountPairCrossover4.hh#L40-L51)
and [`CountPairFunction.cc:28`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/scoring/etable/count_pair/CountPairFunction.cc#L28).

| Heavy-atom graph distance | Crossover 3 | Crossover 4 |
|---|---:|---:|
| 1 or 2 bonds | Excluded | Excluded |
| 3 bonds (1–4) | Factor **0.2** | Excluded |
| 4 bonds (1–5) | Full | Factor **0.2** |
| >4 bonds | Full | Full |

Names alone are misleading: `cp_half` is explicitly `0.2f`, not one half.
[`BaseEtableEnergy::determine_crossover_behavior:703–760`](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/source/src/core/scoring/etable/BaseEtableEnergy.tmpl.hh#L703-L760)
selects behavior according to energy-method mode, residue connectivity/types and
options. Classic intra-residue and xover4 scoring are distinct methods.

In installed `ref2015`, classic `fa_intra_rep` has weight **0.005**, whereas
inter-residue `fa_rep` has weight **0.55**. The native VAL39 CG1–C pair is three
bonds apart and is **not universally ignored**: classic intra-residue repulsion
counts it with 0.2, while the xover4 path excludes it. Do not transplant a lone
exclusion or weight without the surrounding torsion/statistical energy model.

Actual single-pair etable control, moving only CG1 along its C-relative vector:

| Distance | Raw pair repulsion | Derived classic weighted contribution |
|---|---:|---:|
| Native 2.762756 Å | 4.974661 | 0.004975 |
| 0.85×: 2.348343 Å | 41.181150 | 0.041181 |
| 0.65×: 1.795791 Å | 157.296204 | 0.157296 |

The last column is raw repulsion ×0.2×0.005, not the full residue/pose score,
an acceptance threshold, kcal/mol calibration or a justified policy for our tool.
Contracted controls violate local geometry; other energetic/geometric terms
matter. Python's count-pair binding returns its boolean but not mutated C++
weight/path-reference outputs; the factor comes from pinned source inspection.

## 7. Sixteen local ownership probes

Four **retrospective mechanism cases**, not an unfiltered accuracy cohort:
1TIM VAL39 missing CG1 or whole CB/CG1/CG2, and 1EMV THR11 missing CG2 or whole
CB/OG1/CG2. Cross import-time packing off/on with preserve-Cβ off/on: **16 cells**,
each a fresh single-threaded process, seed 7, explicit `ref2015`, H optimization
off. Each records loaded coordinates before explicit target-only repacking.
An additional no-coordinate preflight confirmed that import's default
`get_score_function()` is named `ref2015` and has the same weights as the explicit
score function in this environment; its result is retained with the run.

All 16 executed. Without automatic packing, import preserved every observed heavy
atom in all eight cells, while filling missing atoms. With automatic packing,
the four singleton cells already moved observed target sidechain atoms at import.
Explicit repacking moved observations in **all eight singleton cells**; all eight
whole-sidechain cells preserved observations. Backbone coordinates and residue
identities remained unchanged in every case. Independent target RDKit CIP checks
passed 16/16 (CA=S; additionally THR CB=R), not a global stereo qualification.

Explicit repack after import with packing disabled:

| Singleton | Preserve Cβ | Observed CB movement | Observed partner movement |
|---|---|---:|---:|
| VAL39 missing CG1 | Off | 0.040256 Å | CG2 0.393976 Å |
| VAL39 missing CG1 | On | 0.001234 Å | CG2 0.435373 Å |
| THR11 missing CG2 | Off | 0.115771 Å | OG1 0.246793 Å |
| THR11 missing CG2 | On | 0.015855 Å | OG1 0.118334 Å |

The separate diagnostic candidate set had two VAL candidates and 19 THR candidates
per cell. THR counts include proton-χ variants: only two distinct sampled χ1
values, not 19 independent heavy-sidechain conformers. With import packing off,
only **one candidate per singleton set** preserved original target observations:
the current, already import-completed conformation. IncludeCurrent therefore does
not solve the fixed-partner constraint. Diagnostic sets are not captured annealer
traces; the mover builds its own sets, and no exact-repeat claim is made.

All raw snapshots/diagnostic outputs are retained. No post-hoc observation
restoration disguises an incompatible proposal; RMSD fields are descriptive only.
Whole-sidechain preservation here does not establish native accuracy or general
compatibility across residue types, variants, complexes or source defects.

## Implications for DVBFixer

**Adopt the separation of responsibilities, not the default ownership contract.**

1. Keep source admission, observed anchors, L/THR stereochemistry and rounded
   readback as independent hard constraints. Define feasible candidates first.
2. Introduce topology-aware steric diagnostics before specifying a new gate.
   Rosetta count-pair is evidence against a universal VDW cutoff, not permission
   to suppress all intra-residue or 1–4 clashes.
3. Evaluate richer backbone-dependent candidates and a separately versioned
   score for fixed-context packing, then cache new–new interactions for cohorts.
4. For partial sidechains, constrain generation/closure to observed atoms. Do not
   repack whole residues and splice old atoms back: scoring then described a
   different geometry, and the splice can reintroduce clashes or wrong stereo.
5. Keep ordinary Rosetta whole-sidechain repacking as a **different-contract**
   comparator until an atom-preserving adapter is independently demonstrated.
6. Preserve exact enumeration for small controls; add seeded heuristic search
   only when candidate/cohort growth makes it necessary, with exact small-case
   optima as tests. A lower Rosetta score is not automatically better held-out RMSD.

The pinned [Rosetta license](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/LICENSE.md)
and [PyRosetta license](https://github.com/RosettaCommons/rosetta/blob/1a56185c2592611dec4c9c75ddc9468cd2227c1f/LICENSE.PyRosetta.md)
restrict commercial use and redistribution. Existing local installation is not
authorization to bundle code, databases or derived library exports. No Rosetta
dependency is added to Standard/Desktop and no licensing clearance is inferred.

## Reproduction and limits

`scripts/research_rosetta_repack.py plan` freezes four existing comparison masks
without importing PyRosetta. Each `probe` runs in a fresh process and verifies the
plan's input/harness digests before initialization. Use exclusive new paths:

```bash
python scripts/research_rosetta_repack.py plan --frozen tmp/pdbfixer-research/unfiltered-plan-01 --output tmp/pdbfixer-research/new-rosetta-plan
# Run with the licensed local PyRosetta interpreter; repeat all declared cells.
python scripts/research_rosetta_repack.py probe --plan tmp/pdbfixer-research/new-rosetta-plan --case 1TIM-VAL39-CG1 --output tmp/pdbfixer-research/new-rosetta-cell --preserve-cb
```

`scripts/summarize_rosetta_mechanism.py` independently rechecks saved coordinates,
identities, candidate preservation, hashes and target CIP in the normal scientific
environment, without importing PyRosetta or generating coordinates. The executed
contact helper and full-cell orchestrator remain beside the retained artifacts as
`tmp/pdbfixer-research/rosetta-contact-probe.py` and `run-rosetta-mechanism.py`.
The rejected old flag and an initial vector-binding error are retained separately
as protocol blockers; neither was counted as a scientific failure or silently
overwritten. No broader families/seeds, joint-packing performance benchmark,
exact-repeat guarantee, external packing calibration or complete legacy-pipeline
comparison was executed.

Final verification: **559 focused tests passed in 25.15 s**, targeted Ruff and
seventeen-file private-package mypy passed. Independent snapshot-summary replay
is byte-identical. All previous 48 comparison outcomes, 32 native audits and
six original/eight mixed evidence records still pass read-only validation with
unchanged frozen source/artifact hashes. These are not full-repository tests.
