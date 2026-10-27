# Corrected Apple multi-gap regression: supplied 5J13 input

Date: 2026-10-22. This is a development regression, **not a qualified acceptance
cohort or a successful structure preparation**. No validation threshold was
relaxed and no candidate was selected using reference coordinates.

## Input and confirmed failure

For the source-code comparison with PATCHR TCD/LRD and the separate 9CQD
per-step reinjection experiment, see
[PATCHR versus dvbfixer gap reconstruction](patchr-gap-reconstruction-comparison.md).

### Related 9CQD metadata regression (2026-10-23)

The separate user-supplied `diagnostics.tar.gz` contains a 9CQD derivative,
not another 5J13 run. Its refined coordinates contain the requested M/W/X
27/214/221 residues, but inherited SEQRES declares 49/220/261. The archived
FASTA contains 27/220/221; explicit terminal exclusion reduces W to 214 in
the frozen request. Thus the request, not the raw FASTA, defines output SEQRES.
Inherited foreign-chain DBREF/HELIX/SHEET/CISPEP records are also stale.

Diffusion materialization now reconciles dependent metadata, preserves valid
links and experimental headers, and records changes in a `METADATA_SYNC`
remark. Invalid ranges are omitted rather than guessing DBREF mappings;
relational sheets, sites and biological assemblies are retained only whole.
The digest-bound input remains available for original deposited annotations.
Independent validation rejects supplied inconsistent metadata via
`pdb-metadata-consistency`; absent optional annotations are allowed.

Re-auditing the archived refined candidate after metadata-only reconciliation
removes all SEQRES/ATOM-count warnings with unchanged coordinate, ANISOU,
CONECT and SSBOND records. It **still fails** amide planarity at W28→W29:
ω = −132.2902°, 47.7098° from trans (30° tolerance). Metadata repair does not
repair or qualify these coordinates. No geometry threshold was relaxed.
Only the archived source PDB and FASTA are tracked under
`tests/fixtures/protein_context/9cqd_split.*`, not generated coordinates.

### Original 5J13 input

The user supplied `test/errors_diffusion/selected.pdb`, now retained as
`tests/fixtures/protein_context/5j13_selected.pdb` with SHA-256
`63154130ff3de95ab1a135154370ac91e07d8b776bb16067096c0c62c91e3106`.
It is derived from 5J13 after CIF normalization and alternate-state selection,
not an untouched deposited structure. There are 10 missing regions / 124 target
residues across chains A/B/C, including long terminal tails and partial observed
sidechains. Target sequence lengths total 655 residues.

The installed historical `protpardelle-1c-mps` launcher was implicitly chosen
because the corrected launcher was absent. The raw candidate had 36 D-sign Cα
centres, all generated. These were predominantly malformed geometries rather
than well-formed D-amino acids: 34/36 also had a CA–N/C/CB bond outside the broad
1.2–1.8 Å diagnostic interval. Across all generated residues, 122/124 had such a
bond anomaly. This diagnostic interval is not a replacement acceptance threshold.

Per-gap refinement then treated other generated gaps as fixed. The first call
refused their D centres before minimization. Its precise exception was discarded
in favor of the generic postprocessing-failure code.

## Same-sequence monomer control

Chain C was sampled without partner chains, retaining the same four gaps and
effective seed 9 (the original third chain invocation used seed 7 + 2).

| Raw chain C result | Historical partner context | Monomer ablation |
|---|---:|---:|
| D-sign Cα centres | 17 | 0 |
| Generated residues with the broad CA-bond anomaly | 52/54 | 0/54 |

This implicates the historical context representation in this run, not an
intrinsic inability of MPS or cc89 to produce local L geometry. It does not
separately identify chain phase, position features or context composition as
the sole cause. Removing partners is a scientific ablation, not the production
solution. Subsequent CPU minimization could invert generated LEU13; another
repeat exposed comparison against temporary hydrogen placement. CPU coordinates
are not claimed bit-identical across reruns.

## Implemented correction

- New requests require the corrected launcher; historical replay selects its
  explicit profile directly, without an implicit downgrade.
- Existing corrected MPS sampling is reused per gap under the request option
  `cc89_gap_execution=independent-local-v1`. Each source-only local view retains
  physical positions and neutral chain phase, a 512-residue engine budget, actual
  sampled atom order, context evidence and 1,000 first-layer RoPE calls.
- Other unobserved regions are omitted from the view, never promoted to fake
  context. Partial sidechains contribute only observed atom masks; unrequested
  network predictions for their absent atoms are discarded. Backbone must exist.
- Generated coordinates alone are merged. All generated regions share one
  mutable refinement mask; original heavy atoms remain fixed. Temporary FF
  hydrogens are rebuilt without changing publication's source records.
- Refinement revision v7 no longer constrains raw sampled bond lengths and only
  ranks finite all-L minima. Expected scientific failures retain the raw candidate
  and exact reason. Private workspaces survive runner/validation failure.

Old single-gap request bytes and historical acceptance records are not rewritten.
The earlier v6 refinement cohort does not establish acceptance for this wider
scope or the revised minimization.

## Real corrected run: remaining rejection

Seed 7, all ten gaps, no terminal cropping, corrected MPS sampling and CPU
refinement completed. The supplied input digest remained unchanged. Artifact
directory during development: `test/errors_diffusion/corrected-second/`; private
workspaces and generated outputs are intentionally not tracked fixtures.

| Independent metric | Result |
|---|---:|
| Missing/wrong generated residues and heavy atoms | 0 |
| Maximum fixed-heavy displacement | 0 Å |
| Detectable D-Cα | 0 |
| Generated residues with the broad CA-bond anomaly | 0/124 |
| Non-planar generated/junction amides | 0 |
| Sidechain χ1/χ2 outliers | 0 |
| Bond-length / bond-angle errors | 2 / 1 |
| Ramachandran outliers | 2 |
| Severe steric overlaps involving generated regions | 133 |
| Largest junction C–N distance | 4.811 Å, C173 → C174 |

The bundle correctly remains `validation_failed`; ZBS must not consume it.
Normal intra-residue bonds do not prove valid placement relative to the complex.

An additional 30-iteration refinement diagnostic measured initial/final energies
of approximately 9.84e25 / 2.52e16 kJ/mol and moved only one generated atom by more
than 0.001 Å. The raw merged candidate contains heavy-atom overlaps as short as
0.148 Å. Thus a finite all-L FF minimum attempt is not evidence of adequate
relaxation. The input audit below establishes an upstream mask error that must
be corrected before attributing the remaining defects to sampling or attempting
clash-safe initialization / staged refinement. Larger tolerances, source drift,
tail deletion or relabeling this rejected candidate are not remedies.

## Follow-up: the C168–173 mask is incorrectly placed

A passive rerun of gap 7 saved actual sampling tensors and predicted coordinates
before source restoration. The supplied 66-residue model sequence exactly matches
the declared local view, and the original raw candidate reproduces. Predicted
C173–N174 is 1.553 Å before restoration, but becomes 4.811 Å after rigid fitting
and restoring the observed anchor, whose N position differs by 4.977 Å.

The source explains why this candidate cannot be trusted as a valid mask:

| Original chain C boundary | C–N (Å) | Interpretation |
|---|---:|---|
| CYS139 → SER140 | 1.331 | Existing continuous peptide |
| SER140 → THR147 | 12.313 | Actual missing region, absent numbers 141–146 |

The ambiguous sequence placement maps observed SER140 to reference position 174
instead of 168, putting a six-residue insertion *before* that SER and breaking
the existing CYS–SER peptide. Source numbering, the surrounding unambiguous
placement offset (+28), and actual peptide geometry instead support missing
reference residues **169–174 (`RSTSES`)**, between observed SER140 and THR147.
The current request incorrectly generates **168–173 (`SRSTSE`)**.

This is not a malformed original PDB or proof of a sequence-tensor decoding bug.
It is incorrect ambiguity resolution upstream of renumbering and diffusion mask
construction. A deterministic sequence-only optimum is insufficient here. Fix
shared placement tie-breaking with observed numbering/connectivity evidence and
rebuild the request before further model-quality or refinement experiments.

## Shared placement correction

The affine aligner now ranks equal sequence optima lexicographically: first
minimize insertions into observed peptide bonds, then minimize insertions into
consecutive author numbering, then retain the historical deterministic tie rule.
It never changes primary sequence scores and never assumes an absolute numbering
offset. Missing, duplicate or altloc backbone atoms supply no geometry evidence;
case-sensitive chains, insertion codes and explicit TER boundaries are respected.
Exhaustive small-alignment tests verify the ranking independently of the DP.

`renumber`, direct diffusion placement, terminal trimming and model residue
restoration supply available evidence through this one service. The real source
now maps original SER140 to reference 168 and selects `169–174 / RSTSES`, between
original SER140 and THR147, both before and after renumbering. All ten gaps and
124 generated residues are retained; the source file is unchanged.

A seed-7 sampling-only rerun on the corrected single-gap view measured:

| Boundary after source restoration | Distance (Å) |
|---|---:|
| Existing C167–N168 | 1.3312 |
| New left C168–N169 | 1.1307 |
| New right C174–N175 | 1.2798 |

This removes the false 4.811 Å right-junction failure. It is a raw single-gap
diagnostic, not acceptance of the whole reconstruction or a claim that every raw
bond length is ideal.

## Complete ten-gap rerun after placement correction

`test/errors_diffusion/alignment-fixed-full/` records a new seed-7 run on the
freshly renumbered source, using the installed corrected profile without PATH
overrides. All ten gaps / 124 residues are retained. Sampling and refinement
complete; the bundle remains correctly rejected.

| Independent metric | Corrected-placement result |
|---|---:|
| Fixed-heavy displacement | 0 Å |
| Detectable D-Cα | 0 |
| Generated/junction bond-length errors | 0 |
| Generated/junction bond-angle errors | 0 |
| Non-planar amides / sidechain χ1/χ2 outliers | 0 / 0 |
| Junction C–N range | 1.1307–1.5428 Å |
| Ramachandran outliers | 1 |
| Severe overlaps | 138 |

The false 4.811 Å junction is gone. Remaining connectivity failures are outside
the unchanged 1.20–1.45 Å junction gate: A56→A57 (1.5428), A144→A145 (1.1487),
C168→C169 (1.1307) and C252→C253 (1.1536). These are residual sampling/refinement
quality failures, not justification to relabel the wrong mask as valid or relax
the gate. The next stage remains clash-safe, geometry-preserving refinement and
independent revalidation; this run does not establish a usable prepared model.

Software verification after the alignment change: 1,407 fast Python tests pass,
two optional-dependency skips and 31 slow deselections; Ruff and strict API-surface
mypy pass. Synthetic tests independently enumerate all small placements to check
lexicographic ranking, and the real regression covers original, freshly renumbered
and misleading legacy numbering while preserving the supplied input bytes.

## Captured refinement-system audit

The corrected-placement raw candidate was reconstructed into the exact temporary
FF system under Reference, with seed 7. The private snapshot is retained in
`test/errors_diffusion/refinement-system-audit/` (XML system, full-precision
positions/forces, mobility mask and temporary PDB). This is diagnostic evidence,
not a published candidate or acceptance result.

- Initial potential: approximately `9.854e25 kJ/mol`.
- Largest mobile force: approximately `5.558e29 kJ/mol/nm`.
- The closest non-excluded pair is temporary B/VAL HG11 against generated
  C/MET24 SD, at **0.02128 Å**. The original heavy-atom overlap alone therefore
  understated the singularity after temporary hydrogen placement.
- Another generated/fixed heavy pair is separated by **0.14812 Å**.
- A fixed/fixed temporary-hydrogen pair is separated by **0.22722 Å**; hydrogens
  outside generated residues are currently massless despite being rebuilt for
  this temporary system, not retained scientific input coordinates.
- The production v7 refined PDB changes just one atom by more than 0.001 Å:
  C/MET24 SD moves approximately 10 Å; ordinary minimization does not amount to
  meaningful whole-generated-region relaxation in this case.

A separate staged soft-core initialization experiment is being evaluated before
changing production refinement. It retains bonded and chirality forces and must
return to the original full nonbonded force before independent validation. No
validation threshold or source coordinate is changed by the experiment.

Refinement v8 now freezes only identities of observed, non-generated heavy
atoms. Temporary rebuilt hydrogens and missing-sidechain completions remain
mobile; publication still rewrites only the requested generated coordinates.
An identity-level regression covers case-sensitive chains and insertion codes,
and the existing source-record preservation and seeded-repeatability tests pass
(28 focused refinement/adapter tests). This mobility correction alone is not
claimed to resolve the raw nonbonded singularities. Staged-initialization controls
are being compared with the former and corrected temporary mobility policies.

### First staged control: former temporary mobility policy

The separate `softcore-refinement-experiment/` control completed five soft-core
stages (`alpha=0.5,0.2,0.05,0.01,0.001`, 500 iterations each) followed by 1,000
iterations using the original nonbonded force. Generated/source-heavy mobility
matches the captured v7 system; this control does not include the v8 temporary
mobility correction. No D centres were detected at stage ends, and the fixed
particle displacement was zero. Every published junction is within the unchanged
1.20–1.45 Å gate (range 1.2950–1.3875 Å).

However, independent inspection of the rewritten PDB found **12 bond-length
errors, 17 angle errors, 3 Ramachandran outliers and 347 severe overlaps involving
generated residues**. This control is rejected. It demonstrates escape from the
nonbonded singularity, not geometry-preserving relaxation or usable output.
The script's final `getState()` energy includes the dormant soft-core force group
as well as the original force; it is not a standalone full-FF energy measurement.
The final minimization itself uses only original nonbonded plus bonded/chirality
forces through the integrator group mask. Production refinement has not adopted
this staged scheme.

### Second staged control: mobile temporary completions

Repeating the soft-core schedule with temporary hydrogens and missing-heavy
completions mobile does not qualify it. Fixed observed-heavy displacement remains
zero, but the final original-FF minimization produces **D-Cα at C/PRO7**. The
rewritten candidate also has 10 bond-length errors, 9 angle errors, 4 Ramachandran
outliers and 366 severe overlaps. It is rejected, irrespective of the stage
energies or successful junction distances. All stage-end chirality checks had
been clean until that final full-FF step; production's all-L admission gate must
therefore remain unconditional after the final minimization.

The next private control replaces soft-core electrostatics with bounded
repulsion-only initialization, followed by the original nonbonded force. An
attempt to add force-field-equilibrium X–H constraints failed at Context creation:
OpenMM forbids constraints involving massless scaffold particles. The replacement
control therefore retains the original equilibrium harmonic bonds without adding
constraints; exact scaffold immobilization is not traded for constrained dynamics.
Neither control is a production change or an acceptance claim.

The full fast suite after the v8 mobility correction passes: **1,408 tests**, two
optional-dependency skips and 31 slow deselections. The first direct-Python suite
launch omitted the environment's `bin` from `PATH`, causing six subprocess CLI
lookup failures; repeating under `micromamba run -n dvbfixer` resolves all six.

### Repulsion-only control: also rejected

The unconstrained replacement completed bounded-repulsion stages at
`k=1000,5000,10000 kJ/mol/nm²` (200 iterations each) and then 1,000 iterations of
the original full force field. Stage-end chirality remained L during repulsion,
but the full-FF step produced **D-Cα at A/LEU14** (also confirmed after PDB
serialization). Fixed observed-heavy displacement was zero. Independent output
inspection found **17 bond-length errors, 10 angle errors, 4 Ramachandran outliers
and 459 severe overlaps**; junctions alone were in range (1.3054–1.3636 Å).
The original full-FF final energy was approximately `158451 kJ/mol`, excluding
the dormant repulsion force through the explicit group mask.

This experiment is rejected as well. Neither the soft-core nor repulsion-only
Cartesian initialization tested here preserves enough local geometry to provide
a usable reconstruction. None has been adopted in production, and none is a
reason to weaken all-L, bond, angle or clash gates. The demonstrated software
corrections (runtime selection, independent gap masks, alignment and temporary
atom mobility) must not be described as solving the remaining reconstruction
quality problem. A geometry-preserving proposal needs a new controlled test,
not further blind energy-only minimization or publication of these candidates.

## Opt-in anchor-segment union experiment

At the user's request, a new explicit frozen-request option,
`cc89_context_geometry=anchor-shell-segment-union-v1`, adds a finite CA–CA
segment neighborhood to the existing source-only context criteria. A residue
qualifies if any observed heavy atom is within 10 Å of the segment between the
two observed anchor Cα atoms **or** meets the existing anchor-shell criterion.
Sequence flanks, padding, immutable-source rules and the 512-position budget
remain in force. Projection is clamped to the segment endpoints, not extended
to an infinite line. This experiment requires internal gaps with two anchors;
terminal directions are not invented. Default profiles retain their original
policy unless the request declares the option explicitly.

The control is the seed-7 no-terminal run: four gaps / 27 residues, no geometric
errors or D centres, but seven remaining severe overlaps. Some conflicts were
already present against PHE77/PRO79 absent from that gap's motif; new SER85
conflicts appeared during minimization, despite its correct input representation.
The latter cannot be attributed to missing motif atoms. Model-parser inspection
confirmed all 446 declared fixed atoms and the fixed-position target sequence;
SER85 coordinates match to float32 precision. No input-mapping correction was
made on the basis of those checks.

New request: `test/errors_diffusion/anchor-segment-union-input/request.json`.
The normalized input, sequences, generated mask, fixed atoms, seed and refinement
are unchanged from the control. Only the A108–120 context expands: 74→85 engine
positions, adding 84 observed heavy atoms, including PHE77 and PRO79. All previous
conditioning atoms are retained. The other three contexts remain identical.
Focused context/profile tests: 61 pass; the complete fast suite passes 1,411 tests
with two optional-dependency skips and 31 slow deselections.

The completed seed-7 experiment at
`test/errors_diffusion/anchor-segment-union-seed7/` **passes all existing independent
output gates** for the four internal gaps / 27 generated residues:

| Metric | Original context | Anchor-shell + segment union |
|---|---:|---:|
| Severe overlaps | 7 | 0 |
| Detectable D-Cα | 0 | 0 |
| Bond-length / bond-angle errors | 0 / 0 | 0 / 0 |
| Non-planar amides | 0 | 0 |
| Ramachandran / sidechain χ1χ2 outliers | 0 / 0 | 0 / 0 |
| Maximum fixed-heavy displacement | 0 Å | 0 Å |
| Junction C–N range | 1.3121–1.3606 Å | 1.3135–1.3613 Å |

No threshold was relaxed, and refinement remains v8. This is a successful
single-input, single-seed internal-gap reconstruction, not proof of experimental
conformation, general cohort qualification, terminal reconstruction or justification
to change the default context policy. Long terminal regions remain unmodeled.

### Integration after the successful control

At the user's subsequent request, fresh corrected Protpardelle CLI and Homology
mosaic requests now explicitly declare the segment-union option by default.
Preparation and validation still interpret older frozen requests without the
option using their original shell-only policy; saved request bytes are not
silently migrated. For a terminal per-gap view, the new versioned policy retains
the original shell selection, because a second observed anchor does not exist.
No terminal direction, additional missing residue or hidden conditioning atom is
invented. Engine budgets, source preservation and all output gates are unchanged.

## Terminal-window 20 Å union control: rejected

The user requested a full seed-7 rerun adding, rather than replacing, context for
terminal gaps. The separate option
`cc89_terminal_context=eight-observed-residues-20a-union-v1` requires the existing
anchor-segment union option. It additionally selects residues with an observed
heavy atom within 20 Å of any observed heavy atom in the first/last eight observed
residues, including the terminal gap's anchor. Short observed chains use only
available residues, never invented positions. Existing shells, flanks, padding,
source-only rules and output gates remain unchanged; internal-gap contexts are
identical to the current default. At the time of this control the terminal option
was opt-in, not a default-policy change.

The frozen input is `test/errors_diffusion/terminal20-union-input/request.json`;
the current-rules control is saved alongside as `baseline-current-rules.json`.
All ten gaps / 124 generated residues, source coordinates, sequence placements
and seed are unchanged. Engine-axis sizes, including generated residues:

| Terminal gap | Current default | With terminal-window union |
|---|---:|---:|
| A N-terminal, 20 residues | 50 | 167 |
| A C-terminal, 3 residues | 29 | 267 |
| B N-terminal, 31 residues | 70 | 310 |
| B C-terminal, 1 residue | 45 | 184 |
| C N-terminal, 31 residues | 70 | 309 |
| C C-terminal, 11 residues | 64 | 180 |

Every old conditioning atom is retained, and all views fit the unchanged
512-position limit. The four internal contexts remain 77, 85, 67 and 51 positions.
Seventy focused context/profile tests and Ruff pass. The completed full generation
run at `test/errors_diffusion/terminal20-union-seed7/` is **rejected**. All 124
generated residues/heavy-atom sets are present, and observed-heavy displacement
remains exactly zero. Private evidence is retained at
`test/errors_diffusion/.dvbfixer-diffusion-run.ud15e5oj/`.

| Whole generated/junction region | Raw merged candidate | Refined candidate |
|---|---:|---:|
| Bond-length errors | 38 | 31 |
| Bond-angle errors | 75 | 35 |
| Non-planar amides | 15 | 28 |
| Ramachandran outliers | 6 | 4 |
| Severe overlaps | 125 | 723 |
| Detectable D-Cα | 1 | 0 |

Final junction C–N distances span 1.3013–1.5382 Å, and the longest generated-region
C–N distance is 3.0430 Å. The unchanged connectivity and geometry gates therefore
also fail. Refinement removes detectable D-Cα but does not produce a usable
minimum; fewer bond/angle errors or finite energy are not acceptance.

Most raw bond/angle defects arise before refinement in the 31-residue N-terminal
tails: B has 14 bond-length / 37 angle errors, C has 23 / 35. The raw per-gap
native-conditioning fits also show substantial motif deformation before source
restoration: B N-terminal motif RMSD 2.4376 Å (maximum 8.4071 Å), C N-terminal
RMSD 2.1667 Å (maximum 10.4750 Å). These diagnostics distinguish neural motif
retention from publication's exact source-coordinate restoration. They warrant
further conditioning/representation investigation, not a conclusion that all
remaining failures are caused solely by steric omission or minimization.

This is one seed and one input, not a matched full-complex current-default
rerun/cohort. This failed experiment alone was not used to change a default or
validation threshold.

## User-supplied FASTA control: all requested gaps pass

The user subsequently supplied `test/errors_diffusion/selected.fasta` as the
authoritative target sequence for the same normalized structure. With the same
seed 7, corrected profile, terminal-window union, refinement v8 and unchanged
validation gates, the run at `test/errors_diffusion/terminal20-fasta-seed7/`
**passes every independent output gate**.

The FASTA has chain lengths A/B/C = 134/213/226, instead of SEQRES 147/245/263.
It requests seven gaps / 42 residues, not ten / 124. In particular, the two
31-residue N-terminal B/C regions are absent from this target; A's N-terminal
region is seven residues, and C's C-terminal region is five. This is not silent
terminal cropping: no `--no-terminal` was used, and every missing position in the
user's FASTA is requested. The source PDB and FASTA bytes were preserved. Frozen
inputs, target sequences and SHA-256 provenance are retained in
`test/errors_diffusion/terminal20-fasta-seed7-input/`.

Final independent metrics: zero severe overlaps, zero D-Cα, zero bond-length and
bond-angle errors, zero non-planar amides, zero Ramachandran and χ1χ2 outliers,
zero missing/wrong generated atoms or residues, and zero observed-heavy
displacement. Junction C–N distances span 1.3125–1.3630 Å.

This qualifies the reconstruction against the implemented gates for this one
input/sequence/seed. It neither rescues the longer SEQRES reconstruction nor
establishes the experimental conformation of the newly modeled regions or a
general terminal-backend acceptance cohort.

### Terminal-union default integration

Following the user's approval after the FASTA control, fresh corrected
Protpardelle CLI and Homology mosaic requests now record both the anchor-segment
and terminal-window union options by default. Frozen requests without the terminal
option retain their former terminal shell selection; no saved requests are
migrated. Internal-gap selection, the 512-position budget, refinement v8 and
scientific output gates are unchanged. Terminal/long-gap reconstruction remains
experimental and can still fail, as the longer SEQRES control demonstrated.

The supplied FASTA is now tracked byte-for-byte as
`tests/fixtures/protein_context/5j13_selected.fasta`, SHA-256
`853b97081f53b6dee15082e8d169c2728cbfc9387afbb291cade4872b082044d`.
Regression tests verify default option declaration, exact source-bound context
sizes and frozen-request replay without running an external diffusion engine.

## cc94 MPS internal-gap seed panel (2026-10-27)

Changing the seed is a practical recovery strategy when a diffusion candidate
fails geometry after refinement. This case demonstrates seed sensitivity, not
guaranteed recovery or general cc94/MPS qualification. Do not transfer the earlier
cc89 successes above to the new checkpoint or silently replace failed trials.

### Fixed protocol

The native Apple M3 Pro run used repository commit `5bdeff1`, profile
`protpardelle-1c-cc94-mps-v1`, pinned cc94_epoch3100, Torch 2.6.0 and the
[cc94 ODE/native-relative policy](../../deploy/protpardelle-1c/cc94.md): 500 updates,
scale 1, churn 0, no per-step reinjection and no MPS-to-CPU sampling fallback.
Refinement used separate Python 3.12 / OpenMM 8.5.2 CPU / NumPy 2.4.6,
one CPU thread, with default geometry-gated anchor relaxation and B weighting.
No standalone full-system `minimize` or solvent stage was run.

The source was the tracked `5j13_selected.pdb` (SHA-256
`63154130ff3de95ab1a135154370ac91e07d8b776bb16067096c0c62c91e3106`).
Water and sulfate were removed in a separate protein-only copy, SHA-256
`34e144ddb55f39b83f79b1aedeb44c075c430cac6e8d6b4ba3ae23c7f6874ba9`;
the original file was unchanged. Targets came from complete input SEQRES with
explicit `--no-terminal`, yielding A/B/C lengths **124/213/221** and **four
internal gaps / 27 generated residues**. This is not the longer ten-gap experiment
or the seven-gap user-FASTA control above.

After seed 7 failed, the additional panel **19, 42, 101, 202** was fixed before
execution. Every seed was run in its own bundle with otherwise unchanged settings;
none was replaced and execution did not stop at the first success.

| Seed | Outcome | Refinement finding |
|---|---|---|
| 7 (initial comparison) | Rejected; raw retained | Stage 1 leaves one generated A115 Ramachandran outlier (φ=32.4°, ψ=−89.4°) |
| 19 | **Published success** | Generated-only stage 1 passes all independent gates |
| 42 | Rejected; raw retained | Stage 1 leaves four severe overlaps involving generated C/SER142 |
| 101 | **Published success** | Generated-only stage 1 passes all independent gates |
| 202 | **Execution failure; no published bundle** | Workspace-size check aborts the wrapper; retained stage-1 result separately passes independent validation |

Seeds **19 and 101** restore every requested generated residue/heavy atom, retain
zero observed-heavy displacement and have zero bond-length/angle errors,
nonplanar amides, Ramachandran/χ1χ2 outliers, severe overlaps and detectable D-Cα.
Neither needs stage 2. For seeds 7 and 42, residual findings do not identify an
eligible problematic observed anchor, so the controller does not run stage 2.

The seed-202 execution failure is **not** a geometric rejection or a published
success. Its retained workspace measured about 20 MB after the run, below the
2 GB limit; the checker also reports filesystem traversal errors as limit
violations. The precise trigger was not captured, so actual transient quota
exhaustion versus a filesystem race is unresolved. Post-hoc validation of its
saved result passes all independent gates without further minimization, but does
not retroactively repair the failed execution/publication path. Thus the four
additional attempts produced **two published successes, one geometric rejection
and one infrastructure failure**, not four comparable geometric trials.

Local, untracked evidence is retained below `.artifacts/apple-cc94/`:
`5j13-seed-panel-protocol.json`, `5j13-seed-panel-results.json`, per-seed command
records/logs, the `runs/5j13-internal-seed{19,42,101}/` bundles and
`5j13-seed202-posthoc-validation.json`. Seed 7 is in `runs/5j13-internal-seed7/`.
Private workspaces retain stage reports and coordinate attempts, including the
unpublished seed-202 result. These are operator-local evidence, not CI fixtures.

Use a declared small panel, retain every outcome and inspect per-candidate gates;
do not select on finite energy, CLI exit code or raw appearance alone. A rejected
staged-refinement bundle can contain the **raw** PDB, not its improved but still
failing stage-1 output. See the
[multi-seed CLI guidance](../commands/model.md#retrying-geometric-rejections-with-multiple-seeds).
This case supports an operational recovery practice, not native-loop accuracy,
seed-independent reproducibility or a terminal-gap acceptance claim.

### Subsequent rejected-output retention correction

The seed-panel outcomes above describe the original executions; their bundles
are not rewritten. The shared staged refiner now preserves the last completed
processed candidate even when geometry is rejected: stage 2 if it returned
coordinates, otherwise stage 1. Raw coordinates, stage-1 evidence and the attempt
report remain separate digest-bound bundle artifacts. Raw is retained as the
main rejected PDB only when processing produced no coordinates. If a later stage
raises an expected refinement error, its error remains visible alongside the
earlier retained result. Scientific acceptance, source-motion/chirality limits
and downstream refusal of rejected bundles are unchanged. This is an output
retention correction, not a rescue or relabeling of any historical failed seed.

Native CPU reprocessing of the saved seed-7 raw candidates, without new neural
sampling, verified publication and independent resume validation of both rejected
bundles. 5J13 now publishes stage 1 with one Ramachandran outlier and zero severe
overlaps, instead of its raw candidate with 18 overlaps. The 9CQD internal-gap
control publishes stage 2 with zero severe overlaps, one nonplanar amide and one
Ramachandran outlier; that execution also exceeds the W29 backbone movement limit
(0.5075 Å versus 0.5 Å). Both remain `validation_failed`, with raw/stage-1/report
artifacts preserved separately. These are new CPU minimization executions, not
bit-identical replacements for the earlier stage coordinates. Local verification:
`.artifacts/apple-cc94/processed-rejection-verification.json` and
`runs/{5j13,9cqd}-internal-seed7-processed-rejection/` below that evidence directory.
