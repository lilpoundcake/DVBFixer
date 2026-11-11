# 1HE8: motif portfolios and joint generation experiments

Date: 2026-11-10. **Case-specific experimental evidence and hypotheses, not
scientific qualification or a production/default decision.** This record
separates executed computations, user visual review and proposed explanations.
P3 was not implemented. Historical requests, failed attempts and cohort outcomes
remain unchanged. MODELLER remains the default backend.

## Scope and common conditions

The 1HE8 target has chain A965 and B166. The source contains 749 observed A
residues and 166 observed B residues. Ten gaps in A contain 216 generated
residues / 1712 heavy atoms. These are source/target identities, not renumbered
engine proxy chains. The A target starts at author residue A144.

The original motif-based runs used seeds **7, 101, 102, 103, 104**. New sampling
experiments below used the pinned cc94 checkpoint on MPS, 500 ODE updates per
call, scale 1, churn 0, native relative position/chain features and no CPU
fallback. No checkpoint/config changes or cc89 rotary intervention were made.
CPU refinement is a separate process and does not count as diffusion sampling.
Refinement controls used no solvent.

The public audited token budget remains **512**. Private execution was explicitly
authorized beyond it: initially a 1024 ceiling, subsequently 1131 for complete
1HE8 input. Execution without OOM is not evidence of scientific validity beyond
512. Retained guards were process RSS 8 GiB / MPS driver allocation 6 GiB.

Observed-motion limits were unchanged: per-residue backbone RMSD 0.5 Å,
Cα displacement 0.5 Å, maximum heavy-atom displacement 1.25 Å, relative to the
original source, not re-anchored between optimizer restarts.

## Executed experiments

Counts below are diagnostic findings, not independent defects or accuracy
measurements. Aggregate gate categories must not be added to their detailed
counterparts as if they were additional physical defects. Process completion,
materialization integrity and scientific acceptance are different outcomes.

| ID | Operation actually executed | Observations / status |
|---|---|---|
| E1 | Original independent local-motif generation: ten gaps × five seeds | 50 actual archived motif inputs were later audited against frozen source/request. These are the donors underlying the portfolio, not joint samples. |
| E2 | Bounded multi-parent portfolio search over generated-only donors | 64 evaluations, beam width 4; whole gaps selected from single donors in a common observed frame. The search exhausted its budget; not an exhaustive optimum. The assembled 1HE8 still had defects. |
| E3 | Issue-local post-portfolio repair and a separate bounded-anchor follow-up | Earlier repairs could freeze other generated gaps. The bounded follow-up for 1HE8 remained rejected, including Ramachandran A976 and observed-anchor A895 movement. This was not new sampling. |
| E4 | Pair-joint pilot: A436–456 and A489–497, seed 7, one call, 155 tokens, 201 generated atoms | Independent native reconstruction and generated-only refinement completed. Pair geometry remained rejected around A437; no topology success was established. |
| E5 | All ten gaps jointly generated from the union of local motifs, seed 7, 736 tokens; all-generated + spatial-shell refinement | One call / 500 updates. Final: 40 detailed clash findings, Ramachandran A326/A968, 12 movement violations among 235 observed shell residues. Rejection independently reproduced. |
| E6 | Complete observed-complex conditioning: all ten gaps, seed 7, 1131 tokens; spatial-shell refinement | One call / 500 updates. Final: 17 clash findings, non-planar A343–344, Ramachandran A543, 19 movement violations among 276 observed shell residues. Rejection independently reproduced. |
| E7 | Retained E2 patchwork, all 1712 generated atoms plus all direct anchors minimized together | No new sampling. Nineteen anchor residues / 148 atoms mobile under restraints; other observations fixed. All 19 initial clash findings disappeared; one detailed Ramachandran finding and three anchor-motion violations remained. Rejected. |
| E8 | Source-neighbor components sampled on local motifs, seed 7: nine gaps jointly at 704 tokens, singleton at 44 tokens; assembly then all-generated/direct-anchor refinement | Two calls, 500 updates each. Final: 23 clash findings, two generated-backbone-break findings, two bond-length and two bond-angle findings, one non-planar amide, one Ramachandran finding, five anchor-motion violations. Rejected. |
| E9 | New all-ten-gap local-union multi-seed run: seeds 7, 101–104, 736 tokens per call | All five calls and five independent native-to-PDB exports completed. **No minimization** in this experiment. Outputs are raw proposals, not five accepted structures. Seed 7 was newly executed, not silently reused from E5. |

E5/E6 spatial mobility included every generated atom, observed residues within
5 Å of the initial proposal, direct anchors and one observed sequence-neighbor
layer. E7/E8 used every generated atom plus **direct anchors only**, not that
larger shell. These refinement arms therefore are not matched comparisons.
Their controls used CPU, 500 full-force-field iterations and four configured
restarts, with retained intermediate/search diagnostics.

The reported refined outputs of E5–E8 passed global L-chirality and exact record
preservation outside their declared mobile masks. Those checks do **not** imply
geometry acceptance, correct fold, secondary structure or absence of threading.

### Retained execution failures and preparation controls

- The first E5 attempt failed pinned preflight **before sampling** because an
  earlier pilot had created 17 engine bytecode files. They were retained with
  checksums, not deleted; engine source, weights and preflight were not weakened.
  This attempt had zero model calls/updates and is not a scientific seed result.
- E6 complete-input preparation initially reported 1131 required tokens against
  the then-current 1024 execution ceiling. Prepare-only and actual execution
  refusal before weights/runtime preflight were retained. Only after separate
  authorization was the private ceiling raised to 1131 and E6 executed.
- Test runs are engineering evidence, not scientific replicates. A final earlier
  full fast run passed 2371 tests, with three optional Torch/Biotite skips and
  33 slow tests deselected. Its preceding failed run had six missing-CLI-PATH
  errors and an older loaded version of a subsequently corrected test. Later
  component-focused verification passed 26 tests; do not describe the earlier
  full-suite result as verification of every later script/change.

## Input audit: what was and was not wrong

Independent checks confirmed deposited SEQRES A965/B166, observed residue
placement, exact source coordinates (0 Å difference), packed sequence/physical
axes, fixed/generated masks and motif reader/placement/crop replay. All 50
archived independent motifs matched their expected inputs and archived digests.
No numbering, atom-order or coordinate-frame error was established by that audit.
Replay was before the common random rigid rotation, without model weights or
sampling; it was **not** a captured historical denoiser tensor.

However, E5's 736 tokens were **not the full observed complex**:

| Conditioning coverage | E5 local union | E6 complete complex |
|---|---:|---:|
| Observed A residues | 520 | 749 |
| Observed B residues | 0 | 166 |
| Fixed conditioning atoms | 4218 | 7395 |
| Omitted observed residues | 395 | 0 |
| Generated residues | 216 | 216 |
| Total tokens | 736 | 1131 |

A241/A242/A284 were absent from the local motifs and their union; new A589/A590
collided with this omitted environment. A150 **was** conditioned, yet A150/A323
also clashed. Thus omitted obstacles are a plausible contributor, not a complete
explanation. E6 restored all observations, including OXT B166, but still failed.
OXT was excluded only from alignment against dummy predicted-native coordinates,
not from conditioning or the preserved source output.

E8/E9 also use local motifs. Their poor appearance cannot be attributed simply
to feeding the entire observed complex: E8's main component is 704 tokens and
E9 is the 736-token local union. “All gaps generated together” and “all observed
structure supplied as conditioning” must not be conflated.

## How the patchwork was selected

The E2 selector did **not** score secondary-structure content or visual realism.
It used full assembled-complex assessments, prioritizing admission failures,
then the number of actionable geometric findings. Per-gap issue counts ordered
the bounded beam exploration; they did not replace whole-assembly assessment.
One donor supplied each complete gap; original observed coordinates remained
fixed. The selected seeds for gap indices 0–9 were:

`101, 103, 103, 103, 104, 103, 7, 103, 103, 103`.

Local geometrical cleanliness is not evidence of compatible relative placement,
absence of interwoven loops or correctness against a reference structure.
The earlier description “selected structured fragments” overstated the actual
algorithm: visually structured conformations were present in donor outputs,
but secondary structure was not an explicit selection criterion.

The ordinary whole-candidate refiner already shared the generated mask across
all gaps. The **post-portfolio issue-local repair** was different: it selected
problematic gaps and froze others. E7 explicitly tested the user's requested
post-assembly all-generated/direct-anchor joint minimization. This new private
control is not an already integrated mandatory public pipeline stage.

## User visual review (qualitative, not validation)

The user found original local-motif donors/patchwork more protein-like, with
recognizable secondary structure, than the large joint proposals, described as
“spaghetti.” The user also reported apparent interwoven loops, including
A436–456/A489–497. These are qualitative observations, not detector-established
topological classifications or quantitative secondary-structure measurements.

Ten visualization-only ensemble PDBs were exported from the **generated-only
donors actually used by E2**, each with five models and its observed local motif.
No further fitting/refinement or coordinate changes were applied. Model order:
**1=seed 7, 2=101, 3=102, 4=103, 5=104**. Gemmi readback checked all 50 models.
The user's acceptable-looking model sets, in file order, were:

| Display number / gap index | A residues | Visually acceptable models |
|---|---|---|
| 1 / 0 | 322–352 | 1–5 |
| 2 / 1 | 375–378 | 1–3 |
| 3 / 2 | 436–456 | 1–5 |
| 4 / 3 | 489–497 | 5 |
| 5 / 4 | 528–619 | 1–5 |
| 6 / 5 | 754–760 | 1–5 |
| 7 / 6 | 896–901 | 1–5 |
| 8 / 7 | 968–980 | 2–5 |
| 9 / 8 | 1038–1046 | 1 |
| 10 / 9 | 1085–1108 | 2–3 |

This is **36 visually approved fragment variants / 50**, not 36 scientifically
accepted structures. The previous automatic donor falls outside these choices
for display gaps 2, 4, 9, 10. Joint compatibility of the approved fragments was
not established. Their 75,000 possible combinations were **not** enumerated or
minimized. Unlisted models are not automatically scientific rejections.

## Working hypotheses, not conclusions

| Hypothesis / possible explanation | Supporting observation | Unresolved alternative / needed control |
|---|---|---|
| Smaller local tasks may suit this checkpoint better than large joint tasks | Original local donors were visually preferred; 704/736/1131-token proposals remained problematic | Task size, simultaneous gap count, conditioning coverage and sample selection all changed. Need matched multi-seed/context controls; no general cutoff established. |
| Portfolio selection contributes to apparent patchwork quality | Five donor seeds were available; visually useful alternatives exist | Most early joint comparisons used one seed, versus selected pieces from five. E9 supplies five raw joint seeds, but not matched refinement/selection or a quantitative comparison. |
| Omitted observed obstacles contribute to clashes | Local-union A589/A590 met previously unconditioned A241/A242/A284 | Conditioned A150 also clashed; complete-context E6 still failed. Coverage alone is insufficient. |
| Transitive source-flank grouping can make overly large tasks | E8 groups nine of ten gaps into a 704-token problem | Eight-residue flanks and ≤6 Å heavy contacts are a heuristic, not a learned/generated-path interaction measure. Smaller alternative groupings were not qualified. |
| Independently good fragments can be mutually incompatible | User sees interwoven loops after assembly/minimization | No implemented threading detector or independently established topology labels. Source-neighbor proximity alone cannot diagnose it. |
| Joint minimization can improve local geometry without resolving the fold | E7 removed all 19 reported initial clashes but remained rejected and visually problematic | No validated chain-passage/threading check; no proof that minimization can or cannot repair every such case. |
| Explicit structural ranking and user-guided donor choice could help | Automatic choices disagree with visual choices at four gaps | No DSSP/secondary-structure ranking benchmark, GUI selector or validated compatibility optimizer was implemented. |

These observations do not establish weights-only causality, a defective model
family, motif superiority in general, correct experimental structure recovery,
or that every long joint request fails. Fewer clashes and all-L geometry do not
prove protein-like global organization. No held-out reference accuracy benchmark
or statistically controlled comparison was completed here.

## Proposed next work, not executed here

- Preserve local-motif multi-seed donors and explicit lineage. Compare automatic
  and user-filtered selections with separate geometry, secondary-structure and
  cross-fragment compatibility diagnostics rather than an unexplained score.
- Investigate better spatial grouping without silently splitting an interacting
  component to fit a token limit or calling independent predictions joint.
  Current E8 rule: up to eight observed flank residues on each side, ≤6 Å heavy
  contact, transitive connected components. No guarantee about future loop paths.
- Make all-generated/direct-anchor minimization an explicit versioned
  post-assembly operation; retain pre/post states and unchanged movement limits.
- Implement/qualify the separately planned
  [intersection/threading detector](../plans/diffusion-geometry-repair.md#generated-loop-intersectionthreading-detector-requested-follow-up-2026-11-10):
  distinguish 3D crossings, open-path threading suspicion, legitimate topology
  and ambiguous closure; check raw, assembled and refined structures.
- Consider a GUI portfolio editor for seed switching, overlays, user choices,
  whole-assembly compatibility preview and reproducible Python materialization.
  This is a proposal, not an approved implementation plan or shipped capability.
  Human preference must not bypass independent acceptance.

## Evidence locations and reproducibility boundary

Paths below are relative to the repository root. They are **retained local
artifacts under `tmp/`, not guaranteed to transfer through Git**. Do not claim a
fresh checkout includes the structures/runtimes. Preserve/archive them separately
before cleanup or transfer. This tracked narrative records findings; reports,
native arrays and traces remain authoritative for detailed custody/numerics.

| Evidence | Local path |
|---|---|
| Original cohort, source and managed diagnostic archive | `tmp/real-complex-cc94-2026-11-06/` |
| E2 selection, five hashed donors, 64 trials | `tmp/composed-input-preparation-worktree/tmp/portfolio-publication-v1/retained-panel-01/1HE8/selection/` |
| E3 bounded-anchor follow-up | `tmp/composed-input-preparation-worktree/tmp/portfolio-anchors-v1/retained-controls-01/1HE8-repair/` |
| E4 pair pilot | `tmp/joint-gap-pilot-1HE8-v1/` |
| E5 failed preflight and actual execution | `tmp/1HE8-full-joint-spatial-v2/attempt-01/`, `attempt-02/` |
| E5 input/archived-motif/independent audits | `tmp/1HE8-full-joint-spatial-v2/{input-audit-01,archived-independent-input-audit-01,independent-audit-01}/` |
| E6 preparation/refusal controls and executed run | `tmp/1HE8-complete-context-v3/` and its `attempt-01/` |
| E7 all-generated/direct-anchor patchwork control | `tmp/1HE8-patchwork-all-generated-anchors-v1/attempt-01/` |
| E8 components, native assembly proof, refinement | `tmp/1HE8-neighbor-motif-v4/attempt-01/` |
| E9 all-five-seed raw proposals and exports | `tmp/1HE8-all-gaps-multiseed-local-v2/attempt-01/` |
| Visual ensembles and exact user selection | `tmp/1HE8-motif-seed-ensembles/manifest.json`, `user-visual-selection.json` |

The broader original ten-case cohort's **27/50** outcomes remain historical and
unchanged; the corrected 1N8Z input-preparation arm stays separately identified.
None of the post-hoc 1HE8 trials replaces an original outcome.

Related context: [geometry repair plan](../plans/diffusion-geometry-repair.md),
[portfolio publication and calibration](portfolio-publication-and-anchor-symmetry.md).
