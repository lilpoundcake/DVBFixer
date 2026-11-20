# Experimental structure coordinate completion: current handoff

Updated: 2026-11-19. This is the current continuation entry point for the retained
Mac reconstruction work. It supersedes conflicting next-action instructions in
the older handoffs, not their evidence, frozen requests or recorded failures.
This documentation update does not launch or authorize additional experiments.

**Current priority — latest user decisions:** investigate actual **diffusion-based
gap restoration**, not MODELLER proposals or repeated attempts to minimize them.
Use explicitly prepared **protein-only** inputs for both sampling and subsequent
refinement; heterogen support is deferred. For 1NON, the question is restoration
in the intact four-chain protein scaffold, not a silently cropped local context.
Preserve native observations before any source-coordinate copying. Respect the
benign restoration scope and explicit exclusions in section 2. These
instructions supersede conflicting historical next steps below. This update
does not itself launch a new checkpoint run or approve an unbounded cohort.

## 1NON: which model saw which input?

**Practical objective (user clarification):** replace the MODELLER gap-building
step with a learned proposal that yields a plausible, preparation-ready peptide
of the supplied sequence in the observed protein environment. Evaluate the
output the user actually receives: both junctions, backbone geometry and path,
side-chain chemistry, contacts with the full assembly, preservation of observed
atoms, and passage through downstream preparation. A known-loop RMSD is a
secondary benchmark, not a per-gap acceptance requirement. A plausible model
may differ from the unobserved conformation. For the next prototype, test the
learned backbone separately from side-chain construction: close both ends with
the observed scaffold fixed, pack only the new residues, and run existing checks.
Record raw versus rebuilt atom ownership and all local geometry failures.

| Model on MPS | Input actually seen by the model | Result for 1NON |
|---|---|---|
| Protpardelle-1c (`cc94_epoch3100.pth`), full tetramer | Four 179-residue chains; all 5,040 observed atoms in one frame; residues 73–90 missing in each chain. `motif.pdb` contained the **entire observed scaffold**. | Seed 23, 500 updates. Fixed-scaffold inpainting preserved the observations. Both raw outputs (baseline and placement correction) failed loop geometry. |
| Same cc94, local motifs around the gaps | Four separate MPS calls, seed 23, 500 updates each; A/B/C/D: 114/107/98/83 positions and 745/694/618/496 fixed atoms. | All calls completed. Combined raw output in the tetramer fails geometry: 590 severe contacts vs 357 in the full-input raw baseline; 2 vs 4 failing junctions. |
| cc94, full own chain per loop | Same four seed-23 calls; all 1,260 observed atoms of that monomer plus the identical partner atoms as in the local control. | Own-monomer clashes 528 → 171 (A/B/C/D: 80/74/6/11); total severe contacts 590 → 235. One failing junction (C90–91 at 1.031 Å). Raw geometry rejected. [Receipt](../../tmp/1non-own-chain/attempt01/receipt.json); [inspection PDB](../../tmp/1NON-inspection/07_raw_own_chain.pdb) and `07_input_motif_A` through `D` in the same directory. |
| Protenix | Historical: 354 residues in 14 fragments. New [native feature readback](../../tmp/1non-backends/attempt01/protenix-preflight/result.json): all 716 tokens, 5,596 atom positions, 72 missing residues masked from templates. | Full-tetramer seed 7 actually entered MPS model execution; [device guard](../../tmp/1non-backends/attempt01/protenix-guard-torch/device-limit-exceeded.json) stopped it before an endpoint. Initial launch/configuration failures retained separately. |
| RFD3 | Four chains, 716 positions, 5,036 fixed canonical atoms, 556 generated; four OXT source-only. The [pipeline readback](../../tmp/1non-backends/attempt01/rfd3-preflight/pipeline-readback.json) verifies all 712 peptide bonds including eight junctions. | Full MPS forward attempted using the 1KTZ tiled pair/Transition and cache hooks, seed 101. The first trace used stale 1KTZ atom counts; the corrected [declaration](../../tmp/1non-backends/attempt01/rfd3-mps-axis/declaration.json) binds the 1NON counts. Sinusoidal pair features completed; reference pair features were in progress when the [driver guard](../../tmp/1non-backends/attempt01/rfd3-mps-axis/execution/device-limit-exceeded.json) stopped the process at 6.13 GiB. No denoiser step or final structure. |
| Boltz-2 | The original builder refused four gaps. Updated shared-template connector maps each chain's sequence and native atom/token axis independently. All four chains reached the model on MPS: 5,036 observed + 556 generated canonical atoms; four true-terminal OXT explicitly source-only. | [Final consumer](../../tmp/1non-backends/attempt01/boltz-terminal-run/final-consumer-axis.json) verified the axis. [Guarded execution](../../tmp/1non-backends/attempt01/boltz-terminal-guard/execution.json) stopped on MPS driver memory, before any candidate PDB. |

**RFD3 memory follow-up:** a separately declared seed-101 replay accumulated the reference embedding in 64-row tiles directly into the motif tensor. The 25 CPU/MPS operator controls passed. The [full-size buffer receipt](../../tmp/1non-backends/attempt01/rfd3-mps-combined/combined.json) confirms the combined 5592×5592×16 result occupied the original buffer, at 4.22 GiB driver allocation after cache release. The [device guard](../../tmp/1non-backends/attempt01/rfd3-mps-combined/execution/device-limit-exceeded.json) subsequently recorded 8.97 GiB; no first denoiser step or PDB. The next native operations broadcast the single-atom features over the pair axis and expand token features to atom pairs; the trace does not distinguish which allocation caused that peak.

**RFD3 operation-level measurement:** the [synchronized boundary trace](../../tmp/1non-backends/attempt01/rfd3-mps-operation-probe/pair-operations.jsonl) identifies the first allocation crossing 6 GiB: the sum of the two projected single-atom features broadcasts to a new `[5592,5592,16]` tensor. Driver allocation rose from 4,533,731,328 to 6,546,997,248 bytes before adding that temporary to the pair tensor. The token-to-atom gather was never entered; its actual footprint is still unmeasured. No denoiser step or PDB. Next target: compute the inner sum by query-row tiles and add each tile directly into the owned pair buffer, preserving `pair + (left + right)` order.

**Full initializer now completes on MPS.** The seed-101 replay tiled the single-atom broadcast, token gather and pointwise pair MLP into the owned pair tensor; it tiled the uncontracted key axis of the first pooling contraction while preserving both token reductions. The [boundary trace](../../tmp/1non-backends/attempt01/rfd3-mps-finish/cache-boundaries.jsonl) records initializer return and entry into the diffusion encoder. The first noisy input was retained, but there is no completed denoising step or endpoint. A separate [local bias projection control](../../tmp/1non-backends/attempt01/rfd3-mps-gather/declaration.json) gathered native neighbour indices before the linear projection on eligible sparse calls. Both runs hit the same sampled 7.25 GiB driver peak; the last module boundary logged was entry into encoder block 0. Four CPU/MPS numerical checks compare pair updates, token pooling and local attention outputs to the native expressions. These checks are operator-level, not full-model equivalence. Finishing this on the Mac requires further denoiser-side memory engineering or a separately authorized resource policy; the existing 6 GiB threshold remains active.

**Sparse attention follow-up:** query-row tiling of K/V gather, logits, full-key softmax and gating passed CPU/MPS comparisons with the native gated attention (including repeated indices). In the [full-MPS replay](../../tmp/1non-backends/attempt01/rfd3-mps-sparse-rows/execution/execution.json) the atom encoder completed and the first denoiser reached `diffusion_token_encoder`. The [device receipt](../../tmp/1non-backends/attempt01/rfd3-mps-sparse-rows/execution/device-limit-exceeded.json) records 6.22 GiB at termination. This does not identify a specific allocation inside the token encoder or establish the eventual peak for 200 steps. No complete step or PDB. A higher guarded device budget or moving the job to a different host is a resource decision, not an implicit retry.

**Authorized 8 GiB RFD3 replay completed; geometry rejected.** The same seed-101 sparse-row arm returned a native endpoint in 702.9 s (199 diffusion-module calls, 200 configured timesteps), peak sampled MPS driver 6.30 GiB. Before source copying, 5,036 represented observations matched after one rigid fit to maximum 4.91e-6 Å. Four terminal OXT are source-only. Raw coordinates retain 46 loop–scaffold and two intraloop severe contacts, four failed junctions and inverted Cα chirality at B90. [Native readback](../../tmp/1non-backends/attempt01/rfd3-mps-sparse-rows-8g/native-assessment.json); [inspection PDB with original observed records](../../tmp/1NON-inspection/08_RFD3_source_records.pdb). This supersedes the earlier memory-blocked status for this arm.

**RFD3 seed reproducibility:** [prospective seeds 102/103](../../tmp/1non-backends/attempt01/prospective-102-103/comparison.json) both returned full endpoints with native fixed-atom maximum displacement below 0.000005 Å. Across 101/102/103, loop–scaffold contacts were 46/63/46, additional generated–generated contacts 2/6/0, and failed junctions 4/2/2. All three raw candidates fail geometry. [Native B90 readback](../../tmp/1non-backends/attempt01/prospective-102-103/B90-native-identity.json) verifies the same distinct ASN atom names on all three seeds. B90 Cα–Cβ is **0.0306/0.0314/1.5097 Å**: the negative signed chirality on 101/102 accompanies Cβ collapse, already present before alignment/serialization. Seed 102 also flags D90; seed 103 passes global L. The declared recurrence condition held, so no generated-only refinement was launched. The [8 GiB Protenix trial](../../tmp/1non-backends/attempt01/prospective-102-103/protenix-guard/execution.json) stopped on driver memory without an endpoint.

**RFD3 bond audit and experimental holdout:** all 568 native covalent bonds touching the 556 generated heavy atoms were measured per seed. Counts shorter than 0.9 Å were 13, 15 and 4 for seeds 101, 102 and 103. ASP A73 side-chain collapse recurred on all three seeds. Even globally L seed 103 has four A73 bonds between 0.012 and 0.031 Å; it does not qualify for the conditional refinement. [Bond inventory](../../tmp/1non-backends/attempt01/rfd3-bond-audit/summary.json).

One new seed-101 control additionally withheld experimentally observed residues 137–154 in all four chains. The original four gaps remained unknown. A single fit on retained scaffold atoms, excluding every withheld atom, gave maximum fixed error 4.25e-6 Å. Across the 72 withheld residues (580 heavy atoms), RMSD was 4.37 Å Cα, 4.36 Å backbone and 5.20 Å all heavy. Chain Cα RMSDs A/B/C/D were 5.94/3.27/5.02/2.33 Å. This estimates recovery on this specific held-out segment only. [Scores](../../tmp/1non-backends/attempt01/holdout-137-154/score.json); [masked input](../../tmp/1NON-inspection/09_holdout_input.pdb); [raw inspection output](../../tmp/1NON-inspection/09_RFD3_holdout_137-154.pdb).

**Packed seed-103 minimization and last-step probe:** minimizing the 556 generated heavy atoms with disposable H under ff19SB preserved all 5,040 source records. Severe clashes dropped from 54 to 7; final projected RMS force was 5.58 kJ/(mol·nm). A72:C–A73:N shortened to 1.1296 Å and the stereo shape guard failed. [Inspection PDB](../../tmp/1NON-inspection/11_RFD3_packed_minimized.pdb); [audit](../../tmp/1non-backends/attempt01/packed103-ff/independent-audit/refinement-receipt.json). Separately, on the archived last noisy state of the **raw** seed-103 run, the tiled MPS call reproduced the endpoint after the sampler update to max component error 1.91e-6 Å. The untiled MPS reference exceeded the 8 GiB guard. The native CPU denoiser using the same initializer outputs, coordinates, weights and terminal noise level agreed with the tiled MPS call on the generated atoms (RMS/max displacement 4.51e-5/2.57e-4 Å). Denoised A73 Cα–Cβ: 0.02774/0.02776 Å respectively. This tests one call given the saved state and shared conditioning, not divergence earlier in sampling. [Paired readback](../../tmp/1non-backends/attempt01/prospective-102-103/last-step-cpu-reference/probe.json).

**ASP73 timing:** a seed-103 [199-call diagnostic replay](../../tmp/1non-backends/attempt01/prospective-102-103/asp73-replay/transport/asp-trace.jsonl) reproduced the first noisy tensor exactly and the A73 collapse qualitatively, not the full endpoint bit-for-bit (maximum absolute coordinate difference 0.622 Å). In the denoised predictions A73 Cα–Cβ stays below 0.9 Å from call 119 and below 0.1 Å from call 160. B73 transiently collapses but recovers, as do the normal final C/D controls. All four [input feature slices](../../tmp/1non-backends/attempt01/prospective-102-103/asp73-replay/transport/asp-input.json) encode the same ASP type and real atom identities; reference positions are all zero and reference masks false. Only identifier/index fields differed among the recorded single-residue slices. This does not compare their spatial neighbourhoods. The auxiliary sequence-head argmax for A73 stays at index 7 from call 141 while the final B/C/D argmax is 3. This correlation does not establish feedback from the classification head to coordinates. Next investigate reference-feature construction for known sequence / unknown position, specifically whether the checkpoint expects residue-local chemical reference geometry on this path.

**ASP reference-feature hypothesis resolved for the pinned inference path:** the [live native transform receipt](../../tmp/1non-backends/attempt01/prospective-102-103/reference-policy/policy.json) records `generate_conformers=True`, `generate_conformers_for_non_protein_only=True`, `provide_reference_conformer_when_unmasked=False`. `CreateDesignReferenceFeatures` initializes zero reference coordinates/false masks and excludes protein atoms from conformer generation even when their sequence is known. The defaults in the pinned training dataset configuration also exclude proteins from conformer generation. The native transform emitted zero reference positions and false masks for all 5,036 observed and 556 generated protein atoms. Observed scaffold geometry is provided through separate motif features. Thus no ASP conformer was dropped by the connector; forcing one in would be a new conditioning experiment, not restoration of an expected input.

**Matched junction replay:** starting from the saved phase-3 coordinates, phase 4 with a soft flat-bottom restraint on all eight junction C–N distances moved A72–A73 from the previous 1.1294 Å to 1.2739 Å (native). All eight pass the existing junction check. Seven severe clashes remain, as do two non-planar amides, one Ramachandran outlier and a near-planar B73 Cα. Nearby angle and bond-length warnings at A72–74 worsened slightly rather than resolving. The 5,040 observed atoms stayed fixed. [Inspection PDB](../../tmp/1NON-inspection/12_RFD3_junction_control.pdb); [receipt](../../tmp/1non-backends/attempt01/packed103-ff/junction-replay/receipt.json).

**Local torsion repair on model 12:** three generated windows A73–75, B73–75 and B88–90 were rebuilt with 32 starting torsion sets per window and observed same-residue side-chain conformers. All 5,040 observed records and records outside those windows stayed exact. [Model 13](../../tmp/1NON-inspection/13_RFD3_local_closure.pdb) had 22 severe clashes and one failed gap junction. A second experiment restarted from model 12 with five-residue windows A73–77, B73–77 and B86–90 and tighter closure/amide penalties. [Model 14](../../tmp/1NON-inspection/14_RFD3_local_closure_expanded.pdb) passed the eight junction-length checks and removed the non-planar amide findings, but retained 40 severe clashes and four Ramachandran outliers. Neither improved the seven-clash starting model or passed the preparation geometry gate. The search ranks backbone candidates before full side-chain packing; the final full-atom score needs to enter selection before expanding windows further.

**A73 packed-selection result:** starting again from model 12, each of the same 32 A73–75 backbone proposals was individually packed and diagnosed in the complete tetramer. Eight tied the parent at seven severe clashes; 24 had more (up to 26). None reduced the count. No replacement was selected, and model 12 remains the comparison baseline. The first generated N is placed identically in all proposals by the left-anchor construction; changing downstream torsions cannot alter that N. [All packed scores, parent hash and selection result](../../tmp/1non-backends/attempt01/A73-packed-selection/comparison.json).

**Left-anchor holdout:** the same placement formula (preceding O/CA/C, 1.33 Å, 116.2°, trans-side construction) was checked against the next deposited N on 636 intact pairs of 1NON. Median error 0.019 Å, RMS 0.022 Å, maximum 0.076 Å; reversing the side gives median 2.383 Å. Projecting from A72 puts A73:N 2.205 Å from fixed A115:CB, consistent with the recurring conflict. This diagnoses the fixed planar continuation, not a measured A73 position or proof that A72 is experimentally wrong. [Per-pair errors and gap contacts](../../tmp/1non-backends/attempt01/left-anchor-holdout/receipt.json).

**A72 fixed-anchor clearance:** transferring the 636 observed C–N length/direction triples into the A72 O/CA/C frame gives at most 2.2463 Å from A73:N to fixed A115:CB. Even allowing length, CA–C–N angle and out-of-plane azimuth to vary independently across their entire observed ranges gives an analytic maximum of 2.2545 Å, below the current severe-contact boundary of 2.35 Å. With this empirical geometry envelope and immutable observations the contact cannot be removed by resampling the loop. This is a conditional case-specific refusal, not a universal peptide limit. [Clearance calculation](../../tmp/1non-backends/attempt01/left-anchor-holdout/clearance.json).

### 2026-11-19: approved hybrid reconstruction, first saved-proposal controls

The user approved the proposed investigation order: saved diffusion proposals
through a hybrid backbone/closure/packing route first, then a FrameDiPT
inpainting investigation and a separate terminal/IDR investigation. Observed
coordinate release was **not** approved. No new checkpoint sampling, MODELLER
proposal, source repair, installation, or public capability change occurred in
the first hybrid slice below. New results do not supersede historical outputs.

Private implementation:

- `scripts/hybrid_gap_reconstruction.py`: frozen source/request/ownership
  admission, ideal internal-coordinate backbone construction, optional route
  fitting, isolated worker/resource control and independent production geometry
  assessment. Only canonical, whole-residue internal gaps of 3–25 residues are
  admitted; this does not implement terminal completion.
- `scripts/hybrid_gap_rosetta.py`: the already installed operator-provided
  PyRosetta 2026.06 runtime, GeneralizedKIC and native rotamer packing with
  sequence design disabled. Non-generated heavy atoms are audited before any
  source copying. PyRosetta is a separately licensed research tool, not a new
  package dependency or redistributed runtime.
- `scripts/audit_hybrid_gap_reconstruction.py`: independent saved native/final
  readback, exact source-record custody, generated-coordinate serialization and
  unchanged request geometry assessment.

Evidence root: `tmp/hybrid-gap-2026-11-19/`. Each executed replay records its
own code/input bindings, native coordinates, candidate PDBs, resource telemetry,
diagnostics and failures. The first three engineering attempts failed before a
usable candidate (unsupported Mac address-space rlimit, then two GenKIC adapter
errors); they remain retained. `control-attempt04` established native closure
and exact observed-coordinate preservation on the small control.

Two **saved** proposals were used: Protenix seed 102 for the explicit 1UBQ
A1–72 projection with gaps A10–13/A46–49, and RFD3 seed 103 for the intact
protein-only 1NON tetramer with four gaps 73–90. Four alternatives per replay
use perturbation widths 0/10/20/40 degrees and a shared seeded RNG stream; they
are not four independent diffusion samples. GenKIC permits at most 32 attempts
per gap and stops at the first successful solution batch, packing its solutions
before selection. It is not an exhaustive 32-solution search.

The naive v1 torsion-transfer route badly displaced the long loops. V2 therefore
fits ideal-chain torsions to the **saved diffusion C-alpha coordinates**, using
an analytic Jacobian, endpoint initialization penalties and at most 200 function
evaluations per gap; exact GenKIC closure still follows. The chain is initialized
from both anchors with an internal cut. One 1NON fit reached its evaluation cap;
the receipt records that failure rather than claiming convergence. V2 also
exposed a cutpoint O-C-N angle of 88.9 degrees. V3 explicitly rebuilds only
generated carbonyl O in the closed peptide plane before final joint packing.
Native per-root preselection still precedes that explicit O rebuild, a limitation
of this prototype. No final force-field relaxation was executed.

| Saved-proposal replay | Geometry-passing alternatives | 1NON severe clash count, best of four by clash count |
|---|---:|---:|
| 1UBQ, v1 | 1/4 | not applicable |
| 1UBQ, v2 route fit | 3/4 | not applicable |
| 1UBQ, v3 carbonyl materialization | 2/4 | not applicable |
| 1NON, v1 | 0/4 | 260 |
| 1NON, v2 route fit | 0/4 | 106 |
| 1NON, v3 carbonyl materialization | 0/4 | 105 |

The raw 1UBQ proposal **already passed geometry**; these counts do not establish
an improvement. Its absent terminal OXT remains separately reported as a
preparation/readiness error, not hidden by the loop geometry result. Raw 1NON has
46 severe clash findings plus 10 coincident-atom findings under the same full CLI
diagnostic call; do not compare the sum of those categories with clash count
alone. The lowest-clash v3 alternative removes the bond-angle error but remains
rejected on clashes/Ramachandran findings, with 1.935 A generated-heavy RMS
movement relative to the saved diffusion proposal. **The hybrid does not yet
improve 1NON, let alone establish superiority to MODELLER.**

All 24 alternatives in these six replays passed native observed-coordinate
custody to exactly 0 A before source reinsertion and final observed-record
readback. Each replay has `independent-readback.json`. V3 workers used about
0.94/1.38 GiB peak aggregate RSS for 1UBQ/1NON and completed without guard stops.
Retain the prior conditional A72 clearance refusal; this work did not release
that anchor or weaken any geometry threshold.

Continuation: keep these as negative/engineering controls. Do not start a large
GenKIC sweep or call repeated packing a solution to the underlying proposal
problem. The approved FrameDiPT investigation subsequently executed below.
Terminal/IDR sampling, downstream preparation, topology qualification and any
bounded observed-coordinate repair remain distinct, unexecuted gates.

### 2026-11-19: FrameDiPT CPU controls — no accepted completion

Implemented `scripts/framedipt_research.py`, the isolated native core entry point
`deploy/framedipt/research_worker.py`, and read-only
`scripts/audit_framedipt_research.py`. No TCR annotation, remote structure/MSA
fetch, sequence design, target cropping or production-capability change is used.
The source/known-sequence/whole-gap adapter has tests for case-sensitive chains,
insertion codes, complete target masks and absence of hidden coordinates.

Pinned upstream: `instadeepai/FrameDiPT` revision
`74cf94cb7b4ed31d3a0d4230c5739d92cc061c50`. The official inpainting checkpoint
was downloaded from Hugging Face `InstaDeepAI/FrameDiPTModels` revision
`8a1f6a54676a4131e919e7975b51db16c655557b`, size 209,267,967 bytes, verified
SHA256 `b2d5768abdb0d03cfd40c4bd3fbc99cd0bda272b1f53d2b0d370eb1481790b51`.
Source and model card specify CC-BY-NC-SA-4.0; this is a private research route,
not a new distributable package dependency. A new venv overlay borrows the
existing sampling prefix's Torch 2.6.0/NumPy 2.4.6/SciPy 1.18.1, adding GPUtil
1.4.0, ml-collections 0.1.1, contextlib2 21.6.0 and six 1.17.0 only in the new
overlay. No base-prefix package installation or native upstream-code patch was
performed. Local runtime/assets are below `tmp/framedipt-2026-11-19/` and do not
transfer through Git.

**Representation limitation, not exact native atomic preservation:** FrameDiPT
consumes known residue types and N/CA/C-derived frames plus observed O-derived
psi. Source sidechains and OXT are not directly encoded. All 716 positions and
four chains of 1NON reach the native consumer, including 644 fixed residues and
72 generated residues; this is not all-atom conditioning of its 5,040 source
atoms. Of these source atoms, 2,464 are not direct frame/psi contributors. The
observed-frame-to-ideal-backbone roundtrip alone moves N by up to 0.287 A in 1NON
(0.246 A in the 1UBQ control), while C-alpha error is about 1e-6 A. The reports
explicitly fail the strict atomic-preservation claim. A later source-record
restoration is separately labelled, not presented as native preservation.

Both full-input preflights completed one actual CPU model call and passed
independent sequence/chain/mask/consumer readback. The first small engineering
attempt failed on an adapter dtype mismatch before a model endpoint; its files
remain retained. These preflights are not scientific seeds or complete samples.

Bounded inference was declared at seed 23, 100 steps, min_t 0.01, noise scale
0.1 and one CPU thread, with the existing pressure-v2/RSS-8-GiB guards:

- **1UBQ A1–72, gaps 10–13 and 46–49:** unmodified native inference completed
  101 model calls (self-conditioning plus 100 steps), with independent readback.
  Native junction C–N distances are 1.998/1.678/1.243/1.393 A; internal generated
  C–N ranges from 1.089 to 1.791 A. This is a small two-gap engineering control,
  not qualification of the model's trained single-mask 8–50-residue regime.
- **Intact protein-only 1NON:** unmodified native inference stopped after 39
  recorded model-consumer calls when fixed-frame rotation drift exceeded the
  predeclared 1e-5 component tolerance (1.01924e-5). No native endpoint exists.
  This is a numerical custody refusal, not a scored complete reconstruction.
  `deploy/framedipt/frame_roundtrip_control.py` then reproduced threshold crossing
  at conversion iteration 39 **without any model or diffusion calls**, reaching
  2.54512e-5 after 100 rotation/rotvec/tensor7 conversion cycles. The exact
  executed diagnostic source is retained in its evidence directory. This is an
  operator control, not full-model trajectory equivalence.
- **Separate 1NON fixed-reference-frame arm:** a new, explicitly declared
  `reference-projection` policy keeps the native denoiser and reverse calculation
  but restores fixed tensor7 entries after each reverse step. The tolerance was
  not relaxed, and generated entries are unchanged by the projection. It is not
  called unmodified native sampling or a continuation of the stopped run. A new
  preflight passed; every initial array exactly matches the stopped seed-23
  native arm. The complete 100-step run took 675.8 s of inference (696.2 s total),
  peak aggregate RSS 4.40 GiB, without resource-guard termination. Independent
  readback verifies all 101 calls and all 99 saved before/after projection pairs.
  Final fixed-frame translation maximum is 2.13e-6 A. Native junction distances
  are A 1.083/1.451, B 0.738/1.639, C 1.128/1.368 and D 1.391/1.412 A;
  internal generated C–N ranges from 0.964 to 1.545 A. These are backbone
  measurements, not an all-atom acceptance result.

The 716-residue intact assembly is also beyond the native base configuration's
512-residue training filter; transport success is not training-domain coverage.

`scripts/framedipt_materialize.py` separately restores original observed records
and uses the already installed PyRosetta in **pack-only** mode: no closure,
backbone minimization, observed-heavy motion or sequence design. It distinguishes
frame-derived N/CA/C/O/CB from Rosetta-constructed sidechains. The two outputs
contain 38+17 (1UBQ) and 356+200 (1NON) added atoms respectively; do not describe
all of them as neural all-atom predictions. Independent native packing/final-file
readback preserves every source record and confirms zero observed-heavy motion
inside Rosetta, without claiming zero native FrameDiPT atom drift.

| Endpoint materialization | Severe clashes | Other geometry errors | Accepted |
|---|---:|---|---:|
| 1UBQ, native backbone + packing | 10 | 1 bond-length and 2 bond-angle errors | no |
| 1NON, projected backbone + packing | 103 | 7 bond-angle errors | no |

The pre-existing missing OXT on the cropped 1UBQ source remains a separate
readiness finding. One bounded replay through the existing v3 hybrid (four
alternatives, same seed/attempt budgets, no new diffusion) also failed: 1UBQ
produced three rejected candidates and one closure refusal; 1NON produced four
rejected candidates, with best severe clash count 179 versus 103 before that
postprocess. All seven serialized hybrid candidates passed independent exact
source-record/native-packing custody. No further sweep or relaxation was launched.

Evidence: `tmp/framedipt-2026-11-19/`, especially `1ubq-sample-seed23-attempt01`,
`1non-sample-seed23-attempt01` (stopped native arm),
`fixed-frame-roundtrip-attempt01`, `1non-projection-sample-seed23-attempt01`,
`{1ubq,1non}-packed-attempt01` and `{1ubq,1non}-hybrid-attempt01`.
Per-stage plans, code snapshots, initial arrays, model-consumer traces, native
backbone endpoints, ownership ledgers, failures, diagnostics and independent
readbacks remain separate. Focused regression verification: 41 tests passed and
Ruff passed; this is not a full-suite or scientific-qualification claim.

**Conclusion:** this adapter executes and exposes native limitations, but these
two controls provide no accepted completion or superiority to MODELLER, RFD3 or
cc94. They do not establish that FrameDiPT fails on all eligible loops. Keep
public profiles/defaults unchanged. The separately approved terminal/IDR
investigation remains pending; do not label the above as terminal qualification
or silently release observed anchors to obtain a passing structure.

### 2026-11-19: IDPForge terminal investigation — official weights blocked

Inspected clean `THGLab/IDPForge` revision
`e8da98924e956fd307a6df903432e5a271637076`, retained under
`tmp/idpforge-2026-11-19/upstream`. Its source license permits educational,
research and not-for-profit use; do not treat it as blanket commercial permission.
The official README points to DOI `10.6084/m9.figshare.28414937` for weights.
Figshare's API returns HTTP 403; the landing page/DOI produces an AWS WAF
JavaScript challenge or an empty HTTP 202 response from this environment. No
challenge bypass, unofficial weight substitution or checkpoint deserialization
was performed. DataCite's public metadata is accessible: dataset title
`idp_data_weights`, updated 2025-07-08, total recorded size 2,109,243,080 bytes,
dataset rights MIT. Those metadata do **not** identify local checkpoint bytes,
supersede the source license, or establish compatibility with the pinned source.
The API/landing responses and DataCite record are retained alongside the summary.

The proposed first case remains the existing source/request pair in
`tmp/mac-reconstruction-comparison-2026-11-11/expansion-development-extension-attempt01/runs/1UBQ-A1-72-Nterm1-4-protenix-seed101/`.
Source digest and exact fixed-atom inventory were checked: 72 target positions,
68 observed residues/538 fixed atoms, missing A1–4/36 requested generated atoms.
No reference coordinates, alternative sequence, source repair, template cropping
or observation release entered this inspection. This four-residue terminal
control is not evidence that the absent residues are biologically disordered,
nor a benchmark of long-IDR ensembles.

**Native input restoration is not an output-preservation guarantee.** In
`IDPForge.recon`, each network call receives template backbone and chi torsions
at masked residues with their timestep set to zero. The denoiser aligns its
prediction and excludes the motif from updates. However, `recon` returns the
last network prediction, not the restored noised state; `sample` exports that
prediction's atom14 coordinates without a final template substitution. Source
sidechains are represented through chi features, not an exact atom-coordinate
constraint. Actual checkpoint output motion remains **unmeasured**.

`scripts/idpforge_contract_control.py` extracts the unchanged pinned `recon`
method via AST without importing the full native dependency stack. Using the
existing Torch CPU runtime, synthetic forward/denoiser consumers verified 40
calls at timesteps 199 through 4: template coordinates/torsions are restored
before every forward, while a deliberately shifted surrogate prediction remains
shifted in the return value. Exact executed code and upstream method bytes are
retained in `synthetic-recon-control-attempt01/`. This is a synthetic orchestration
control: **zero trained-model calls, zero checkpoint samples**, no chemical
acceptance, and no claim that the real network moves atoms by the chosen shift.

No IDPForge environment or sampling adapter has been provisioned. Full native
imports/strict loading, source-only torsion/secondary-structure preparation,
actual model-consumer preflight, the bounded raw CPU terminal sample and
independent geometry/observation audits remain unexecuted. Do not run the stock
`sample_ldr.py` success-count loop as a single declared sample: it may generate
replacements, relax and filter. Resume only after obtaining the official weights
and their available file metadata, recording the local checksum and verifying
strict loading. Summary: `tmp/idpforge-2026-11-19/summary.json`.

### 2026-11-19: local IDPForge weights, native terminal smoke and SS-units control

The user provided the extracted official archive at
`tmp/idpforge-2026-11-19/idpforge_dw/`. Its `weights/mdl.ckpt` is 6,525,205 bytes,
local SHA256 `a623756262138376672afcb17f2f8bd9e6780b265b57aa1b5fbd7e8b0f5d4769`.
This is user-supplied provenance and a local custody digest, not an independently
retrieved publisher checksum. Safe tensor-only `torch.load(weights_only=True)`
loaded its EMA dictionary; all 1,617,090 parameters strictly match the pinned
IDPForge architecture after the four documented OpenFold-2 IPA point-key renames.
No archive training structures, ensemble data or reference coordinates were read.
The earlier blocked summary is historical and remains unmodified; current state
is `tmp/idpforge-2026-11-19/continuation-summary.json`.

`deploy/idpforge/prepare_runtime.py` copies the pinned IDPForge packages and
OpenFold `be2ec1841f16c966c65ae0e7599ebbadc725757d` into a separate local source
tree, retaining licenses. The compatibility manifest records exactly three
patched files: Python 3.12's dataclass default-factory requirement in ESM's trunk,
and conditional CUDA-extension imports in two OpenFold modules. No attention
math, weights, sampler/reverse diffusion or native coordinate restoration was
patched. CPU uses the existing default Torch attention; requesting the absent
CUDA extension fails rather than substitutes. Clean upstream clones remain
unchanged. A new venv overlay borrows the sampling prefix's Torch 2.6.0/NumPy
2.4.6/SciPy 1.18.1 and locally adds ml-collections 0.1.1, contextlib2 21.6.0,
six 1.17.0 and MDTraj 1.11.0; no base-prefix package installation occurred.

`scripts/idpforge_research.py` admits only a single complete canonical protein
with one anchored N-terminal whole-residue gap, refusing multichain/nonprotein
context, explicit links, insertion codes/discontinuous numbering unsupported by
the native template parser, partial occupancy and incomplete atom ownership.
The source-bound control uses the existing 1UBQ A1–72 request with missing A1–4,
all 538 observed heavy atoms, and 36 missing heavy atoms. Native template
preparation prepends only the request's known missing sequence, with one stock
seed-position offset and no cropping, sequence variant or automatic repair.
An initial adapter indexing failure before a completed forward is retained in
`1ubq-terminal-preflight-attempt01`; it is not a scientific sample.

The succeeding preflight performs one actual native CPU forward. Independent
readback reconstructs the entire source-to-template atom14 transport and checks
known sequence, all 72 target positions, index axis, masks, and the exact
float32 template backbone/torsion entries received by the network. A declared
seed-42, 40-step raw sample then completed with identical preflight initialization.
All 40 consumer states (timesteps 199 through 4) and the full raw neural atom14
endpoint remain retained. Sampling took 7.83 s (31.52 s including setup), peak
aggregate RSS 570 MiB; no guard fired, no success-count loop, relaxation,
reflection, source repair, or resampling occurred.

**Concrete native template bug:** `mk_ldr_template.py` converts MDTraj phi/psi
to degrees before calling `assign_rama`, which calls `np.degrees` again. A
source-only units control found 34/68 final observed SS codes differ from the
radian-consistent preparation. This does not change observed coordinates.
The stock result remains an engineering control, not clean model qualification.
One separately declared `radian-consistent` arm corrected only this source SS
preparation, passed a new actual consumer preflight, and executed one seed-42
40-step native sample (7.52 s inference). Its source/request bytes and all five
initial coordinate/torsion/mask/noise arrays match the stock arm. The full reverse
noise stream is **not** claimed identical: the native SS encoder's coil-token
sampling/RNG consumption also differs when corrected A/B codes are present.
This is not an isolated causal estimate of the units bug. No further seeds or
cohort expansion were run.

**Global frame versus deformation:** raw network endpoints are in a different
global pose. Direct unfitted observed-heavy RMS displacements of 23.02/16.43 A
(stock/corrected) must not be described as local deformation. Separate extraction
arms retain (1) the native mapped pose and (2) one proper rigid registration of
every atom to the observed C-alpha correspondence. Neither fits individual
residues, reflects, changes internal geometry or reads held-out coordinates.
Residual observed C-alpha RMSD after registration is still 1.687/1.353 A. The
network thus does not supply exact atomic preservation even after accounting for
the arbitrary pose.

For each raw endpoint, separate **explicit source-restoration** files copy the
538 original atom records exactly and insert only the 36 neural additions; this
is not claimed native preservation. Independent readback reconstructs all four
extraction arms from the untouched native endpoints and the declared transforms,
checks inventory and generated coordinates to PDB precision, preserves all source
records, and reproduces the rejected geometry decisions.

| Raw sample, after global registration then exact source restoration | C4–N5 | Severe clashes | Other geometry errors | Accepted |
|---|---:|---:|---|---:|
| Stock native SS preparation | 3.379 A | 7 | 1 bond-length + 6 bond-angle | no |
| Radian-consistent source SS | 5.550 A | 14 | 1 bond-length + 2 bond-angle | no |

The raw full neural structures themselves also fail: stock has 41 severe clashes,
30 angle and 2 length findings; corrected SS has 73 clashes, 1 coincident-atom
finding and 2 length findings. Unregistered source-restoration controls fail as
well (13/5 geometry errors) and are not alternatives selected by the smaller
error count: their C4–N5 separations are 17.71/18.84 A. The original observed
source has no geometry ERROR findings; its missing terminal oxygen is a separate
pre-existing readiness finding in all arms. These diagnostics used the declared
Python steric fallback because Probe was unavailable, not a claimed Probe pass.

No generated-only minimization, packing, closure, automatic repair or further
checkpoint call followed these raw controls. **No accepted terminal completion
or advantage over MODELLER is established**, and this short 1UBQ case does not
qualify or rule out longer, genuinely disordered tails. Public defaults/profiles
remain unchanged. Focused verification: 56 tests passed and Ruff passed; not a
full-suite or native CUDA-compatibility claim.

Evidence includes `1ubq-terminal-sample-seed42-attempt01/`,
`1ubq-radian-ss-sample-seed42-attempt01/`, both successful preflights,
`template-rama-units-control.json`, `source-baseline-diagnose.json`, the native
endpoints/consumer states, code snapshots, and each extraction's declaration,
diagnose/assessment and independent readback under `tmp/idpforge-2026-11-19/`.

### 2026-11-19: approved source-anchored terminal reconstruction comparison

The user approved a bounded follow-up: reuse the saved terminal proposal, grow
ideal backbone backwards from fixed A5, pack sidechains and minimize only new
heavy atoms, comparing with MODELLER under identical processing. No new diffusion
sampling or observed-coordinate release was authorized or performed.

Implemented `scripts/terminal_reconstruction_research.py`, the source-only
`scripts/terminal_modeller_worker.py`, read-only
`scripts/audit_terminal_reconstruction.py`, and inspection/plotting in
`scripts/inspect_terminal_reconstruction.py`. The existing isolated Rosetta worker
adds an explicit `pack-and-cartmin` operation; earlier closure/packing modes are
unchanged. Evidence is retained in
`tmp/terminal-reconstruction-2026-11-19/attempt01/`.

**Declared comparison, no best-of selection:**

- One saved radian-consistent IDPForge seed-42 proposal from the preceding section;
  its native endpoint/extraction was independently re-audited before use.
- One new seed-42 MODELLER initial model plus one fast LoopModel derivative, with
  generated-only selections, chirality restraints, `final_malign3d=False` and no
  full-model `to_iupac` atom swaps. Its native mapped observation displacement is
  zero at PDB precision. No score-based selection or initial-model fallback.
- For **both** methods, extract proposal-internal phi/psi and the A5 phi defined
  by proposal C4/N5/CA5/C5. The undefined phi of residue 1 is explicitly unused.
  Backwards ideal N/CA/C/O construction uses only the original fixed N5/CA5/C5
  frame and those torsions, with trans peptide bonds and C4–N5 = 1.33 A before
  PDB rounding. Proposal global pose is irrelevant; no reference fit, C-alpha
  route fitting, GenKIC, local anchor movement or sequence changes occur.
- Same Rosetta seed 42, `ref2015_cart`, generated-residue-only sidechain packing,
  then one Cartesian LBFGS call capped at 500 iterations, no restarts. Explicit
  per-atom MoveMap variables are the 36 generated heavy atoms plus disposable,
  unobserved hydrogens. All 538 observed heavy atoms remain fixed. Rosetta's
  auxiliary C-terminal OXT and proline ring-closure virtual NV atoms remain fixed
  and are excluded from serialized candidate ownership together with temporary
  H. Native metadata confirms NV is virtual, not a physical heavy addition. No
  solvent is added.
- The unprocessed MODELLER proposal is separately retained and assessed so that
  worsening it with common processing cannot create an apparent IDPForge benefit.

Both workers completed without resource-guard stops (Rosetta about 10.0/5.2 s,
peak aggregate RSS 760/869 MiB for IDPForge/MODELLER). Native packed/minimizer
states, the complete Cartesian atom/mobility axis and before/after coordinates
are retained. Independent audit replays torsion extraction and backbone
construction, verifies the native atom masks and source fixation **before** file
materialization, reproduces final coordinates and geometry assessments, and
checks global L-chirality. Source-heavy maximum movement is exactly zero in both
Rosetta native records; source PDB records remain byte-identical. Both stage
budgets are limits, **not optimizer-stationarity/convergence claims**.

| Branch | Geometry ERROR findings | C4–N5 | Other geometry concern |
|---|---:|---:|---|
| Raw MODELLER, source-record projection | 0 | 1.359 A | none reported |
| IDPForge torsions, ideal growth + packing | 36 clashes | 1.330 A | rejected |
| IDPForge branch after bounded minimization | 0 | 1.348 A | cis MET1–GLN2, omega 0.4 degrees |
| MODELLER torsions, ideal growth + packing | 3 clashes | 1.330 A | rejected |
| MODELLER branch after bounded minimization | 0 | 1.348 A | none reported |

All serialized stages retain zero D-C-alpha. The original missing OXT on observed
ARG72 remains a separate terminal-readiness ERROR in every output; it was not
silently repaired or counted as a new backbone failure. Diagnose used its
reported Python steric fallback, not a claimed Probe run.

**A geometry-error pass is not a scientific win.** Although minimization removes
the clashes, the IDPForge-guided branch moves generated heavy atoms 3.980 A RMS
(maximum 6.886 A) from packing, with large torsion changes, and crosses from the
constructed trans state to an unusual cis non-PRO bond. MODELLER's corresponding
movement is 1.249 A RMS (maximum 1.956 A), without that warning. No trans wall was
added retrospectively and no further minimizer run was selected to hide this
finding. The raw MODELLER baseline was already clean. Error-level geometry gates
therefore pass after processing, but **better plausibility or superiority over
MODELLER is not established**; do not label the cis-containing result qualified.

`comparison-inspection.json` retains the warning, motion and torsion changes.
`matched-view-comparison.png` shows raw MODELLER and both processed candidates in
the same two views with the identical observed C-alpha context and marked cis
bond; it is a source-only inspection, not experimental-accuracy scoring. Final
structures are `idpforge/minimized.pdb` and `modeller/minimized.pdb`; all earlier
native/packed states remain separate. The first independent-audit invocation
found a tuple-versus-JSON-list comparison bug; wire-normalized comparison fixed
the reader, and reconstruction then reproduced exactly without rerunning science.
Both audit receipts remain distinguishable when additional readback fields are
added; no scientific output is overwritten.

Stop condition reached for this small control: chemical cleanup is feasible with
exact observations, but no benefit of the neural torsion proposal is shown. Do
not expand sampling or change public defaults on this evidence. This is still
only a four-residue 1UBQ terminal control, not qualification of longer tails,
IDRs or internal gaps. Focused verification: 60 tests passed and Ruff passed;
no full-suite claim.

### 2026-11-19: 1NON ordinary MODELLER, common gate, then authorized local release

The user reported a clashing MODELLER result marked successful, then authorized
ordinary runs first and a separate local-rearrangement experiment afterwards.
Investigation distinguishes **process completion** from **scientific acceptance**:
the GUI MODELLER Homology path currently derives `success` from exit code zero
(`gui/server/homology-api.ts`); it does not apply the diffusion bundle gate.
These research runs now report execution and quality separately. They do **not**
change GUI/production status handling or enable a new public capability.

Implemented private orchestration in `scripts/one_non_local_research.py`, an
ordinary MODELLER wrapper in `scripts/one_non_modeller_worker.py`, the isolated
`scripts/local_gap_rosetta_worker.py`, and independent
`scripts/audit_one_non_local.py`. Evidence is under
`tmp/1non-common-validation-2026-11-19/attempt02/`; `attempt01` retains a PIR
preparation failure before model generation (wildcard endpoints read only the
first chain). The existing explicit first/last-chain PIR adapter corrected that
engineering failure; no extra successful MODELLER proposal was selected.

**Ordinary stage completed before local work.** One fresh ordinary seed-101
MODELLER initial model plus one fast loop derivative saw all four chains, 716
target positions and 5,040 observed protein heavy atoms; 556 atoms belong to
the four missing 73–90 segments. Native output and logs remain unmodified.
It exited zero but has **44 severe clashes and 10 bond-angle ERROR findings**.
The same validator rejects it; `execution_status=completed` is not acceptance.

MODELLER's output also has a different global pose. Direct unfitted observed
maximum motion is 17.027 A, not a pure local-deformation metric. A separately
declared proper rigid registration on all 644 observed C-alpha atoms of the
intact assembly (no per-chain fits or held-out atoms) preserves native geometry
and gives residual C-alpha RMS 0.986 A/max 14.305 A and observed-heavy max
15.790 A. It still has the same 44 clashes/10 angle errors. Native raw, mapped,
registered and source-restored files/reports remain distinct; independent
registration readback reproduces every registered native coordinate. After
explicitly restoring exact source observations, MODELLER's proposal has 199
clashes, 6 length and 20 angle errors. The earlier unregistered source projection
is retained as a frame-transport control, not used as a primary refinement input.

The diffusion comparator is the **already retained RFD3 seed-103 packed proposal**
`tmp/1non-backends/attempt01/packed103-ff/case/proposal.pdb`, not a new sample or
an unprocessed neural endpoint. It binds the same source/request/ownership and
starts with 54 severe clashes plus one coincident-atom finding. Different native
generation histories/seeds are retained explicitly: this is matched refinement
and validation, not a statistically matched backend-superiority benchmark.

**Prospective common local policy:**

- Select two consecutive observed residues outward from each gap flank and all
  observed residues with a heavy atom within 4 A of those flank heavy atoms in
  the **original source only**. Both methods use the identical frozen mask:
  58 residues/477 observed heavy atoms, including A115; the other 4,563 observed
  heavy atoms remain strictly fixed. Generated ownership stays exactly 556 atoms.
- Per observed residue: backbone RMS <=0.5 A, C-alpha movement <=0.5 A, any heavy
  atom <=1.25 A in the original frame. Native and serialized states are both
  audited, with no clipping or global fitting of refinement displacements.
- One generated-sidechain packing plus one `ref2015_cart` Cartesian LBFGS call
  capped at 500 iterations per cell, seed 101, no solvent, no random restarts.
  Both fixed/local cells start from identical source-restored proposals, not
  from the fixed cell's result. Independent native readback confirms **identical
  packed requested-atom coordinates within each pair**.
- Both policies use the same trans-omega restraints on generated/junction
  peptides (CircularHarmonic pi, sd0.15 rad). The local cell additionally uses
  source-coordinate harmonic restraints (sd0.2 A) and flat-bound wall penalties
  (sd0.05 A). These are **soft walls, not hard inequalities**: an exceeded
  movement bound causes rejection, never an increased bound or coordinate reset.
  Disposable H may move; virtual atoms and outside-mask physical heavy atoms
  remain fixed. No heavy-atom rebuilding, chirality reflection or new sampling.
- Same gate for both methods: request-scoped expanded geometry, full-context
  chemistry/sterics at unchanged 0.6/0.9 A overlap thresholds, explicit junction
  lengths, global L-chirality, no cis non-PRO, source metadata and movement bounds.
  A low clash count or engine exit code alone cannot pass. Nonfatal warning
  counts remain visible; sterics used the reported Python fallback, not Probe.

Fixed controls ran and were assessed first. Both first local launches then failed
**before packing/minimization**, because this Rosetta binding exposes `BoundFunc`
under `core.scoring.constraints`, not `core.scoring.func`. Those failed directories
are retained. Only the local cells were retried under `local-retry/`, with a new
code-bound declaration correcting that namespace, unchanged mask/seeds/budgets,
and no repeated fixed run, MODELLER generation or diffusion sampling.

| Completed refinement cell | Severe clashes | Other full-context geometry errors | Movement bounds | Final quality |
|---|---:|---|---|---|
| MODELLER proposal, all observations fixed | 87 | 30 angles, 8 lengths, 1 chirality | fixed exactly | rejected |
| Same MODELLER proposal, bounded local policy | 41 | 15 angles, 8 lengths | exceeded | rejected |
| Retained packed RFD3, all observations fixed | 7 | 3 angles | fixed exactly | rejected |
| Same RFD3 proposal, bounded local policy | 1 | explicit B90–B91 junction still fails | within limits | rejected |

The MODELLER-fixed D-C-alpha is a failure of this common refinement branch, not
a claim that MODELLER itself emitted that centre; it is retained privately and
never accepted or silently reflected. MODELLER-local is globally L but exceeds
bounds at A111, A72, B91 and C68: maximum backbone RMS 0.719 A and heavy movement
1.398 A (maximum C-alpha 0.471 A). It fails both quality and movement gates.

**The local RFD3 result is a real improvement but still rejected.** Its maximum
native observed backbone RMS is 0.347 A, C-alpha motion 0.239 A and heavy movement
0.565 A; outside-mask movement is exactly zero. Angle/nonplanar-amide issues are
removed, but ASP74:O–VAL77:CG1 still has a severe 1.39 A vdW overlap (distance
1.83 A), and B90:C–B91:N is 1.5066 A and fails the unchanged junction gate.
Generic gate-summary issues are not counted as extra distinct junctions. It has
no D-C-alpha or cis non-PRO warning. This does not establish broader loop
qualification or authorize publishing this candidate as accepted.

All four completed workers exited zero; **all four scientific results are
rejected**. Native outside-mask maximum movement is zero in every cell.
`independent-readback-v2.json` independently checks original masks, packing,
native before/after states, source metadata, PDB-coordinate transport, bounds and
recomputed request geometry. Each completed cell took about 51–58 seconds under
the one-thread/RSS-8-GiB/pressure-v2 guard, with peak aggregate RSS about 0.9–1.2
GiB and no guard stops. Iteration caps are not stationarity/convergence claims.

Primary outputs: `summary.json`, `registered-baseline/ready.json`,
`{modeller,rfd3}-fixed/result.json`,
`local-retry/{modeller,rfd3}-local/result.json` and their rejected `candidate.pdb`
files. No further mask expansion, threshold change, minimizer continuation or
sampling was performed. Focused verification: 68 tests passed and Ruff passed;
no full-suite or public GUI-fix claim.

### 2026-11-19: one additional local minimization — zero severe clashes, junction rejected

At the user's request, performed **one** additional minimizer call on the best
bounded RFD3-local result above. This is not another packing proposal or diffusion
sample. `scripts/polish_one_non_local.py` declares the continuation and independent
audit; the existing native worker has an explicit resume-only 1,000-iteration cap
while the original fixed/local modes retain their 500-iteration caps.

Evidence is under
`tmp/1non-common-validation-2026-11-19/attempt02/polish-rfd3-local-1000-attempt01/`.
The original source positions, 58-residue/477-observed-heavy local mask, 556
generated atoms, restraint strengths, trans-omega restraints, `ref2015_cart`
objective and all quality/movement thresholds are unchanged. Positional
restraints remain referenced to the **original source**, not the already moved
parent. No ordinary whole-complex `dvbfixer minimize`, solvent, new neural call,
repacking, coordinate clipping, chirality reflection or global observation
release occurred.

The worker reconstructs the native particle axis and restores all **11,436**
coordinates (including disposable H and virtual atoms) from the parent's retained
full-precision returned state. Coordinates and weighted score match exactly
before optimization: parent and resumed score are both 206.90000422132846 REU.
This excludes a hidden new H/packing state as the cause of improvement. The new
LBFGS instance resets optimizer history; it is not a restored Hessian checkpoint.
One call capped at 1,000 additional iterations took 100.8 s total, peak aggregate
RSS about 1.09 GiB, with no resource-guard stop. No stationarity claim is made.

| Metric | Parent local result | After extra minimization |
|---|---:|---:|
| Severe clashes | 1 | **0** |
| Full CLI geometry ERROR findings | 1 | **0** |
| Explicit B90:C–B91:N junction | 1.5066 A | **1.5127 A**, still rejected |
| Maximum observed per-residue backbone RMS | 0.3465 A | 0.3511 A |
| Maximum observed C-alpha movement | 0.2394 A | 0.2435 A |
| Maximum observed heavy movement | 0.5648 A | 0.5674 A |
| Native outside-mask movement | 0 | 0 |
| Weighted Rosetta score, same objective | 206.9000 | 98.7364 |

All movement bounds still pass, global C-alpha chirality remains L and there is
no cis non-PRO finding. Sixteen steric WARNING findings and four bond-length/
three angle warnings remain visible. The request-scoped explicit junction check
is stricter than relying on full CLI ERROR counts alone; its unchanged C–N
range is 1.20–1.45 A, so it still rejects B90–B91. The generic junction-summary
issue is not a second physical defect.
The candidate therefore remains **`quality_status=rejected`**, despite zero
severe clashes and process exit zero. No threshold was relaxed to accept it.

Read-only native topology inspection confirms all eight gap junctions have their
polymer connections, including B90–B91: this is not an omitted covalent bond at
import. Lowering the current energy did not enforce the required junction length;
blindly increasing iteration counts is not established as a remedy.

`resume-audit.json`, `independent-readback.json` and
`native-junction-topology-readback.json` retain state/objective, mask, native and
serialized custody, global L, independently recomputed sterics and final request
geometry. `candidate.pdb` is an inspection artifact, not an accepted publication.
The earlier four-cell summary stays unchanged. Focused verification: 70 tests
passed and Ruff passed; production/GUI defaults remain unchanged.

### 2026-11-19: targeted B90–B91 length restraint — first passing local 1NON control

The user authorized the proposed explicit junction-length control. Starting from
the immediately preceding polished RFD3-local native state, one additional
Cartesian minimizer call (cap 1,000 iterations) added **only** an AtomPair
Harmonic restraint on B90:C–B91:N, target 1.33 A, sd0.02 A, weight 1. This is a
strong soft restraint, not an exact holonomic distance constraint. It is not an
adaptive force/weight/seed sweep. The native polymer connection and both
endpoints' existing mobile permissions were verified before adding it.

`scripts/polish_one_non_local.py --junction-control` and the isolated native worker
retain the original 58-residue/477-observed-heavy mask, original-source positional
references, backbone/CA/heavy movement limits, trans restraints, complete 556-atom
generated ownership and base `ref2015_cart` objective. No observation re-anchoring,
packing, new H state, new diffusion/MODELLER proposal, source repair, clipping or
outside-mask release occurred. Other modes' constraints/budgets are unchanged.

All 11,436 native coordinates matched the polished parent exactly at resume. The
base score also matched exactly at 98.73644083254612 REU **before** the new term
was introduced. Independent readback verifies that augmented-minus-base energy
both before and after optimization is precisely the declared harmonic penalty;
an objective change is not disguised as an ordinary energy decrease.

**This case now passes the unchanged common geometry and movement gate.**

| Metric | Polished parent | Targeted junction control |
|---|---:|---:|
| Native B90:C–B91:N | 1.5124 A | **1.3334 A** |
| Serialized B90:C–B91:N | 1.5127 A | **1.3327 A** |
| Allowed junction interval | 1.20–1.45 A | unchanged |
| Severe clashes / full CLI geometry ERROR findings | 0 / 0 | **0 / 0** |
| Request geometry issues / admission failures | junction failure | **empty / empty** |
| Maximum observed per-residue backbone RMS | 0.3511 A | **0.3501 A** |
| Maximum observed C-alpha movement | 0.2435 A | **0.1919 A** |
| Maximum observed heavy movement | 0.5674 A | **0.5674 A** |
| Native outside-mask movement | 0 | **0** |
| Base score, excluding the new distance term | 98.7364 | 85.5312 |

No D-C-alpha or cis non-PRO finding remains. Sixteen steric warnings, four
bond-length warnings and two angle warnings remain visible; passing the declared
gate does not mean a warning-free structure. Native augmented energy is
181.9394 before and 85.5595 after the call; it must not be directly compared with
the parent's unaugmented energy without separating the new penalty. The worker
completed in 54.6 s, peak aggregate RSS about 1.06 GiB, without a guard stop.
Iteration caps and successful geometry checks are not stationarity claims.

Evidence:
`tmp/1non-common-validation-2026-11-19/attempt02/junction-rfd3-local-1000-attempt01/`.
`candidate.pdb` is the passing **private local-motion control**; `summary.json`,
`result.json`, `junction-length-restraint.json`, full diagnostics,
`resume-audit.json` and `independent-readback.json` separate execution, objective
change, native state, metadata/custody, sterics, global L and final acceptance.
The full request geometry was recomputed from serialized bytes. All ancestors,
including the zero-clash-but-rejected parent, remain unchanged.

Scope is deliberately narrow: one retained packed RFD3 proposal, one protein-only
tetramer, an adaptively targeted standard peptide length, and explicitly allowed
local observed motion. This does **not** restore the original all-observations-
fixed claim, qualify diffusion generally, prove superiority over MODELLER or
establish stability after removing the extra restraint. No unrestrained follow-up
or additional minimizer call was executed. Production defaults and GUI behavior
remain unchanged. Focused verification: 72 tests passed and Ruff passed.

### 2026-11-19: larger unaugmented minimization budget does not repair B90–B91

At the user's authorization, returned to the saved **pre-junction-restraint**
polished parent, not the accepted length-restrained candidate. Declared three
successive Cartesian MinMover calls, each capped at 3,000 iterations, with the
same `lbfgs_armijo_nonmonotone` tolerance 1e-6. All 11,436 starting native
coordinates and the base score matched the parent. Original source anchors,
atom ownership, local mask and movement/geometry limits remained unchanged;
there was no new atom-pair restraint, packing or sampling. Coordinates carry
between calls, while LBFGS history restarts. This is a bounded restarted-budget
control, **not** one uninterrupted 9,000-iteration trajectory.

| State | Base energy (REU) | Native B90:C–B91:N (A) | Maximum atom movement during call (A) |
|---|---:|---:|---:|
| Parent | 98.736440832546 | 1.5124312711 | — |
| Pass 1 | 98.736436273764 | 1.5124312040 | 0.000002438 |
| Pass 2 | 98.732180665146 | 1.5124505099 | 0.042475643 |
| Pass 3 | 98.732179758311 | 1.5124505076 | 0.000000153 |

All three independently validated serialized checkpoints remain rejected by the
same 1.20–1.45 A junction gate. Final serialization still gives 1.5127 A.
There are zero severe clashes/full CLI geometry errors, global L, no cis non-PRO,
and 16 steric/four length/three angle warnings. Native outside-mask movement is
zero; maximum observed backbone RMS/CA/heavy movement is
0.351066/0.243546/0.567464 A, within the original limits.

Total worker time was only 8.23 s, without a guard stop. **3,000 is a cap, not a
measured iteration count.** The installed MinMover binding reports score change
and MS_SUCCESS but not actual iteration count or a native termination reason.
Those fields are explicitly unavailable, and no native convergence/stationarity
claim is made. The tiny final-pass energy change (9.07e-7 REU) and coordinate
change support an operational plateau; increasing the allowed budget/restarting
did not repair this bond. This weakens the simple budget-exhaustion explanation,
but does not prove the exact balance of force-field and positional-restraint
terms, global optimality, or behavior under a different optimizer/tolerance.

Evidence is in
`tmp/1non-common-validation-2026-11-19/attempt02/budget-rfd3-local-3x3000-attempt01/`:
`pass-{1,2,3}.json` retain native states and progress, corresponding
`pass-*-candidate.pdb`/`pass-*-validation.json` retain common validation,
`summary.json`, `resume-audit.json`, and `independent-readback-v2.json` with its
archived auditor bind the unchanged source, native fixation, movements, global L,
sterics and rejection. No extra optimization occurs during readback. The passing
targeted-junction candidate and all earlier evidence remain unchanged. Focused
verification: 73 tests passed and Ruff passed; no production/default changes.

**Why cc94:** its working MPS connector let us test complete-tetramer transport
and native scaffold fixation. No head-to-head result established it as the best
model for these loops. At the [first denoiser call](../../tmp/1NON-inspection/details/connector/report.json),
hiding chain B's crop coordinates at unchanged noise altered predictions for
the gaps on A, C and D too. Sensitivity to an individual late-stage steric
contact such as B64–B90 remains untested.

**Inspect the coordinates:** [local inspection folder](../../tmp/1NON-inspection/).
`01` is the protein-only input; `02` and `03` are separate raw cc94 outputs;
`04` and `05` are two refinements of **`02`**. See `sources.tsv` there for the
exact filenames, original paths and checksums.

**Local comparison completed:** the four source-defined components were sampled separately,
so context size and joint versus separate sampling both differ. Each local observation
passed the native fixation audit before assembly (maximum error 5.18e-6 Å or less;
1,501 projection payloads per arm). All 556 generated atoms were scored together;
all 5,040 observed PDB records are exact in that output. Of 563 severe loop–scaffold
contacts, 349 involve scaffold atoms present in the corresponding motif and 214
involve omitted atoms. The other 27 severe contacts involve only generated atoms,
including nine cross-chain contacts. Full-input baseline: 328 loop–scaffold and
29 generated–generated. [Local readback and all eight junctions](../../tmp/1non-local-motif-2026-11-18/attempt01/assessment/receipt.json);
[assembled raw inspection PDB](../../tmp/1NON-inspection/06_raw_local_motifs.pdb).

**Open backend comparison:** compare models on the
same protein-only input after verifying the actual encoded chains, positions
and atom mask for each. Score raw outputs first: native observed-atom drift
before any copying, generated–scaffold clashes and C–N distances at both ends
of each gap. Report omitted observations in a local arm as omitted, not natively
preserved. Keep refinement outcomes separate.

**Historical calibration status:** the finite anchor-mobility study in section 3
superseded the older autonomous next-case continuation. Do not resume
unrelated reconstruction batches. The frozen **100-cell grid was launched on
2026-11-15** after four-case full-context preflight, adapter verification and a
green fast suite. All 100 cells and independent readback are now accounted for;
the first saved-panel trade-off map is complete. The original protein-adapter
gate flags have a scope mismatch and are not a scientific success-rate table;
use the separately bound corrected readback below. **The subsequent failure
investigation found polar-H collapse, nonstationary returned states and a steric
exclusion bug. The map is descriptive evidence of that executed setup, not an
unconfounded calibration of anchor-policy effectiveness.** Do not start a new
cohort or choose movement limits until those numerical/diagnostic issues have
separately declared controls.
**Diagnostic fix is now available:** `219ba8b` (PR #27, main merge `9cdd09b`)
replaces proximity-created exclusions with supplied/template connectivity and
reports unresolved chemistry as partial. Research integration is `17b46b4`.
The user authorized a read-only v2 revalidation of the retained inputs/results;
its status and evidence are recorded below. This does not authorize new
minimization, H placement, source repair or a new cohort.
**Subsequent user instruction:** begin implementing the next handoff step.
Private immutable-context/returned-state guards and six saved-state controls are
now implemented below. Direct OpenMM X–H constraints remain incompatible with
massless observed parents. The user subsequently authorized a reduced-coordinate
route that eliminates fixed observations and covalent X–H lengths from the
optimization variables; its bounded implementation controls are declared below.
No observation release, soft-penalty substitution or production-default change
is permitted. Launch is not evidence of completion or successful reconstruction.
3C7D/1HE8 remain blocked, not executed or replaced. Use run/readback receipts,
not process success or this status line, to establish individual outcomes.

## 1. Purpose and source data

We work with **experimentally determined structures deposited in public structural
databases**, principally the Protein Data Bank. The task is to complete absent
coordinates for residues whose identities and order are already known from the
deposited structure's sequence records and associated experimental construct
information. We are not designing new amino-acid sequences or new proteins.

The distinction is between a known sequence and incomplete experimental
coordinates: an experimentally studied construct can contain residues that were
not resolved in the deposited coordinate model. We estimate coordinates for those
residues while preserving the observed scaffold and relevant structure context.

Before completion, bind the accession, original file digest, assembly, chain and
insertion-code identities, construct sequence, sequence-to-coordinate mapping and
exact missing-residue mask. A database-reference difference, engineered deletion
or omitted domain is not automatically an unresolved coordinate gap. If the
deposited construct does not establish the requested residue identities, report
the ambiguity rather than inventing sequence or borrowing an unverified canonical
sequence.

**Coordinate completion is not sequence generation.** ML backends do generate
candidate atomic coordinates; these are computational proposals for the fixed
known sequence, not new experimental observations. Do not describe predicted
missing coordinates as experimentally measured or claim that no coordinate
generation occurs. Existing fields such as `generated_atoms` retain their spelling
and mean computationally completed atom ownership.

## 2. Tools and boundaries

The workflow uses existing structural-biology tools and publicly documented
modeling methods, including open-source tools where applicable. Research
diffusion backends include RFD3, Protenix, Boltz and Protpardelle. MODELLER appears
only in retained historical comparator records and is **not an active arm**.
Record the
actual tool, revision, checkpoint and license for each executed arm; not every
backend is applicable to every structure.

Open source code, downloadable weights and permission to redistribute them are
different properties. MODELLER has separate licensing requirements; checkpoint
terms must also be checked independently. Do not label the entire stack freely
licensed or treat a public structure accession as evidence of model accuracy.

Current scope:

- Complete known-sequence missing residues in retained experimental structures,
  focusing on loops/local omissions in a sufficiently observed scaffold.
- Preserve observed protein coordinates and original full-source provenance;
  under the latest protein-only policy, remove nonprotein context explicitly
  before a new request is frozen and do not attach it back after reconstruction.
  Keep computational additions distinguishable from observations through
  atom-level provenance. Historical full-context contracts remain unchanged.
- Check structural consistency, geometry, source preservation and reproducibility.
- Do not introduce sequence design, mutations, novel sequence combinations,
  binding/activity optimization or biological-function engineering into this task.
- Domain-scale omissions, 7N8X and the four already deferred focused submodels
  remain deferred. Do not silently crop them into eligible tasks.
- Do not expand to synthesis, expression or biological testing. None is part of
  this coordinate-restoration handoff.

### Benign coordinate-restoration scope; excluded biological-threat objectives

**The user explicitly reaffirmed the project purpose:** work with public
experimental structures and already known sequences to restore absent structural
coordinates. We are **not designing or generating new protein sequences**, new
biological entities or enhanced biological functions. Diffusion generates proposed
coordinates for the existing sequence; those coordinates are modeled, not new
experimental observations. Do not misdescribe that operation as novel-protein
design or assume harmful intent from the use of a generative coordinate model.

Biological-threat development, harmful pathogen/toxin enhancement, sequence or
function engineering, synthesis, expression, dissemination and weaponization
are explicitly outside the task. Protein-only cleanup is solely an adapter/input
scope decision; it does not alter the source object's identity or biological role.
Preserve accession, sequence, organism/function metadata and original provenance.
Do not relabel targets or delete metadata to bypass a concrete concern.

These exclusions describe **non-goals**, not a finding of a biological threat in
the current work, and do not create another mandatory case-review queue. Do not
invent a case-specific threat or block ordinary structural/numerical work merely
because proteins or diffusion are involved. Conversely, public availability and
an unchanged sequence are not an absolute proof of zero risk for every possible
future operation. If a concrete threat-enabling concern actually arises, state
it specifically, pause that operation and keep the continuation non-enabling;
do not claim blanket clearance or silently change the scientific task.

## 3. Current direction and settled state

### Active continuation: diffusion only, protein only, explicit whole-scaffold input

The user clarified that the objective is **diffusion**, not MODELLER. Do not run
new MODELLER candidates, tune its reconstruction behavior or keep trying to
rescue its retained proposals as a substitute for a diffusion study. Keep the
old experiments as historical numerical/comparator evidence only. In particular,
the old **18/100** map contains MODELLER-derived starting proposals for 1NON/2BE1
and common-refiner derivatives; it is not a diffusion-only success rate. Its
18 passing cells derive from one 1KTZ proposal, not 18 independent diffusion runs.

**Protein-only policy, explicitly authorized by the user:** remove waters,
nonprotein ligands and ions **before** generation, request freezing and any
subsequent refinement. Do not attach them back to an output from a model that
never saw them. Keep the original input plus an exact removal ledger separately;
claims now concern protein-only restoration, not full-complex compatibility.
Recognized protein residues written as HETATM are not removable just because of
the record tag. Ambiguous modified protein or protein–nonprotein covalent links
need an explicit chemistry decision or an input refusal, not silent bond cutting.
Retain case-sensitive chain IDs, observed atom records, occupancy and ANISOU.
Freeze a **new source-bound request** for the new projection; do not mutate or
reinterpret an old full-context request, proposal or force-field snapshot.

`scripts/prepare_protein_only_reconstruction.py` is an initial, **not yet
production-qualified**, private preparation helper. It now has seven focused
tests and was executed on 1NON as recorded below. Do not infer a production
feature or a completed cleaned cohort from this single-case preparation.
Production CLI/API/GUI defaults remain unchanged.

**What actually happened for diffusion on 1NON:** the retained Protenix seed-7
attempt used `anchor-shell-v1`, local partner context. The full target was four
179-residue chains (716 target residues), with one 18-residue internal gap
(73–90) per chain. The encoded native input represented **354 residues as 14
fragment entities**, taken from all four chains: 2,193 fixed observed atoms were
represented conditioning and **2,847 fixed observed atoms were preserved-only**,
outside the model's conditioning. It was **not the intact tetramer as four whole
chains**. The run hit the declared MPS driver-memory limit and produced **no
candidate**; its geometry is not assessable. Do not confuse this resource outcome
with failure to reconstruct the complete tetramer or with the later MODELLER
geometry failures.

**Historical clarification:** the statements above refer specifically to the
seed-7 attempt. A later seed-11 Protenix run on the **same incomplete 354-residue
context** did retain a raw sampler return before a downstream resource stop;
its independently audited raw geometry failed (junction C–N up to 3.9705 Å).
See the retained
[seed-11 raw readback](../../tmp/mac-real-gap-continuation-2026-11-12/sampler-return-study-seed11-readback-attempt01/raw-readback.json).
Do not erase that actual diffusion endpoint, or treat it as an intact-tetramer
test. MODELLER-derived failures remain a different lineage.

Evidence:
- [Actual encoder readback: 14 entities / 354 residues](../../tmp/mac-real-gap-continuation-2026-11-12/retained-loop-protenix-encoder-attempt01/1NON/encoder-readback.json).
- [Source-bound request and represented/preserved-only atom roles](../../tmp/mac-real-gap-continuation-2026-11-12/retained-loop-protenix-encoder-attempt01/1NON/context-request.json).
- [Independent historical sampling outcome: resource stop, no candidate](../../tmp/mac-real-gap-continuation-2026-11-12/retained-loop-context-sampling-final-readback-attempt01/readback.json).

**Next bounded implementation task:** identify a supported diffusion route that
actually consumes the complete protein-only 1NON tetramer, carries all four gaps
and the full intended observed scaffold, and respects native fixed-coordinate
semantics. Record the real model axes/chain representation and any size/resource
limits before launching. If full input is unsupported, report that blocker;
do not silently replace it with local fragments, monomers, MODELLER output or a
larger resource allowance. Preserve/read back raw model observations **before**
any source-coordinate restoration; copying source atoms back is not evidence of
native fixation. Any intentionally permitted observed-anchor movement still
needs its own explicit scope and movement audit.

Within that benign restoration scope, after adapter/data preflight freeze a finite
backend/checkpoint/input/seed/budget manifest before actual diffusion sampling.
Report raw diffusion proposals separately from optional generated-only refinement
and all failures. No new source recruitment, broad cohort, silent fallback,
terminal-chemistry work or production capability change is implied. The private
reduced-coordinate refiner is supporting numerical infrastructure, not the active
research objective or evidence that diffusion succeeded on 1NON.

### Executed full-tetramer diffusion smoke: 1NON, cc94/MPS, seed 23

The user's instruction to proceed was implemented as one finite, source-bound
protein-only **diffusion** control, not another MODELLER/refinement study.
`scripts/run_1non_diffusion_full.py` declares the run and
`scripts/audit_1non_diffusion_full.py` performs native-first readback.

- Protein-only preparation removes **358 HOH atoms** and retains all **5,040
  source protein atom records** exactly, with the full original and removal
  ledger separate. No ligand/water context is attached back. A fresh request
  preserves the same known target sequences and protein fixed/generated masks.
- Actual complete-context cc94 input: **716 positions, four complete 179-residue
  chains, four 18-residue gaps, 556 generated heavy atoms**, all 5,040 fixed input
  atoms included in conditioning and zero output-only observations. The explicit
  private complete-context route is used; the public 512-token profile is not
  widened, and the intermediate public local plans are not the sampling input.
- Native reader transport is tested before sampling, without loading weights:
  four chains / 644 observed residues / 5,040 fixed atoms, maximum input decode
  error **4.20e-6 Å**. This is input representation, **not proof of native output
  fixation**. The existing device guard was narrowly extended to count actual
  `protpardelle.*` MPS tensor calls, retaining all original resource limits.
- Exactly **one cc94 model call, seed 23, 500 ODE updates**, scale 1, churn 0,
  native chain/relative positions, no fallback. Execution completed in about
  **35.5 seconds** without a resource stop (aggregate RSS about 0.83 GB,
  native MPS-driver peak about 1.43 GB). No MODELLER or refinement ran.

**Result: actual full-tetramer sampling executed, but native observation
preservation failed.** After one global rigid-frame fit, the **5,036 atoms
represented in the native canonical residue output mask** have RMS displacement
**0.2750 Å**, maximum **3.2131 Å**, versus the declared 1e-4 Å float-transport
tolerance. Backbone RMS is **0.2215 Å**; Cα maximum displacement is **2.5012 Å**.
This is not a rounding issue or solely a terminal-oxygen problem. Full scaffold
conditioning in this adapter is not sufficient evidence of hard fixation.

The four input terminal **OXT atoms are read by the motif reader but absent from
the native canonical residue output mask**, despite atom37 having an OXT column.
They are recorded as missing native observations, not interpreted as predicted
coordinates or silently restored. The first audit mistakenly included those
unused native columns in its displacement statistic (0.2968/4.3940 Å); preserve
that report. The mask-aware read-only audit corrects the measurement and still
rejects the observation contract. No resampling occurred, and the meaningful
non-OXT/backbone motion alone independently causes refusal.

The run is therefore **not an accepted restoration**. No source-restored candidate
or generated `.dat` was written after the failed native gate, and no minimization
or extra seed was launched. It does not show that 1NON is intrinsically hard or
unreconstructable; it establishes a concrete missing guarantee in this native
conditioning/adapter route. The next implementation issue is native fixed-atom
handling (including actual residue output-mask semantics), not more MODELLER
minimization or claiming success after copying observations back. Any sampling
continuation needs a separately frozen bounded declaration.

### Private fixed-scaffold inpainting: implementation and declared control

The user authorized the next adapter implementation and bounded control. Source
inspection confirmed that native `replacement_guidance` only substitutes the
denoiser input; it does not enforce the predicted clean coordinates or the ODE
state. Merely enabling that option is not a fixed-scaffold contract.

`deploy/protpardelle-1c/fixed_scaffold_inpainting.py` implements a separate private
VE/Euler intervention. For known atom indices, maintain
`x_known(sigma) = clean + sigma * epsilon`, taking `epsilon` from the model's
existing initialized Gaussian tensor without another RNG draw. Project the
denoiser's known clean prediction to `clean`, correct its normalized structural
self-conditioning, and enforce the known noise path on initialization, each
denoiser input and each accepted ODE update. At zero noise, the native known
coordinates equal the internal motif reference exactly. Intermediate noisy
coordinates are **not** misreported as fixed clean observations: their expected
noise-dependent state and the invariant clean reference are audited separately.

Only declared fixed atoms are projected. Generated coordinates and unowned dummy
channels remain bit-identical across each projection. The four existing input
OXT atoms are carried as fixed input data through the native process; they are
not newly generated atoms or a claim that cc94 learned an OXT output capability.
The native output mask explicitly records canonical predictions **plus fixed
input atoms**, rather than interpreting unused atom37 columns as predictions.

This intervention is compiled as a scoped **model-instance method override**
after the existing pinned runtime preflight. The exact original/patched method
hashes and executed sampler source are retained. No installed engine file,
checkpoint, public profile, frozen original request or production default is
edited. Unsupported schedules, sequence-design modes, churn, recurrence, masking
or tensor-axis changes fail rather than receive an implicit fallback. The scope
is fixed-sequence, no-churn, uniform VE/Euler, one recurrence, scale 1.

`scripts/audit_cc94_inpainting.py` independently checks every retained input,
clean-prediction and accepted-step projection in NumPy: actual atom ownership,
the reference/noise relation, unchanged generated values, a decreasing noise
schedule ending at zero, and the raw native endpoint before any source copying.
Compact before/after arrays are retained at every such boundary; source-frame
and final geometry checks remain separate in the native-first auditor.

Pre-execution controls: five stdlib-unittest cases passed in the isolated real
Torch runtime, including CPU and MPS numerical paths, exact known-state and
generated-bit checks, noise consistency, OXT input-only semantics, mask/order
refusals, scoped installation/restoration against the actual pinned sampling
method, and an independent audit that detects a corrupted frame. No checkpoint
weights were loaded by those tests. A first import needed the existing model
parameter directory declared in the test child environment; missing unused
ESMFold/MPNN/Foldseek assets remain visible warnings, not installation changes.
Seventeen host regression tests and Ruff also pass.

Exactly **one** new full-1NON run is declared, with the same protein-only input,
four chains / 716 positions / four gaps, checkpoint, **seed 23**, 500 updates,
scale 1 and churn 0 as the preceding control. The changed factor is the explicit
known-state inpainting process; this is a new method intervention, **not an
independent seed, old-run restart or silent retry**. No refinement, new proposal
backend, extra source, wider geometry threshold or larger resource limit is
authorized. Public capabilities remain disabled.

- [Frozen inpainting control manifest](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/plan.json).
- [Native reader verification for this control](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/native-reader.json).
- [Isolated native unit-test record](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/software/native-unit.json).

**Control completed, including independent trajectory/native-first readback.**
Exactly one new cc94/MPS call completed all 500 updates in about **56.8 seconds**,
without a resource stop (aggregate RSS about 0.96 GB; MPS-driver peak about
1.43 GB). The installed engine files and frozen input/checkpoint remain unchanged.

The independent audit checked **1,501 compact projection payloads**: initialization
plus all 500 denoiser inputs, 500 clean predictions and 500 accepted ODE steps.
There are also 500 self-conditioning events. All 556 generated atoms are
bit-identical across each projection; the NumPy recomputation agrees with the
declared known-noise path exactly for this run. At zero noise, the returned
native known coordinates equal the internal reference **bit-for-bit**. This
is verified from retained arrays, not inferred from an optimizer flag or a
source-restored output.

All **5,040 fixed input atoms**, including the four input-only OXT atoms, are
present under the explicit native output-mask contract. Before any source-record
copying, a single global rigid-frame fit gives known-atom RMS/max error
**2.588e-6 / 7.285e-6 Å**, below the unchanged 1e-4 Å transport tolerance. The
Cα maximum error is **6.253e-6 Å**. Missing native observations: zero. The known
scaffold-fixation mechanism therefore passes this case-specific implementation
control; the earlier plain-conditioning refusal remains separately retained.

**Raw loop geometry still fails.** There are no observation-admission failures,
and the full CLI reports global L-Cα, but mask assessment retains seven bond-length
findings, one amide finding, one Ramachandran outlier and **357 severe pairwise
clash findings**. Full CLI additionally exposes four severe heavy coincidences,
including D82:CA–D82:CE1 at **0.274 Å**. These are real raw generated-geometry
problems, not OXT bookkeeping or a failure to supply the tetramer. Pairwise clash
counts are not counts of independent residues or independent experiments.

All eight serialized junction C–N distances (Å):

| Chain | C72–N73 | C90–N91 |
|---|---:|---:|
|A|1.3835|1.3008|
|B|1.4973|1.5547|
|C|1.4466|1.4060|
|D|1.5327|1.5028|

The existing junction gate flags four of those junctions (B/D); neither its limits
nor the other geometry thresholds were widened. A private raw protein-only PDB
and generated-only `.dat` (556 atoms) were written **only after** the native
fixation gate passed, retaining all observed source records exactly. They remain
**rejected raw outputs**, not accepted/public structures or donor approvals.
At that inpainting checkpoint, no refinement, new seed, extended schedule or
further model call had been executed. The subsequent authorized refinement is
recorded separately below.
Post-execution checks: **82 focused host tests passed**, in addition to the five
isolated native Torch tests above; agent documentation and all 62 local handoff
references validate. This is not a new complete scientific-stack suite result.

Evidence:
- [Native-first outcome and source-frame fixation metrics](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/independent-readback/receipt.json).
- [Independent full trajectory/projection audit](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/independent-readback/inpainting-audit.json).
- [Raw mask geometry findings](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/independent-readback/assessment.json).
- [Complete raw diagnostic report](../../tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/independent-readback/full-cli.json).

**Continuation scope:** assess generated-only geometry repair while retaining this
proven fixed-scaffold mechanism and the failed raw candidate. Do not return to MODELLER,
move observations, call fixation success a successful reconstruction, or conceal
raw defects with source copying. One fixed-scaffold seed is not backend
qualification or evidence that the native loop coordinates are accurate.

Retained evidence (local payloads must be transferred explicitly between hosts):
- [Protein-only preparation and exact removal ledger](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/preparation/preparation.json).
- [Complete-tetramer axis/data admission](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/planning/data-admission.json).
- [Frozen one-call execution manifest](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/study/plan.json).
- [Native reader transport](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/study/native-reader.json).
- [500-update native sampler trace](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/study/native-sampling/trace.json).
- [Final mask-aware native-first refusal receipt](../../tmp/1non-full-tetramer-diffusion-2026-11-17/attempt01/study/independent-readback-final/receipt.json).

Preparation/wrapper/guard/native-decoder checks pass (15 focused tests); the
combined focused regression set passed **80 tests**, and Ruff/agent-documentation
checks passed. Original
execution code is archived under the frozen plan; the corrected auditor is
archived beside its final receipt. Historical reports and model artifacts are
unchanged. This case-specific research does not enable a public capability or
change production defaults.

### Generated-only refinement of the inpainted 1NON candidate

The user's next authorization is for **one** CPU refinement of the saved diffusion
candidate, not another diffusion proposal. `scripts/refine_1non_generated_only.py`
binds the successful native-fixation receipt, the exact raw candidate, complete
four-chain/716-position request, source, generated-only `.dat`, implementation and
installed OpenMM force-field XML. It refuses changed observations, foreign atoms,
incomplete canonical heavy atoms or inconsistent ownership before execution.

All **5,040 observed heavy atoms** are excluded from optimization variables;
only the **556 diffusion-generated heavy atoms** and disposable auxiliary H are
mobile. OpenMM's `Modeller.addHydrogens` is H preparation, **not external MODELLER**.
Its output must retain every input heavy coordinate bit-for-bit, including the
generated coordinates. No PDBFixer heavy rebuild, OXT completion, water/ligand
reattachment, source repair, chirality reflection or anchor release is allowed.
The reduced-coordinate method enforces FF-equilibrium X–H lengths without OpenMM
constraints involving massless observed atoms.

The declared ff19SB/NoCutoff/no-solvent job has four phase caps of
**100 + 75 + 75 + 1,000** iterations, no restart/rechart and a **1,800-second total
limit including H preparation**. The existing serial host-pressure/RSS guard is
reused; its historical module name does not invoke the external MODELLER worker.
Existing chirality, amide, peptide-angle and geometry-conditioning forces are
reused unchanged. Solver termination alone is not acceptance: independent sparse
Cartesian tangent-force RMS must be at most 10 kJ/(mol·nm), alongside fixed/native
and serialized custody, X–H equalities, global L and unchanged geometry gates.

`scripts/audit_1non_generated_refinement.py` adds a read-only raw-to-prepared heavy
coordinate check and independently recomputes candidate geometry, reusing the
separate sparse-Jacobian and native/serialized custody auditor. The new tests cover
case-sensitive ownership, observation/context refusal, phase limits and detection
of a 1e-10 nm heavy-coordinate change. **86 focused host tests pass; Ruff passes.**

- [Frozen one-job refinement declaration](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/plan.json).

**Completed; refined candidate remains rejected.** One CPU job finished in
**222.8 seconds** (about 407.6 MB peak RSS), with no resource stop, repeat or new
diffusion. H preparation took 9.1 seconds and added 5,812 disposable H; the
independent audit verifies that all **5,596 input heavy atoms** remained
bit-identical through preparation. All four phases reached their declared
iteration caps (100/75/75/1,000); the zero process exit status means technical
completion, **not optimizer convergence or scientific acceptance**.

The independent sparse-Jacobian audit confirms all **5,040 observed atoms**
remain exact in retained native states and serialized output, all 556 generated
atoms alone own the `.dat`, and serialized global L-Cα passes. The largest final
X–H equality error is **9.44e-16 nm** (tolerance 1e-12 nm). Auxiliary H are not
included in the heavy-atom output. Final tangent-space DOF-normalized RMS force
is **25.745 kJ/(mol·nm)**, failing the unchanged **≤10** stationarity gate.

The declared refinement improves but does not solve the generated geometry:

| Full CLI finding | Raw | Refined |
|---|---:|---:|
| Severe pairwise clashes |357|37|
| Severe heavy coincidences |4|0|
| Error-level bond lengths |7|1|
| Error-level bond angles |0|1|

The mask assessment still reports an amide finding and a Ramachandran outlier.
The stereocentre guard additionally rejects near-planar generated Cα at **B90
(4.10°)** and **D76 (9.55°)**; these are geometry warnings, **not D assignments**.
Generated-heavy RMS/max movement from the raw diffusion candidate is
**4.231 / 13.333 Å**. This substantial repair motion is not evidence of native
accuracy or a good raw diffusion proposal.

Final junction C–N distances (Å):

| Chain | C72–N73 | C90–N91 |
|---|---:|---:|
|A|1.2474|1.4746|
|B|1.3225|1.4712|
|C|1.3341|1.3774|
|D|1.3503|1.3674|

The unchanged junction gate still rejects the A/B right junctions. The output
remains a **private rejected candidate**, retained beside the raw proposal and
all phase states. No anchor release, chirality reflection, solvent/context
reattachment, threshold relaxation, automatic retry or budget extension occurred.
The independent auditor made zero optimizer calls and recomputed geometry from
the serialized candidate. This result establishes custody through refinement,
not a successful restoration. Further work requires a separately declared scope;
do not silently continue the optimizer or substitute another proposal.

- [Worker result and phase termination records](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/execution/1NON-generated-only-seed-23/result.json).
- [Independent preparation and geometry receipt](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/independent-audit/refinement-receipt.json).
- [Independent sparse projection and native/serialized custody](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/independent-audit/receipt.json).
- [Independent final geometry assessment](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/independent-audit/independent-assessment.json).
- [Raw/refined findings and all junction distances](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/independent-audit/raw-refined-comparison.json).

### Read-only failure localization (2026-11-18)

The next user authorization was for analysis of **saved states only**, not an
iteration extension, another refinement or another diffusion call.
`scripts/analyze_1non_refinement.py` binds the retained inputs and checks native
fixation and exact reconstruction from saved reduced coordinates. It inventories
clash ownership, generated Cα geometry, junction distances, sparse-projected force
localization, callback trends and stereographic chart scaling. No OpenMM Context,
new force-field evaluation, optimization, sampling or coordinate repair is used.
The retained proposal, candidate, request and execution evidence remain unchanged.

**Clash localization:** all **37** remaining severe pairwise contacts are
**generated–fixed**, and every pair is within one explicit chain. There are no
remaining severe generated–generated, fixed–fixed or cross-chain contacts in the
final CLI report. These are pair counts, not independent defects. Generated B90
participates in 12 contacts, A73 in eight, D77 in eight and D78 in six. The worst
overlap is B64:C–B90:CB (**1.275 Å**). This does not implicate missing tetramer
conditioning or justify releasing the observed scaffold.

**When the distortions appear:** the raw B90/D76 Cα out-of-plane angles were
48.82°/51.47°, without near-planar generated Cα findings. At the ends of phases
1 and 2 these angles were about 54°/53° and 52°/53°; all eight junction C–N
distances were **1.3338–1.3373 Å** after phase 2. After the stronger soft-steric
phase 3, B90/D76 were **19.01°/37.65°**. Full-FF phase 4 further flattened them:
B90 is first seen below 10° at saved callback 50 (8.96°), and D76 at callback
100 (8.85°). The last preceding saved observations are callback 5 (B90 10.21°)
and callback 50 (D76 22.31°). Missing intermediate states were **not replayed**;
these are onset intervals, not exact first failing iterations. Final right
junctions A/B lengthened to 1.4748/1.4713 Å in native coordinates; the final
serialized distances are separately retained above.

**Residual force is not mostly auxiliary H on observed atoms:** final squared
projected-force fractions are **65.47% generated heavy**, **33.61% H on generated
heavy**, and **0.92% H on fixed heavy**. B90 alone contributes **34.89%** of the
total squared norm, D87 9.75%, D76 7.58%, D78 5.19% and A90 5.13%. These are
fractions of squared force, not force magnitudes or energy attribution. The
independent sparse projection agrees with the previous audit. Final H-chart
local contraction ranges from 0.0257 to a median 0.9905, with no H below 1e-3;
this provides no evidence for a near-antipode chart collapse as the principal
explanation at the returned states, but is not a general conditioning proof.

**The chirality guard was present, not accidentally disabled.** Read-only
inspection of the captured system shows the group-0 positive-volume wall at
`k=10000`, `qmin=0.2`. Final B90/D76 relative signed volumes are 0.1926/0.1807;
the wall is active, with analytically computed penalties only 0.274/1.869 kJ/mol.
It is an anti-inversion wall referenced to the raw positive volume, **not a hard
constraint on a healthy Cα out-of-plane angle**. Passing global L can coexist
with failing stereocentre shape. This is a verified objective/acceptance
distinction, not proof that changing this term alone will resolve the pocket.

**Budget versus basin remains unresolved.** None of phase 4's 1,000 accepted
callbacks passes the force threshold. During its last 100 callbacks, tangent RMS
ranges **15.47–69.77** (median 22.32), while energy still decreases by 206.8 kJ/mol.
The run stopped at its declared iteration cap, not line-search failure or a
demonstrated stationary minimum. Ongoing energy descent does not establish that
extra iterations would repair the persistent clashes, Cα shape or junctions.
Energies from different phases must not be compared as the objectives differ.
Without a controlled intervention or force-term analysis, a causal attribution
to budget alone, the FF alone, or an inescapable loop basin would be premature.

**Recommended next experiment, not executed/authorized by this analysis:** first
test a private generated-only Cα-shape protection compatible with immutable
observations and X–H elimination; then, if its numerical tests pass, separately
declare one matched-budget replay of the **same saved diffusion proposal** with
that single objective change. Keep the existing source, ownership, geometry
thresholds and phase caps; no new diffusion, anchor release or production-default
change. Compare all raw/intermediate/final states, not just the final clash count.
This is preferable to a blind iteration extension, but its success is not
inferred from the current negative result.

- [Source-bound saved-state analysis and full contact list](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/failure-analysis-2026-11-18/analysis.json).
- [Captured chirality-wall parameters and analytic endpoint penalties](../../tmp/1non-generated-only-refinement-2026-11-17/attempt01/failure-analysis-2026-11-18/captured-chirality-wall.json).

Verification: **88 focused host tests pass**, including new contact-ownership and
callback-summary tests; Ruff and agent-documentation validation pass. These checks
do not change the rejected scientific outcome or qualify a new refinement method.

### Matched generated Cα-shape control (2026-11-18): rejected

One authorized replay reused the exact prepared coordinates (including 5,812 H),
baseline system and four phase caps. A group-0 one-sided signed-elevation sine
penalty with a 20° target was added to 68 generated non-glycine centres. Three
OpenMM numerical tests checked energy and all Cartesian force components against
finite differences in inverted, flattened and healthy geometries; all seven focused
tests passed. There were no new diffusion calls. All 5,040 observations remained
exact natively and on serialization; generated-only `.dat` and global L passed.

The final angles improved from 4.10° to 17.72° (B90) and from 9.55° to 18.77°
(D76); neither triggers the 10° warning. At phase ends 1–3 the B90/D76 angles
were 53.84/53.30°, 52.08/53.05°, and 20.40/39.57°. Yet 37 severe generated–fixed
contacts remain, with **36 of the same atom pairs** as in the previous 37. A/B
right junction lengths still fail (1.4740/1.4609 Å). Final projected RMS force
fell from 25.745 to 17.574 kJ/(mol·nm), above the limit of 10. Generated-heavy
RMS movement from the common raw proposal increased from 4.231 to 4.661 Å.
The restraint resolves the final shape warnings in this case; it does not resolve
the placement of the loops against the scaffold. Further tuning of the angular
penalty alone is poorly supported by this comparison.

- [Frozen declaration](../../tmp/1non-ca-shape-control-2026-11-18/attempt01/plan.json)
- [Independent preparation and geometry checks](../../tmp/1non-ca-shape-control-2026-11-18/attempt01/independent-audit/refinement-receipt.json)
- [Independent phase forces and fixation](../../tmp/1non-ca-shape-control-2026-11-18/attempt01/independent-audit/receipt.json)
- [Raw, intermediate and final geometry and clash identities](../../tmp/1non-ca-shape-control-2026-11-18/attempt01/saved-states/analysis.json)

### Diffusion-time residue placement control: rejected

One full-tetramer seed-23 cc94 sample used the same fixed-scaffold schedule and
input. In the last 250 of 500 updates, each generated residue's clean prediction
was translated away from its closest eligible fixed atom within 3 Å, capped at
0.35 Å per update. Direct same-chain sequence neighbours were excluded. The
shift was propagated to self-conditioning. Readback checked the 250 before/after
arrays and the common translation within each residue; 35–47 residues moved per
guided step. Native preservation passed for all 5,040 observations before source
copying (maximum source-frame error 7.29e-6 Å). No refinement was performed.

Compared to the raw inpainting baseline, severe generated–fixed contacts fell
from **328 to 156**, but generated–generated contacts rose from **29 to 51**.
Six of eight junctions now fail (baseline four); bond-length errors doubled from
7 to 14 and bond-angle errors rose from 0 to 29. The four severe coincidences
remain in count. Translating residues independently improves scaffold clearance
in this case at the cost of backbone geometry. A next control should couple
neighbouring residue displacements and protect junction lengths during diffusion.

Attempt 01 failed pinned-source preflight because the unit-test invocation left
17 untracked bytecode files; it made zero sampler calls. The bytecode was moved
to local evidence, and the unchanged intervention was declared again. Attempt
02 completed the sole sampling call in 58.6 s. Use `PYTHONDONTWRITEBYTECODE=1`
for isolated tests against the pinned tree.

- [Retained preflight refusal](../../tmp/1non-placement-control/attempt01/native-sampling/trace.json)
- [Native-first receipt](../../tmp/1non-placement-control/attempt02/independent-readback/receipt.json)
- [Raw placement comparison including all eight junctions](../../tmp/1non-placement-control/attempt02/independent-readback/placement-comparison.json)

### Connector check at the actual cc94 network input

The first denoiser call was reproduced with seed 23 and the pinned checkpoint.
Observed coordinates reached the three crop channels of the network input after
division by 10.38: exactly 1,260 atoms per chain and 716 active sequence rows.
A single proper rigid fit from the unscaled crop to the deposited protein gave
7.29e-6 Å maximum error; all four chains shared the frame. Generated predictions
at this call matched the archived first step exactly. With noisy coordinates and
all other inputs unchanged, zeroing just chain B's crop changed generated-atom
RMS predictions on A/B/C/D by 4.57/33.12/7.28/4.43 Å, versus zero in an identical
repeat. This tests actual cross-chain sensitivity at high noise, not late-step
steric accuracy. The diagnostic stopped before the first ODE update. Independent
atom-name readback of the 556 native generated atoms into each inspection PDB
(02/03) agreed within 0.00082 Å, consistent with PDB precision.

- [Live network inputs and perturbation report](../../tmp/1NON-inspection/details/connector/report.json)
- [Independent output-frame check](../../tmp/1NON-inspection/details/connector/output-frame.json)

### Test-set work is data curation, not a modeling run

Keep these operations distinct: reading source metadata, defining a test manifest,
downloading a deposited file, and executing coordinate completion. A request to
inspect or organize the test set does not authorize the subsequent operations.
For the retained list, reuse local metadata and source files by default.

The inspected collection helpers (`scripts/collect_real_gap_candidates.py`,
`scripts/acquire_real_gap_sources.py` and
`scripts/collect_additional_long_gap_candidates.py`) use experimental method,
observed/missing residue counts, gap spans and sequence-cluster grouping, with
deterministic selection and a pinned historical accession list. Their inspected
selection logic does not rank biological activity, pathogenicity or binding
performance. This code inspection does not establish suitability of every
individual accession or explain a particular assistant refusal.

For further test-set curation, document structural coverage, duplicate/family
grouping, construct-sequence evidence, data quality and reproducible split rules.
Preserve accession titles, organism/function annotations and source provenance;
do not conceal identities or relabel biological objects to avoid review. Do not
expand curation into selection for enhanced biological properties. If a specific
entry or operation needs separate assessment, state the concrete issue rather
than treating all public experimental structures as either automatically safe
or categorically prohibited.

**Use the existing retained list.** Source recruitment and model-size calibration
remain paused. The old 200-entry ambition, length quotas and large job matrices
are historical planning targets, not current execution requirements.

The user has removed a separate case-admission review phase as a prerequisite to
continuing restoration on that list. The 177 pending review labels are provenance,
not a new mandatory queue. This does not remove per-task sequence/identity,
context, supported-chemistry, numerical, resource or output-validation checks, nor
does it reinstate excluded/deferred structures.

At handoff creation, previously declared calculations were settled; this document
did not itself queue another restoration job. The user subsequently authorized
autonomous implementation; the additional retained outcomes are recorded below.
Input preparation accounted for 233 source tasks:
194 ready, 31 refused and eight time-blocked. These are preparation outcomes,
not 194 accepted restorations or a completed independent evaluation cohort.
No final evaluation split is frozen.

### Latest retained example: 1KTZ

- The private four-chain candidate contains 468 residues and **710 completed
  heavy atoms** within its frozen mask. Those additions are modeled, not measured.
- Independent readback confirms 2,986 original protein observation records and
  326 original water records unchanged, global L-Cα geometry and the reviewed
  missing THR stereocentre. The separately confirmed native target has 21 S–S
  pairs; the original 18 explicit source records remain preserved.
- Consistent full diagnostics report no new geometry findings in the current
  export, eight unchanged source warnings and **two real missing OXT errors at
  receptor B136/D136**. Those OXT atoms are outside the frozen 710-atom mask;
  the errors are not baseline findings and must remain visible.
- Receptor terminal chemistry is deferred to future user-facing preparation,
  not a blocker for the original restoration task or continuation on the retained
  list. Do not ask again about the already confirmed TGF-β3 redox target.
- Do not silently add OXT, caps, hydrogens, charges or protonation choices.
  The candidate is not chemically complete or MD-ready. Its result does not
  establish missing-coordinate accuracy, broad backend qualification or topology/
  trajectory correctness.

Evidence: [final candidate and provenance](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/README.md),
[independent readback](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/readback-receipt.json),
[consistent validation reproduction](../../tmp/mac-real-gap-continuation-2026-11-14/validation-consistency-review-attempt01/report.json).
These local artifacts are not guaranteed to arrive through Git; report missing
evidence rather than recreating it or assuming a pass.

### Autonomous full-context continuation: four additional cases

Historical authorization: the user explicitly requested continued implementation without confirmation at
each routine step. Select bounded next tasks from the retained list, declare
their inputs/backend/seeds/budget, execute, independently audit, and continue to
other eligible tasks after a case-specific refusal. Do not stop merely because
an old handoff did not name the next accession. This authorization does not
remove the chemistry/identity/context/resource checks or permit silent scope
changes, new recruitment, broad job matrices or public capability changes.

Four additional single-protein full-context MODELLER controls executed on Mac:

| Entry | Retained result |
|---|---|
|1CC6|Native fixation passes for all 3,374 prepared-protein/FAD/PHB/water atoms. Raw and generated-only derivatives remain rejected: four angle warnings after the original CG continuation, five new geometry findings after the separate standard-stereo continuation.|
|9M3K|Native fixation passes for all 2,670 prepared-protein/water atoms. A serial-alias materialization failure was retained and recovered from saved native bytes without another model. The separate standard-stereo derivative retains one new boundary-angle warning and is rejected.|
|1X0X|One seed-104 native proposal passes the heavy-coordinate contract without extra refinement: 23 completed atoms, all 2,917 prepared-protein/NAD/sulfate/water atoms fixed, global L and zero new full-geometry findings.|
|5HJI|One seed-104 native proposal passes the same contract without extra refinement: 18 completed atoms, all 3,116 prepared-protein/adenosine/water atoms fixed, global L and zero new full-geometry findings. Its 24 previously completed fixed-conditioning atoms are separate from those 18 new gap atoms.|

The four cases are development-exposed, not an independent evaluation cohort.
Source preparation and fixed conditioning remain distinct from the new gap atoms;
the canonical `.dat` sidecars own only the declared gap-completion atoms.
The successful candidates retain 190/253 source warnings respectively, with no
CLI errors; neither result establishes coordinate accuracy or MD readiness.
No threshold was relaxed, observed anchor released, new source downloaded or
public capability enabled. Failed raw stages and failed derivatives remain failed.

The private adapter carries all retained context into the native model, selects
only the exact generated atoms and skips documented full-model `to_iupac()` atom
swaps. Native method-call evidence and direct pre-materialization coordinate
audits are retained. This is not a production MODELLER default change. An offline
`to_iupac()` control on the older 6Z38 initialization produced 94 changed named
atoms rather than its final 90; it does not establish the exact historical call
sequence or retroactively pass 6Z38. Source serial aliases are resolved without
changing coordinates or any other observed record field; topology-reload metadata
comes from the saved native proposal, never by normalizing experimental occupancy.

Evidence:
- [Four-case proposal declarations and retained stages, first pair](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt01/plan.json),
  [second pair](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt02/plan.json).
- [Independent negative-case readback](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-stereo-independent-readback-attempt01/receipt.json).
- [Independent successful-case readback](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-next-independent-readback-attempt01/receipt.json).
- [Private review exports and provenance](../../tmp/mac-real-gap-continuation-2026-11-14/restored-heavy-review-attempt01/index.json).
- [Consolidated outcomes and development-exposure overlay](../../tmp/mac-real-gap-continuation-2026-11-14/autonomous-full-context-continuation-summary.json).
- [Software verification](../../tmp/mac-real-gap-continuation-2026-11-14/full-context-modeller-controls-attempt02/software-verification-final.json):
  2,866 fast tests passed, three optional Torch/Biotite skips and 33 slow tests
  deselected; final focused tests, Ruff and agent-documentation checks pass.

Review exports are deliberately not public diffusion bundles or manual donor
approvals. Public CLI/API/Homology/GUI integration and broader qualification
remain separate work. Local evidence must be transferred explicitly with a
continuation; missing evidence is not a pass.

### Historical agreed stage: finite anchor-mobility calibration

At that stage, the task was **calibrating observed-anchor mobility and collecting
evidence from unsuccessful restorations**. Pipeline modifications, a required
reconstruction-success percentage and broader qualification are deferred until
those data are available. The first deliverable is a map of trade-offs, not a
preselected winning threshold or a production-policy change.

Start with retained proposals and controls, then assess dependencies on a separate
prospectively selected cohort from the existing inventory. The later cohort has
no approved membership, size or split and is not automatically a held-out
validation set. The old 200-entry/80–40–80 targets remain superseded.

#### Saved-proposal panel

Use 4–6 structures and exactly one retained base proposal per structure. The
provisional expanded panel is:

| Entry | Internal gaps, residues | Context and outstanding input checks |
|---|---|---|
|3C7D|10|NAD and pyruvate: verify actual force-field participation, not just retained output records.|
|6Z38|9|Context includes 3,115 ANISOU records. Bind exact proposal/native transport lineage; the historical 90 changed named observations remain a failure.|
|1NON|4 × 18|Four target chains, 358 waters; verify exact joint proposal and full context.|
|2BE1|2 × 8, 2 × 10, 2 × 20|Target chains A/B each have 339 residues and D has eight; preserve D and all 26 waters.|
|1HE8|Historical ten-gap stress case; exact internal/terminal breakdown not yet bound|Verify exact retained proposal, source mapping and whole-entry metal/ligand context. Earlier protein-only experiments are not whole-context qualification.|
|1KTZ|2 × 18|Four protein chains, 326 waters, already confirmed 21-pair disulfide target; bind exact retained 710-added-atom candidate and lineage. Do not re-ask the redox decision.|

This is **not an executable file manifest** or a guarantee of six admitted cases.
Do not substitute a protein-only projection, delete context or silently replace
a blocked structure. Record demonstrated blockers and coverage limitations;
membership changes require an explicit amendment. Bind original/prepared source,
context, request, base proposal, generator, seed and prior continuations before
execution. One proposal per case does not establish generator-independent effects.

Main scientific cases must have meaningful internal gaps. 1X0X/5HJI are short
N-terminal technical controls, not main successes; 1CC6/9M3K also lack internal
gaps. 6ZH0 has a one-residue internal gap, but its remaining A1–A2 defect is at
the N-terminal reconstruction junction, not evidence about internal-loop mobility.
Preserve these outcomes without inflating the main panel.

Each cell starts from the **same base proposal**, never another cell's optimized
derivative. All declared generated atoms, including terminal additions already
in the proposal, remain jointly mobile. Report internal and terminal junctions
separately. Only internal-gap flanks define the observed mobile windows; there
is no spatial-window arm initially.

#### Agreed grid and movement contract

| Factor | Agreed values |
|---|---|
|Observed windows|Zero motion; one observed flanking residue on each side of each internal gap; two consecutive observed residues outward on each side. Union overlaps once; never bridge another gap or chain boundary.|
|Anchor atoms|Observed heavy atoms in declared residues, uniformly weighted; no B-factor/ANISOU/occupancy weighting. All observed heavy atoms outside the window remain fixed.|
|Envelopes, Å|Per residue: (N/CA/C/O backbone RMSD, Cα displacement, maximum heavy-atom displacement) = (0.10, 0.10, 0.50), (0.25, 0.25, 0.75), (0.50, 0.50, 1.25), (0.75, 0.75, 1.50).|
|Positional strength|1,000 / 10,000 / 100,000 kJ/(mol·nm²), with energy `0.5*k*displacement²`. Historical common-refiner base: 10,000.|
|Envelope mechanism|Existing one-sided penalty walls, fixed coefficient 1,000,000 kJ/(mol·nm²), separate from positional strength; strict independent final movement audit.|
|Reference frame|Original source observations, no fitting or ratcheting between phases. Disclose earlier source preparation separately.|
|Count|Two mobile windows × four envelopes × three strengths + one zero-motion control = **25 per proposal**, 100–150 for 4–6 admitted cases.|
|Successful zero control|Still execute all 25 cells to measure unnecessary-motion costs; research only, not production escalation behavior.|

Envelopes participate in optimization, not just post-hoc classification. The
walls are soft penalties, **not mathematically hard constraints**. Independently
check per-residue movement both before and after PDB serialization. Exceeding a
limit is a failure: no clipping, rounding away excursions or widening thresholds.
These are experimental engineering controls, not physical uncertainty estimates.

Preserve occupancy and ANISOU. Where an anchor moves, disclose its unchanged
source tensor as retained metadata, not a newly refined experimental tensor.
No tensor transformation or fabrication is authorized.

#### Agreed physical and numerical setup

Use one common refiner across proposals; retain proposal-generator identity as
a separate factor instead of tuning each backend independently.

| Setting | Agreed setup |
|---|---|
|Engine|OpenMM CPU, one thread; record actual versions/properties, no silent fallback.|
|Base force field|`amber19/protein.ff19SB.xml` + `amber19/tip3p.xml`; bind XML resources and all additional context parameters. These files alone do not parameterize arbitrary ligands/metals.|
|Environment|Nonperiodic `NoCutoff`, no added solvent box or implicit solvent; preserve original waters and relevant heterogens.|
|Bond treatment|No bond-length constraints, `rigidWater=False`; preserve declared covalent/disulfide graph.|
|Temporary chemistry|One disclosed chemical model and identical initial temporary hydrogens per case across cells. Existing helper uses pH 7.0; record actual variants, charges and H inventory. Do not infer ligand microspecies or metal coordination from that heuristic.|
|Preparation|Up to 100 bonded-geometry + 75 first soft-steric + 75 second soft-steric iterations, following the versioned existing conditioning schedule.|
|Main optimization|Up to **1,000 full-force-field iterations**; up to 1,250 across four phases. No restarts or automatic extensions.|
|Tolerance|Existing 10 kJ/(mol·nm), identical across cells. Normal API return is not proof of convergence.|
|Other forces|Declare/freeze existing chirality, trans-amide, peptide-angle terms and soft-core schedule; no tuning per cell or backend.|
|Initialization|One declared preparation seed and identical pre-optimization temporary particle state per case, reused across all cells; minimization, not MD sampling.|
|Resources|Sequential execution, 30 minutes per cell including preparation/all phases; retain existing RSS 8 GiB, host-pressure and 10 GiB evidence-cap safeguards. Check available storage before launch.|

Temporary H atoms are auxiliary model particles, not observed coordinates or
publicly owned gap atoms. Record temporary heavy-atom completions and their
force-field role; never silently extend output ownership or call them observed.
In particular, 1KTZ B136/D136 OXT chemistry stays deferred and visible. Preserve
all source records outside authorized coordinate changes in materialization.

#### Execution gates (protocol; executed readiness updates below)

1. Bind exact files/digests for sources, contexts, requests, proposals, native
   transport audits, generated masks and anchor lists. Keep prepared conditioning
   distinct from original observations. Freeze a finite machine-readable manifest.
2. Correct the private window selector: historical `sequence2` in
   `scripts/run_1he8_traced_mobility.py` takes a ±2 axis neighborhood around each
   direct anchor, wider than the agreed two flanking residues. Do not reuse it
   unchanged for this arm.
3. Before research optimization, preflight actual full-context topology and
   parameters: ligands, ions, waters, connections, particle identities, charges,
   variants and fixed/mobile masses. Output record preservation alone does not
   establish force-field participation. Verify usable PDBFixer resources, not
   merely successful import.
4. Reuse valid existing chemistry where supported. A new ligand-charge, metal or
   chemical-state model needs a separately agreed specification; never choose it
   based on reconstruction success. Record a blocker without deleting context.
   Supported cases need not wait for every other case, but each must pass its own
   frozen preflight before its grid starts.
5. Test the private adapter: exact windows, identical initial temporary states,
   unchanged references, phase budgets, wall participation, all-generated joint
   mobility, native fixation, context preservation and independent serialized
   readback. Preserve production defaults.
6. Freeze executed-code and resource hashes and provide a ready/blocked ledger
   before launch. This document does not authorize running unverified inputs or
   unsettled chemistry. At protocol freeze, no study variant had executed.

Source inspection covered `scripts/run_1he8_traced_mobility.py`,
`src/dvbfixer/model/diffusion/boundary_refinement.py`,
`src/dvbfixer/model/diffusion/geometry_conditioning.py` and
`src/dvbfixer/model/diffusion/anchor_policy.py`. Existing 1KTZ native-context
evidence provides a route to inspect, not a pass for the new anchor adapter.
Source inspection did not establish whole-panel parameterization readiness.

#### Required reporting and first-stage closure

Retain for every cell:

- Exact input/code bindings, ownership and anchor identities; full temporary
  system/topology, initial particle state, chemical choices and force definitions.
- Actual phase calls, iteration observations where available, elapsed time,
  resource stops, energies, final forces and active force groups. Do not invent
  an API-unexposed termination reason. Different restraint Hamiltonians do not
  provide directly comparable energy-based quality scores.
- Native fixation before source copying; unchanged outside-window records;
  per-residue backbone RMSD, Cα and maximum heavy displacement before and after
  serialization, in the source frame.
- Whole-complex/junction geometry, sidechain stereochemistry, global L, declared
  links, Ramachandran/chi12, steric findings including relevant context/waters,
  and changes to source contacts/secondary structure. Bind analysis methods and
  settings before execution; do not invent acceptance thresholds from outcomes.
- Unchanged source findings, introduced defects, terminal-chemistry limitations,
  generated motion and observed-anchor motion as separate fields. Retain raw and
  final failures; no accepted-only reporting or public publication.
- Distinct outcome flags for geometry failure, envelope violation, numerical
  failure, resource exhaustion, input/chemistry blocker and implementation error.
  Multiple flags may coexist. A timeout or code bug is not a physical mobility
  failure.

Report paired effects relative to the zero control and same base proposal:
junction improvements versus anchor motion, lost/new contacts, secondary-structure
changes and introduced defects. Show all cells and per-case/cross-case trade-offs.
Twenty-five derivatives are not 25 independent structures; repeated chains are
not independent biological samples. Geometric acceptance is not coordinate
accuracy, proof of an unthreaded path or MD readiness.

**First-stage completion means a complete map of outcomes**, not a minimum success
rate: every declared variant has completed or has a documented specific blocker;
independent reports and failures are retained; numerical/resource stops are
classified; no unresolved implementation errors remain; and the paired summary
states coverage limitations. A blocked variant is not an executed variant.
Finding no suitable regime is a valid negative result, not authorization for
undeclared retries, wider limits or changed geometry thresholds.

Review the map with the user before selecting a policy, modifying the pipeline or
starting the separately declared prospective cohort. No final validation cohort
or required reconstruction-success fraction is agreed.

#### Executed preparation: 2026-11-14

The user authorized starting the agreed work. Private tooling is now in
`scripts/anchor_calibration_preflight.py`, with an independent saved-particle audit
in `scripts/audit_anchor_calibration_preflight.py` and focused regressions in
`tests/test_anchor_calibration_preflight.py`. Production code/defaults are unchanged.

The first declaration exposed historical prepared-source differences; it did not
execute any worker. The second declaration explicitly replays those already-bound
preparations, without performing a new source repair: 52 previously moved protein
atoms in 6Z38 and 31 in 2BE1, with original coordinates retained separately.
Per-atom original/prepared coordinates and fixed-conditioning additions are in
each case's `historical-preparation.json`. These changes are not anchor motion
produced by this study. Native observation-fixation failures in the selected
source-restored MODELLER derivatives remain disclosed and are not retroactively
passed.

Five guarded force-field preflight workers completed; **zero calibration variants
executed**. These workers may perform the existing temporary PDBFixer/hydrogen
preparation, but stop before the first research minimization call. Process success
is not force-field readiness or scientific acceptance.

| Entry | Historical attempt02 preflight outcome |
|---|---|
|6Z38|Temporary system constructed; independent saved-coordinate/mass/nonbonded readback passed. 6,443 particles, 3,082 fixed heavy atoms in the two-residue preflight window. One auxiliary OXT at A432 is present only in the temporary FF model, not output ownership.|
|1NON|Temporary system constructed; independent particle readback passed. 12,482 particles, 5,250 fixed heavy atoms in the two-residue preflight window. Very large finite starting energy remains visible; this is not a geometry/convergence pass.|
|3C7D|OpenMM template matching stopped at terminal HIS (temporary residue index 403), missing an external C connection. NAD/pyruvate participation has not yet been established; this is the first encountered blocker, not proof they are parameterized.|
|2BE1|Template matching stopped at terminal VAL (temporary index 709), missing an external C connection. No automatic cap or additional output atom was introduced.|
|1KTZ|Template matching stopped at terminal ASP (temporary index 294), missing an external C connection. Deferred B136/D136 chemistry remains unchanged. The new mobile-anchor model has not passed its 21-pair graph check.|
|1HE8|Input assembly retains four LINK records and GNP/Mg context. The reused unlinked-context merger refuses it; no FF worker executed. This is an adapter/context limitation, not a failed reconstruction or demonstrated physical incompatibility.|

At that stage, only **two temporary systems** passed particle readback. No case was
declared grid-execution ready: shared-initial-state reuse across all cells,
cell-specific auxiliary forces, full movement/geometry/contact reporting and the
independent grid adapter remain to be implemented and verified. Do not promote
these two preflight snapshots to completed 25-cell studies.

Before changing terminal chemistry, ligand charges or the metal model, agree an
explicit temporary-force-field specification separately from output completion.
Do not bypass the failures with omitted waters/ligands, unreported terminal atoms,
`ignoreExternalBonds`, relaxed geometry or hidden topology changes. Continue
supported implementation work without treating blocked cases as executed seeds.

**Subsequent explicit user decision:** in the temporary FF copy only, model the
uncapped target-chain C termini as free carboxyl termini under the agreed pH-7
setup and allow the auxiliary OXT atoms needed for that model. Record exact
residue identities/additions. This does not extend source/output/`.dat` ownership
or settle final user-facing terminal chemistry, including 1KTZ B136/D136. Existing
terminal C covalent bonds must fail closed rather than be removed to fit this
model. Ligand microspecies/charges and GNP/Mg remain separate unresolved models.

Preparatory attempt03 repeated the same preflight with an immutable `protocol.md`
copy and pre-bound XML resources, preserving attempt02 and its readback. This
avoids making future verification depend on a living handoff file that is updated
after execution. No research variant was executed in either attempt.

**Attempt04, after the explicit temporary-terminal decision:** independent
particle and base-force-field readback now passes for **6Z38, 1NON, 2BE1 and
1KTZ**. 2BE1 has one auxiliary OXT at D26; 1KTZ has auxiliary OXT at B136/D136
and its actual bond graph matches the exact approved 21 S–S pairs. None of those
auxiliary atoms extends output ownership. 3C7D now reaches the distinct missing
NAD-template blocker; its NAD/pyruvate model is not established. 1HE8 remains
blocked at the linked GNP/Mg adapter/model gate.

`scripts/anchor_calibration_grid.py` implements cell-specific force/mass setup
from the captured **unfrozen base System plus one initial NPZ per case**, without
rerunning hydrogen placement for each cell. It retains all generated atoms jointly
mobile, writes phase traces, compares native and serialized anchor movements,
checks native fixation and outside-window physical records, records full
diagnostics and emits generated-only `DatRecord` sidecars. Implementation errors
stop execution rather than masquerading as scientific refusals.

Descriptive environment metrics are declared before optimization: 4.5 Å observed
heavy-atom residue proximity, excluding covalently connected residue pairs, and
MDAnalysis C3 assignments on continuous peptide segments of at least six residues,
using auxiliary FF backbone H (missing-H donors disabled). The latter is an
**intrasegment comparison**, not full interchain-sheet DSSP qualification; its
reference is the same initial completed proposal, not experimental missing atoms.
Neither metric introduces an acceptance threshold or selects a policy.

The serial executor now accepts an optional **shorter** per-job timeout. The grid
charges the measured shared-preparation time against each variant's 1,800-second
ceiling, so reusing preparation does not silently grant extra optimization time.
Historical/default jobs still have 1,800 seconds; extensions are refused.

Four supported cases permit a finite **100-cell** first series, with 3C7D/1HE8
retained in the blocked ledger, not replaced or counted as executed. Final
software verification and the frozen execution declaration precede launch; no
research variant had executed at that preparation checkpoint.

**Execution launched 2026-11-15:**
`tmp/anchor-calibration-2026-11-14/grid-attempt01/plan.json` freezes 25 cells each
for 6Z38, 1NON, 2BE1 and 1KTZ (100 total), in that order. The background serial
executor uses the existing host/RSS/disk safeguards and shortened per-cell time
ceilings accounting for shared preparation. No alternate proposal, resampling,
restart, threshold change or production-default change is authorized.

Software gate: the first fast-suite invocation had 2,896 passes, three optional
skips and six missing-CLI `PATH` failures. With the environment bin directory
added only for the child invocation, the final complete fast suite passed
**2,903 tests, three optional skips, zero failures/errors**. The initial failed
report remains retained, not relabeled. Final frozen-adapter verification also
reconstructed all four initial states, checked finite forces and agreement with
the original common-refiner force states without executing minimization.

The live run retains per-phase states/reporters, exact native/serialized mobility
and `.dat` readback, full diagnostics and descriptive environmental metrics.
`summary.json` is written on exit, including partial/resource/code-error exits;
unlaunched cells remain unlaunched. Implementation/invariant errors stop the
executor and require an explicit retained correction, not a hidden retry. Do not
claim the 100 cells completed from the existence of the plan or launch receipt.

**Execution completed (before independent reanalysis):** all 100 guarded jobs
returned with retained reports and no resource stop. Ninety-four reached the
final reporting stage; six stopped at a global Cα-chirality gate (one 1NON and
five 2BE1), with rejected phase states retained. No repair, repeat optimization
or additional sampling was launched.

Post-run inspection found an **evaluation-scope implementation error**: the
protein-only `assessment_adapter(request, source)` was given full-context PDBs.
This produced `outside-generated-identity-mismatch` for waters,
`pdb-metadata-consistency` for the full ANISOU-bearing materialization, and
`retained-explicit-link-unexpected` for the separately approved additional 1KTZ
S–S annotations. Do not interpret the original all-failed mask flags as 0/100
scientific restorations. Do not erase those original reports either.

`scripts/audit_anchor_calibration_grid.py` now performs a separately declared,
read-only audit with **zero optimization calls**:

- Bind retained outputs and preserve the original reports/coordinates.
- Independently audit the full inventory, native fixed coordinates, physical
  outside-window records, ANISOU/occupancy/annotations, actual serialized mobile
  coordinates, generated-only `.dat`, global L and the exact approved 21-pair
  1KTZ link target. Independently reproduce rejected native D-Cα states.
- Build a protein-request evaluation view with every fixed/generated/anchor
  coordinate equal to the actual final full PDB. Source-record restoration,
  fitting and hidden anchor-motion erasure are forbidden. The frozen request is
  unchanged; context is assessed separately rather than silently omitted.
- Recompute protein-mask geometry and native/serialized movement under unchanged
  limits, distinguishing internal-gap and terminal findings.
- Run complete CLI diagnose on the prepared baseline and each retained final PDB,
  with waters and the unchanged 0.6/0.9 Å steric thresholds. The worker's earlier
  chemistry/Ramachandran/chi12/steric subset is not a substitute for full CLI
  postflight; in particular, retain deferred 1KTZ B136/D136 OXT errors explicitly.

The original grid runner must not be reused for another cohort without fixing
its assessment scope. Its private current implementation is now corrected to
use a coordinate-identical protein evaluation view and explicitly calls the
chemistry/Ramachandran/chi12/sterics result a subset, not full CLI postflight.
Future declarations use revision `anchor-calibration-grid-v2-scope-corrected`.
**That revision has not run new optimization.** The executed v1 code archive,
inputs, all original failures and result JSON remain unchanged. No production
pipeline/default was modified. Focused post-readback tests: 48 passed; Ruff and
agent documentation checks pass. The earlier complete 2,903-pass suite is a
separate pre-execution result, not a new full-suite claim for the readback edits.

#### Completed first saved-panel map

Independent full-context/request-view readback completed without resampling or
optimization. All 94 returned final candidates preserve native fixed observations,
outside-window records/ANISOU/occupancy and generated-only ownership. The six
rejected phase states independently reproduce D-Cα; a complementary audit also
confirms their native outside-window observations remained fixed. No rejected D
state was reflected, exported as accepted or relabeled a pass.

| Entry | Declared/reported cells | Returned final candidates | Protein-mask geometry + motion passes | Independently reproduced D-Cα refusals |
|---|---:|---:|---:|---:|
|6Z38|25|25|0|0|
|1NON|25|24|0|1|
|2BE1|25|20|0|5|
|1KTZ|25|25|18|0|

Interpretation boundaries:

- **6Z38:** all cells have zero local findings touching internal generated
  residues; the local failures are dominated by ten N-terminal clash findings.
  Do not present this as 25 failures to repair the internal nine-residue junction.
  Full CLI also retains a terminal OXT error and other introduced findings.
- **1NON/2BE1:** returned candidates retain real internal-gap geometry/steric
  defects. Local internal-residue finding counts range 24–40 / 15–20 across the
  returned cells; none clears the restoration mask and movement gates. One saved
  proposal per structure does not establish that every alternative proposal or
  backend would fail.
- **1KTZ:** the zero-motion control already passes the declared protein mask.
  Seventeen of 24 mobile cells also pass mask+motion; seven exceed the strict
  movement envelopes. Observed proximity-contact loss increases by up to two
  relative to zero, with no measured intrasegment C3 changes. This is evidence
  about unnecessary-release cost on this case, not a universal anchor threshold.
- **No cell has zero introduced full-CLI findings.** For 1KTZ the two deferred
  B136/D136 OXT errors remain, plus two CLI `chain_break` warnings explicitly
  describing ASP136 → HOH201 boundaries. Those protein-to-water warnings are not
  evidence of unclosed peptide junctions; retain the raw diagnostic and track the
  diagnostic limitation rather than silently removing warnings or changing
  production code. Mask success is not whole-entry chemical completeness.

The finite first-stage deliverable is complete: **a 100-row paired trade-off map,
failure evidence, two blocked cases and explicit coverage/metric limitations**.
It is not a selected/validated optimal policy, successful recovery of all panel
members, coordinate-accuracy validation or MD qualification. Twenty-five cells
from one retained proposal are not independent biological samples. No pipeline
change, extra seed, wider envelope, new source or prospective cohort has been
launched on the basis of these results.

**Next decision belongs with the user:** review this map and choose the separately
declared prospective cohort and/or bounded pipeline experiments for the retained
internal-gap failures. Freeze that next scope and analysis contract before new
execution; do not revive the old large cohort or replace blocked cases silently.

#### Failure investigation: 2026-11-17, read-only findings

The user requested a refusal investigation, **not new reconstructions or loosened
thresholds**. `scripts/analyze_anchor_failures.py` scanned all 100 saved endpoints
and 11 witness histories (initial state plus four recorded phases). It creates no
OpenMM Context, performs no minimization, repairs no coordinates, and preserves
all original reports. Four focused regression tests pass; production code and
thresholds remain unchanged.

| Case | Concrete saved-state evidence | Interpretation |
|---|---|---|
|6Z38|The initial proposal has seven severe generated-heavy contacts. The fixed-anchor endpoint retains ten; generated-heavy RMS force is 13.34 kJ/(mol·nm), while the massless-zeroed all-component diagnostic is 4.48. THR A1 Cβ orientation is +23.90° initially and +58.85° after phase 3, then opposite canonical R after the full-FF phase (−22.35°).|A bad terminal proposal remains in a geometrically bad endpoint, and the refiner introduces a sidechain stereocentre defect. The Cα wall does not protect THR Cβ. This is not an internal-anchor-window failure.|
|1NON|Initial actual-topology control finds 82 severe generated-heavy contacts; the closest is generated THR C76:CG2 to observed HOH C254:O at 0.269 Å. Fifteen of 24 returned final candidates have masked returned-force RMS above the configured scale of 10. Several show OH HG1 collapsing toward the same THR's CG2, including a 1–4 electrostatic exception with zero LJ epsilon. Some other endpoints meet the global force scale but retain real loop defects.|Mixed failure: deeply colliding initial proposals, numerical polar-H pathologies, and residual bad geometry at the aggregate force criterion. Do not treat every reported endpoint as a meaningful local minimum or infer that more anchor motion is the missing ingredient.|
|2BE1|All 25 last full-FF states have very large mobile residual forces; all 20 reported final candidates are above the configured force scale. SER B161:HG is the dominant mobile-force site in every cell and approaches HOH B12:O to 0.124 Å or less (0.002877 Å in `w1-e3-k10000`). HG has charge +0.4275e and zero LJ epsilon; the water oxygen has −0.834e and no exception for this pair.|Unphysical polar-H collapse in the selected unconstrained full-FF setup, not evidence that the internal loops are inherently unreconstructable. A normal minimizer API return did not establish convergence.|
|1KTZ|No newly exposed severe final contacts in the topology-only control; returned-force diagnostic is below 10 in every cell, as in 6Z38. The existing 18 mask+motion passes remain unchanged by this particular control.|Useful zero-motion positive control, not a universal policy calibration or chemical-completeness/MD pass.|

The force diagnostic zeros massless-particle components and reports their RMS
over the full component axis, separately from active-particle and generated-heavy
RMS. Comparing it with the configured force scale is not an exposed OpenMM
termination code or proof of convergence; a future numerical guard must match the
engine's actual normalization. The enormous 2BE1 gradients are not marginal
rounding/tolerance excursions.

**Concrete fixed-context contradiction in 2BE1.** SER B161:OG and HOH B12:O
are 1.658382 Å apart, and both are outside the declared mobile windows. Direct
Gemmi readback of the retained deposited mmCIF confirms those exact coordinate
vectors in model 1, occupancy 1.00, no alternate location; neither the historical
31-atom preparation nor this grid created the contact. This establishes the
contact's provenance, not a new chemical interpretation of it. In the declared
SER + TIP3P-water model there is no O–O covalent bond and that contact cannot be
removed while both observed heavy atoms remain fixed. Any source-context change
needs its own explicit preparation/chemistry policy; do not delete or move the
water silently.

**Confirmed diagnostic exclusion bug.**
At the executed study revision, `diagnose.steric._build_neighbor_graph` augmented its covalent graph with
every sufficiently short distance (1.90 Å for ordinary heavy-atom pairs), even
between canonical protein and water. It therefore invents the SER–water O–O edge,
then `_expand_exclusions` hides that pair and related 1–3/1–4 contacts. This is not
a valid H-bond-envelope exemption: 1.658 Å is below that envelope, and the actual
FF graph/annotations do not contain the bond. The existing 0.9 Å error overlap
would flag 1.382 Å overlap if the fake exclusion were absent.

A private **same-threshold, actual-FF-topology control** exposed at least one
hidden severe fixed contact in all 25 2BE1 endpoints and additional hidden severe
generated contacts in six states. It exposed no hidden severe final contacts in
the other three cases. In raw witness proposals it exposes 59 previously hidden
contacts in 1NON and three in 6Z38. This control retains actual 1–2/1–3/1–4,
same-residue, water-inclusion and heavy-H-bond rules; it does not tighten the
0.6/0.9 Å thresholds or silently replace the historical reports. Production clash
logic has **not** been changed by this investigation.

**Proposal/context mismatch, independently bound.** The retained seed-101
MODELLER template files for 1NON and 2BE1 contain no HETATM/water records and
match their historical declared source digests. They were not conditioned on the
explicit crystallographic waters later attached for common refinement. The very
large initial generated–water overlaps are therefore not evidence that the
generator already solved the full-context problem. This is a property of these
retained proposals/routes, not a claim about all MODELLER capabilities.

**D-Cα localization.** All six recorded D-Cα failures first occur in the full-FF
phase: LYS B78 in 1NON and VAL B261 in five 2BE1 cells. The initial and preceding
recorded phase states are L. All six centres are present in the active Cα wall;
the wall is finite, not an omitted force or a hard feasibility constraint. The
H collapse and inversions coexist in that phase, but endpoint inspection alone
does not establish that one caused the other. The post-phase rejection correctly
kept all six states rejected; no reflection or repair has been performed.

**What this changes about the earlier interpretation:** the original 18/100
count is retained, not inflated by dropping failures. But 2BE1 and many 1NON
cells are confounded by numerical/model failures, while 6Z38 is a terminal
counterexample and 1KTZ already passes at zero motion. These data cannot select
a physically justified universal anchor threshold. In particular, merely raising
the iteration budget or widening envelopes is not a supported remedy for an
unphysical near-singular H attraction.

Proposed next work **at the failure-investigation checkpoint**: first correct topology-based clash
exclusions and reject/report incompatible immutable source contacts; then add
explicit returned-state numerical/H-geometry guards and test a separately
declared H-safe refiner setup (e.g. FF-equilibrium constraints on covalent X–H
bonds, not fixed noncovalent H-bond distances or frozen raw sampled heavy-bond
lengths). Retest a small predeclared saved-proposal control
set before another 100-cell grid or a new cohort. Proposal regeneration with the
actual water context and THR sidechain stereochemistry protection are separate
changes to evaluate, not automatic post-hoc repairs of the frozen study.

#### Production diagnostic fix and authorized v2 revalidation

Main now includes `219ba8b` through PR #27 (`9cdd09b`). Python steric engine
version 2 never creates bonds from proximity. It uses supplied topology bonds
and unique intra-residue template names, preserving true 1–2/1–3/1–4 exclusions;
unknown/ambiguous intra-residue connectivity is a warning and partial status,
not an invented bond or a clean pass. Supplied parser bonds remain trusted;
this is not independent qualification of arbitrary ligand connectivity.
Successful Probe reports follow their own unchanged engine route.

The user authorized the next **read-only** step: rerun diagnostics on the retained
original/prepared sources, proposals and outputs, retaining prior reports and
the original 100-cell denominator. No optimization, proposal regeneration,
automatic repair or change of thresholds is included. Research/main integration
is `17b46b4`; README keeps the research scope link and main's release 0.10.5,
and the license declaration inventory was regenerated after the merge.

`scripts/reaudit_anchor_diagnostics.py` declares code/input bindings before
execution, reuses the independent inventory/native-fixation/movement/mask audit,
and adds explicit Python-v2 provenance independent of Probe availability.
Serialized heavy-atom candidates and rejected native FF states (with auxiliary
H/OXT) remain separate scopes. Original reports are hashed before/after; no
historical count is overwritten. Saved residual-force diagnostics are attached
without creating a new Context or claiming an OpenMM termination reason.

The historical `scripts/analyze_anchor_failures.py` proximity-control comparison
now pins its old 1.30/1.90/2.25 Å v1 inference constants locally because production
removed that helper. This preserves the meaning of that historical control; it
does not reintroduce proximity inference into production. Executed archived code
and previous evidence remain unchanged. Focused pre-execution checks: 42 tests
passed (v2 revalidation helpers, historical analysis, independent grid audit and
production steric exclusions).

Execution declaration:
`tmp/anchor-steric-v2-readback-2026-11-17/attempt01/declaration.json`.
**Revalidation completed:** all 100 cells are accounted for: 94 serialized
candidates independently re-audited plus all six rejected native FF states.
Original/prepared contexts and the four proposals were checked separately with
explicit Python-v2 provenance and unchanged water/0.6/0.9 Å policies. No source
repair, H placement, Context construction or optimization was performed. The
final hash check confirms historical reports and bound analysis code unchanged.

| Case | v2 mask + motion passes / declared cells | Additional v2 findings versus historical final reports | Remaining interpretation |
|---|---:|---|---|
|6Z38|0/25|None in the final PDBs.|All 25 retain severe generated-atom clashes, dominated by terminal reconstruction defects. This is not a demonstrated internal-junction failure.|
|1NON|0/25|One additional clash WARNING across the 24 final candidates; no additional severe final clashes.|All 24 retain severe generated-atom clashes; 15/24 have the previously recorded high force diagnostic. One rejected native D-Cα state remains rejected.|
|2BE1|0/25|2–11 additional clash findings per final PDB (34 ERROR and 39 WARNING findings across 20 candidates).|All 20 finals expose immutable severe clashes, including the deposited SER–water conflict; five finals also expose previously hidden severe generated contacts. All 20 retain severe generated contacts and high saved force diagnostics. Five rejected native D-Cα states remain rejected.|
|1KTZ|18/25|None in the final PDBs.|No severe final generated or immutable clashes; all 25 have the saved force diagnostic below the configured scale. The seven movement-envelope failures remain failures.|

All 94 serialized candidates and six native FF checks have complete connectivity
under the supplied/template policy (no `partial` outcomes). This is a resolved
check scope, not independent chemical verification of every supplied bond.
The result stays **18/100**, all from one 1KTZ proposal, including the minimized
zero-observed-motion control. It is not 18 independent reconstructions or a new
minimization result. No candidate is newly accepted by dropping a failure.
Full 1KTZ reports still retain the deferred B136/D136 OXT errors and two legacy
protein-to-water boundary warnings; those are separate from actual internal
peptide-junction closure and restoration-mask geometry.

Original-source and prepared-baseline findings differ: the 6Z38 deposited/context
view has three severe clashes and its historically prepared baseline none;
2BE1 has four and one respectively. These are historical conditioning effects,
not repairs performed during this readback. The retained 1NON proposal has 82
severe clashes versus none in its prepared observed baseline, independently
confirming the proposal/context problem with the production-v2 engine.

Evidence: [v2 summary and source reports](../../tmp/anchor-steric-v2-readback-2026-11-17/attempt01/summary.json),
[independent 100-cell receipt](../../tmp/anchor-steric-v2-readback-2026-11-17/attempt01/readback/receipt.json).
Each cell directory retains its full CLI report (where a final PDB exists),
deterministic v2 report and `v2-comparison.json` with generated/observed-anchor/
immutable classification and saved force diagnostics. Native H/OXT states are
explicitly labelled and not pooled with heavy-PDB counts.

Post-merge software checks: 42 focused revalidation/audit/exclusion tests and
80 additional diagnostic/preflight/grid tests passed; the latter retain one
NumPy binary-size runtime warning. Ruff, agent-documentation checks and the
regenerated license-inventory check pass. These are focused checks, not a new
full-suite result. Next work remains a separately declared numerical/H-safe
refiner control and immutable-context admission policy, not a wider grid or
silently modified historical source.

**User-confirmed priority:** deferred OXT additions are lower priority than
correct geometry of completed regions and their junctions. Keep OXT findings in
full-entry reports, but report them separately from restoration-mask success;
do not extend the generated mask or spend this step changing terminal chemistry.

#### Private guard implementation and saved-state controls

The user's next instruction authorized implementation. The first bounded piece
is `scripts/anchor_refinement_guards.py`, not a production-default change:

- **Immutable-context admission:** use current Python-v2 supplied/template
  connectivity with retained waters and unchanged 0.6/0.9 Å thresholds. Refuse
  unresolved connectivity or severe heavy-atom clashes whose two particles are
  both fixed. Report exact atom identities rather than moving/deleting a source
  water or reinterpreting chemistry. This does not assert that all admitted
  proposals are good; generated clashes and ordinary geometry still need repair.
- **Returned-state guard:** finite positions/forces, byte-exact numeric fixation
  of massless particles, independent active-Cartesian component force RMS <=10,
  X–H compression/extension screening, H coincidence checks across residues as
  well as within them, global L-Cα, and generated-residue stereocentre checks.
  Generated-heavy RMS and the historical all-axis diagnostic remain separately
  reported; no new generated-only force threshold is invented. A guard pass is
  not a geometry/accuracy pass, exposed engine termination reason, energy audit
  or proof of convergence. Full geometry/movement gates remain additional.
- **Covalent X–H targets:** extract each unique H-parent pair and its positive
  equilibrium length from the captured `HarmonicBondForce`, not sampled lengths
  or noncovalent H-bond distances. Missing/ambiguous FF terms fail rather than
  guess. The geometry screen uses the existing 20% bond-warning fraction against
  those equilibria; H collapse uses the existing 0.05 nm coincidence floor.
  These are detection gates, not newly installed constraints or repair forces.

**Confirmed refiner blocker:** OpenMM rejects a constraint involving a massless
particle (`A constraint cannot involve a massless particle`). Here the observed
heavy parents are massless specifically to enforce exact coordinate fixation.
A two-particle CPU regression reproduces the actual error; this is not a
hypothetical concern. The inspected 1KTZ control alone has 3,490 such X–H bonds.
Do not silently make observations massive/mobile, replace exact fixation with
strong restraints, or describe a finite X–H penalty as a hard constraint.
An alternative constrained optimizer/coordinate parameterization preserving
exact observations needs its own bounded design and numerical validation.
The guard v1 explicitly refuses already constrained systems because plain
Cartesian forces are not a valid constrained/projected stationarity test.

`scripts/check_saved_anchor_guards.py` froze six retrospective witness controls
before execution, reused the serial executor's RSS/host-pressure/storage guards,
and completed each within a 180-second limit. These controls only read saved
initial/full-FF states; there were **zero Context constructions, minimization
calls, H-placement calls or source repairs** in the scientific workers. The
separate two-particle regression does construct the intentionally invalid test
Context. Original systems are verified unchanged in memory, and bound code/input
hashes are rechecked at closure.

| Saved control | Immutable-context admission | Returned-state guard |
|---|---|---|
|1KTZ fixed|Admitted.|Pass; active force RMS 7.032 kJ/(mol·nm). Still subject to the existing geometry/movement gates and deferred OXT reporting.|
|6Z38 fixed|Admitted.|Reject generated stereocentre geometry; active force RMS 6.240 does not hide the THR-Cβ defect.|
|1NON fixed|Admitted.|Reject active RMS 12.070, covalent X–H geometry and generated stereocentre geometry. Its historical massless-zeroed all-axis RMS was only 9.093.|
|1NON w1-e1-k10000|Admitted.|Reject large active force RMS, X–H geometry, H collapse and generated stereo.|
|2BE1 fixed|Refuse incompatible immutable SER–water context.|Retained state also fails force/H/D-Cα/stereo checks; it was not reoptimized past the new refusal.|
|2BE1 w1-e3-k10000|Same immutable-context refusal.|Retained state also fails force/X–H/collapse checks; no new scientific run.|

These controls establish detection of the selected failures, **not repaired
reconstructions or a working constrained refiner**. Historical 18/100 results
remain historical; the guard suite does not replace that denominator or qualify
a new completion policy. The initial severe fixed-context refusal for 2BE1 is
separate from its numerical failures and does not authorize removing the water.

Evidence: [frozen six-control declaration](../../tmp/anchor-refinement-guards-2026-11-17/attempt01/plan.json),
[completed control summary](../../tmp/anchor-refinement-guards-2026-11-17/attempt01/summary.json).
Per-control `result.json` retains admission findings, native-state failure detail
and the exact FF X–H inventory, with archived executed code in `executed-code/`.
Focused checks: **54 tests passed**, including 12 new guard/constraint/stereo
tests; Ruff passed. Production defaults and public capabilities are unchanged.

**Next implementation gate:** choose and test an H-constrained numerical route
that preserves exact fixed observations (with correct projected-force checking)
before a new saved-proposal minimization control. Until then, do not start another
100-cell grid or silently substitute soft X–H penalties. The separate policy for
incompatible immutable source contacts remains refusal/reporting, not automatic
chemical/source edits.

#### Reduced-coordinate X–H route: implementation and bounded controls

Following the user's instruction to proceed, the private route now lives in
`scripts/anchor_constrained_coordinates.py` and
`scripts/anchor_constrained_control.py`. This supersedes the earlier statement
that no constrained numerical route was implemented; it does **not** change
OpenMM's massless-constraint limitation or claim scientific qualification.

- Fixed observations are excluded from optimization variables. Mobile heavy
  atoms remain Cartesian and jointly mobile; no whole-gap atom is frozen to make
  an example easy. Each auxiliary mobile H is represented by two stereographic
  direction coordinates about its actual heavy parent, at the FF covalent
  equilibrium length. H follows a mobile parent's translation analytically.
- OpenMM evaluates the original captured bonded/nonbonded and protective forces
  on assembled Cartesian coordinates; no constraint is inserted into its System,
  no new LJ/penalty force is introduced, and fixed coordinates are not replaced
  by strong springs. The analytic energy-gradient pullback includes H forces on
  mobile parents. SciPy L-BFGS-B optimizes the reduced coordinates.
- The initial temporary H vectors are radially projected to FF equilibrium while
  keeping their original directions. Both original and projected temporary states
  and the exact H inventory are retained. All heavy observations and context are
  unchanged. Inputs with observed/non-auxiliary H or undefined parent directions
  are refused rather than rewritten; no new hydrogen-placement call occurs.
- Every FF evaluation checks fixed-coordinate equality and X–H length residuals
  against a 1e-12 nm floating-point equality tolerance. This is coordinate
  elimination to floating-point precision, not a tunable soft wall. It does not
  constrain noncovalent H-bond distances or guarantee absence of H collisions;
  the existing returned-state collapse/geometry checks remain mandatory.
- Convergence is not inferred from SciPy success or small chart gradients.
  Independently project Cartesian forces onto the joint X–H tangent space after
  removing fixed DOFs, solving the per-parent `J Jᵀ` blocks. The reported RMS is
  `||F_projected|| / sqrt(3*N_mobile - N_XH)`; the unchanged scale is 10
  kJ/(mol·nm). Radial reactions and fixed-source reaction forces do not masquerade
  as mobile gradients. The v1 geometry guard receives explicitly labelled
  projected-force input; this DOF-normalized audit remains separately required.
- Each stereographic chart excludes the direction opposite its initial H vector.
  This is a local-optimizer limitation. There is no automatic rechart/restart;
  a small ill-conditioned chart gradient near that direction cannot pass without
  the independent Cartesian tangent test. No global-minimum claim is made.

Small tests verify fixed coordinates, FF length equalities at arbitrary variables,
analytic gradients by finite differences, the block projection against an
independent dense Jacobian, radial reactions, an ill-conditioned near-antipode
counterexample, and a real CPU OpenMM/reduced-coordinate synthetic minimization.
Seven new coordinate tests plus the existing guard/readback/exclusion checks
passed (61 tests total) before scientific control launch; Ruff passed.

**Frozen control declaration:**
`tmp/anchor-constrained-control-2026-11-17/attempt01/plan.json`.
Three retained zero-observed-motion cases are declared: 1KTZ positive control,
1NON numerical/geometry control, and 2BE1 expected immutable-context refusal.
Each uses its exact saved initial FF state, full context and generated mask.
CPU/one thread, NoCutoff, no added solvent, serial existing RSS/host/storage
guards, no retries, with historical preparation time charged against the original
30-minute ceiling. Phase caps are **100+75+75+250**, deliberately shorter than
the historical 1000-iteration full-FF phase. These are bounded implementation
smokes, **not matched-budget anchor-policy comparisons or a new cohort**.
2BE1 must refuse before Context creation; the raw source contact is not repaired.

`scripts/audit_anchor_constrained_control.py` is the separate read-only auditor:
it constructs a sparse global constraint Jacobian and projects saved raw forces
via `Jᵀ (J Jᵀ)⁻¹ J`, independent of the optimizer's per-parent projection code.
It checks retained callback/endpoint equalities, native and serialized fixed
coordinates, native-to-PDB generated coordinates, complete records/ANISOU,
generated-only `.dat` and serialized global L. It performs no optimization.
**Execution and independent audit completed.** All three jobs returned, with
two actual constrained-coordinate optimizations and one expected pre-Context
admission refusal. Both actual optimizations used all 100+75+75+250 iterations;
SciPy reported the iteration limit, **not solver convergence**. The independent
physical tangent-force test, geometry and custody audits establish the following
bounded outcomes:

| Control | Final result | Independent physical checks |
|---|---|---|
|1KTZ fixed|Passes restoration-mask geometry and returned-state checks; no new public acceptance.|Projected DOF-normalized force RMS **8.0065** kJ/(mol·nm), below 10. No X–H geometry failure, H collapse or generated-stereo finding. Fixed observations and all outside-mask records unchanged; generated-only `.dat` and serialized global L verified.|
|1NON fixed|Rejected, with retained candidate and full reports.|Projected force RMS **52.0365**, above 10; remaining severe heavy clashes, junction/angle/amide defects and four generated-stereocentre findings. No final H collapse or X–H geometry failures; absence of those pathologies does not repair or accept this proposal.|
|2BE1 fixed|Refused for incompatible immutable context before Context creation or optimization.|Original fixed SER–water conflict retained and reported; no chemistry/source edit or numerical bypass.|

The independent sparse-Jacobian projection agrees with the optimizer's per-parent
block calculation in all eight returned phase states. Equality/custody checks also
cover retained callback snapshots. Final X–H residual maxima are **1.33e-15 nm**
(1KTZ) and **8.88e-16 nm** (1NON), versus the declared 1e-12 nm numerical audit
tolerance. Native fixed coordinates remain numerically identical, and all
unselected serialized records, occupancy, waters, links and ANISOU remain exact.
The largest initial *temporary-H-only* radial projections were 0.0778 Å and
0.6731 Å respectively; original and projected states are retained. These are not
changes to observed source heavy atoms or additions to generated ownership.

This establishes a working private coordinate-elimination route despite the
direct OpenMM massless-constraint limitation. It does not establish a general
H-safe solution, superiority over the old optimizer, a converged global minimum,
successful recovery of 1NON, coordinate accuracy or MD readiness. In particular,
the budget/optimizer/constraint model differs from the historical calibration;
do not attribute the changed H outcomes to a single cause from this smoke.
The historical **18/100 stays unchanged**. The new 1KTZ success is a separate
positive implementation control from the same retained proposal, not a new
independent reconstruction. OXT remains separately reported and lower priority.

Evidence:
- [Frozen reduced-coordinate control plan](../../tmp/anchor-constrained-control-2026-11-17/attempt01/plan.json).
- [Worker results and full per-phase numerical outcomes](../../tmp/anchor-constrained-control-2026-11-17/attempt01/summary.json).
- [Independent sparse-projection/custody receipt](../../tmp/anchor-constrained-control-2026-11-17/attempt01/independent-readback/receipt.json).

Additional objective/cache/nonfinite-state and independent sparse-auditor tests
pass. No extra iteration extension, restart, source repair or new proposal was
launched after these outcomes. Any 1NON continuation, window-release comparison
or broader qualification needs a separately declared scope; its failure remains
in the ledger rather than silently extending this control's budget.

Failure evidence:
- [All-100 force/contact scan summary](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/summary.json).
- [Exact per-cell forces, attractive pairs and hidden contacts](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/endpoint-scan.json).
- [2BE1 full phase history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/2BE1-w1-e3-k10000-history.json).
- [1NON fixed-anchor history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/1NON-fixed-history.json).
- [6Z38 terminal/stereochemistry history](../../tmp/anchor-failure-analysis-2026-11-17/attempt01/6Z38-fixed-history.json).
- [Fixed-contact diagnostic mechanism and preparation provenance](../../tmp/anchor-failure-analysis-2026-11-17/source-contact-provenance.json).
- [Deposited mmCIF contact proof](../../tmp/anchor-failure-analysis-2026-11-17/deposited-contact-proof.json).
- [Historical protein-only proposal-template proof](../../tmp/anchor-failure-analysis-2026-11-17/generator-context-provenance.json).
- [Exact saved-state THR A1 Cβ sign-change proof](../../tmp/anchor-failure-analysis-2026-11-17/thr-beta-proof.json).

Retained local evidence (inside the repository working directory; transfer these
artifacts explicitly for another host):

- [Finite input/cell inventory](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/inventory.json).
- [Guarded declaration and exact implementation/input hashes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/plan.json).
- [Worker outcome summary](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/preflight-summary.json).
- [Independent readback and OpenMM XML resource hashes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt02/independent-readback/receipt.json).
- [Explicit-terminal preflight outcomes](../../tmp/anchor-calibration-2026-11-14/preflight-attempt04/preflight-summary.json).
- [Four-case independent particle/base-force-field readback](../../tmp/anchor-calibration-2026-11-14/preflight-attempt04/independent-readback-v2/receipt.json).
- [Frozen 100-cell execution manifest](../../tmp/anchor-calibration-2026-11-14/grid-attempt01/plan.json).
- [Final no-minimization adapter verification](../../tmp/anchor-calibration-2026-11-14/adapter-final-verification/receipt.json).
- [Successful complete fast-suite report](../../tmp/anchor-calibration-2026-11-14/fast-tests-path-fixed.xml).
- [Original complete grid ledger (scope caveat above)](../../tmp/anchor-calibration-2026-11-14/grid-attempt01/summary.json).
- [Corrected independent 100-cell readback](../../tmp/anchor-calibration-2026-11-14/grid-readback-attempt01/receipt.json).
- [Rejected-state native fixation and D-Cα confirmation](../../tmp/anchor-calibration-2026-11-14/rejected-state-readback-attempt01/receipt.json).
- [Paired trade-off map, 100 CSV rows](../../tmp/anchor-calibration-2026-11-14/tradeoffs-attempt01/tradeoffs.csv).
- [Bound map summary and limitations](../../tmp/anchor-calibration-2026-11-14/tradeoffs-attempt01/summary.json).

## 4. Continuation contract

1. Inspect repository changes and retained evidence. Preserve uncommitted work,
   original source bytes, frozen requests and rejected attempts. Start from the
   retained list, not a new collection campaign or a compulsory review backlog.
2. For each next task, identify the experimental accession, exact known construct
   sequence, missing residues, complete relevant context and computationally
   movable atoms. Resolve task-specific ambiguity without rewriting history.
3. Bind any new execution to an explicit input/backend/seed/budget declaration
   under the user's actual execution authorization. Do not infer blanket approval
   from an old autonomous-run note or from this documentation request.
4. Preserve native/raw outputs before source-frame materialization. Independently
   audit actual observation fixation before restoring source records; copying
   observations back does not prove that a backend kept them fixed.
5. Change only the declared missing-coordinate atoms. Any permitted observed-anchor
   motion requires a separate explicit scope and movement audit. Keep source links,
   waters and other relevant context; do not silently crop a frozen request.
6. Apply consistent full diagnostics to the original source and completed result:
   inventory, junctions, stereochemistry including sidechain centres, global
   L-chirality, declared links, Ramachandran/chi12 and steric findings. Separate
   unchanged source findings, newly introduced defects and later preparation needs.
   A narrow-validator pass is not a clean full-diagnostic result.
7. Retain native/raw, materialized and refined stages and all failures. Use
   `DatRecord` for exact completed-atom ownership. Do not tune historical thresholds
   to accept a result or infer accuracy from a geometry pass.

Continue on Mac; transfer to Linux/CUDA is not the current objective. Preserve
declared host/resource limits, one-ML-process concurrency, native-MPS labeling and
no silent CPU fallback. Use `--no-solvent` for development `minimize`/`zbs` runs.
Do not start new downloads, remote sequence uploads, resource-limit calibration,
checkpoint changes or broad job matrices under this handoff alone.

Public CLI/API/Homology/GUI capability changes, production-default changes and
broader scientific qualification remain separate work. The original handoff
recorded 2,848 fast tests passed, three optional skips and 33 slow tests deselected;
retain that historical result separately from the current continuation verification
linked above.

## 5. Historical evidence and implementation references

- [Real-gap history and work packages](mac-reconstruction-real-gaps.md): retained
  inventories, experiments and older proposals; current directions above take
  precedence over contradictory historical next steps.
- [Original multi-backend portfolio design](mac-reconstruction-portfolios-handoff.md):
  comparison methodology, candidate lineage and earlier bounded studies.
- [Private RFD3 source-link custody](rfd3-source-link-custody.md): native graph
  transport, link provenance and derivative checks. A custody pass does not erase
  the original native graph-transport failure.
- [Dated execution log](../research/mac-real-gap-continuation-2026-11-12.md):
  chronology, negative outcomes and readback references.
- [Structure-preparation context](../agent/contexts/structure-preparation.md):
  code ownership, contracts and entry points.

Historical approval, refusal, launch and pending-decision statements describe
their own stages. Keep them readable; do not treat each as a new current command.
