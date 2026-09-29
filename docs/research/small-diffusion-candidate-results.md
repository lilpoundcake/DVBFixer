# Small Diffusion Candidate Results

- Status: active exploratory evaluation.
- Evaluation date: 2026-09-29.
- Scope: compact pretrained models for one-chain, two-anchor internal gaps.
- Platform order: Linux CPU/CUDA viability first; Apple Silicon portability only
  for candidates that pass the Linux constrained-gap pilot.

## Strict All-Atom Classification

For this study, `all-atom` means that the selected generative model directly
produces every canonical non-hydrogen atom of each generated residue. A backbone
model followed by ProteinMPNN, Rosetta, `cg2all`, ESMFold, or another side-chain
builder does not qualify. Hydrogens remain a downstream preparation concern.

### Retained Compact Candidates

| Model | Approximate size | Native output | Current disposition |
|---|---:|---|---|
| Protpardelle-1c `cc89_epoch415` | 25.17M | Canonical atom37 heavy atoms with supplied sequence | Passed the first prescribed-sequence insertion pilot |
| Protpardelle-1c `cc91_epoch383` | 22-25M plus MiniMPNN | Canonical atom37 heavy atoms with sequence codesign | Retain only if the target sequence can be fixed |
| Protpardelle-1c `cc91_tip_epoch480` | 22-25M plus MiniMPNN | Canonical atom37 heavy atoms with side-chain-tip conditioning | Specialized reserve |
| Protpardelle-1c `cc94_epoch3100` | 22-25M plus MiniMPNN | Canonical atom37 heavy atoms, including multichain training | Multichain reserve |
| Original Protpardelle all-atom | approximately 22M | Canonical atom37 heavy atoms | Legacy same-length inpainting reference |

The Protpardelle all-atom representation excludes hydrogens. `cc91`, `cc91_tip`,
and `cc94` require MiniMPNN during sequence/structure codesign, and released
sampling configurations normally use full ProteinMPNN for the final sequence.
That dependency is scientifically unacceptable for DVBFixer unless the known
FASTA/SEQRES identities can be imposed instead of redesigned.

### Excluded Backbone Or Hybrid Models

| Model | Approximate size | Reason for exclusion |
|---|---:|---|
| Protpardelle-1c `cc58`, `bb81`, `bbmd`, `cc78`, `cc83`, `cc95` | 22-25M | N/CA/C/O backbone only |
| FoldingDiff | 14.5M | Backbone internal coordinates only |
| FrameDiff | 17.4M | Residue frames/backbone only |
| FrameDiPT | 17.4M | Backbone inpainting; requires external `cg2all` completion |
| Genie | 4.1M | C-alpha only |
| Genie 2 | 15-16M | C-alpha only |
| Genie 3 v1 | 23.27M | Generated scaffold residues are C-alpha tokens; atom14 is used for conditioned motifs or a separate side-chain-packing task |
| FoldFlow | 17.4M | Residue frames/backbone only |
| RFdiffusion | approximately 60M | Backbone only and above the compact limit |
| Chroma | unverified | Backbone diffusion plus a separate sequence/side-chain network; restricted weights |

Protenix Tiny and Mini are native all-heavy-atom predictors but are not compact:
the official counts are 109.50M and 134.06M parameters respectively. Their
released compact checkpoints also lack documented exact-coordinate conditioning,
so they remain quality references rather than candidates in this track.

## Protpardelle-1c

### Artifact Record

- Source: `https://github.com/ProteinDesignLab/protpardelle-1c.git`.
- Source revision: `ee378400f25b801fa481028000f9060183d7fb4c`.
- Package version: `1.3.2`.
- Code license: MIT.
- Model record: Zenodo `10.5281/zenodo.16817230`.
- Weight license: CC-BY-4.0 in the Zenodo record.
- Published archive bytes: `4,294,232,741`.
- Published archive MD5: `59fc043a9d981d6c74baa2d051ea164a`, verified.
- Evaluated checkpoint: `cc58_epoch416.pth`.
- Checkpoint SHA-256:
  `e92922d29336ca40f93685d547880fd76b0043f49415783de900a66ba258ed37`.
- Config SHA-256:
  `e7e7757e7aa1e7f1b4eb2a512919d3715db0a98fb87fc0fc73198d3757111f81`.
- Evaluated all-atom checkpoint: `cc89_epoch415.pth`.
- All-atom checkpoint SHA-256:
  `dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`.
- All-atom config SHA-256:
  `e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d`.
- Exact training membership is not published. The released config names an
  augmented AI-CATH/Ingraham CATH dataset, so leakage-sensitive selection remains
  blocked until membership or a conservative temporal/cluster exclusion is available.

### Runtime Record

- Host: Linux x86_64, NVIDIA A100-SXM4-40GB.
- Environment: Python 3.12.14, PyTorch 2.6.0+cu124, CUDA available.
- Installed isolated environment: 6.2 GB. This includes unnecessary upstream
  Jupyter, Transformers, W&B, and CUDA development/runtime dependencies and is not
  a minimized deployment measurement.
- Instantiated `cc58` parameter count: `25,169,439`, inside the declared 15-50M
  `small` class but larger than the paper abstract's approximate 22M description.
- CPU one-step forward, batch 1, length 32: 0.324 seconds; finite output.
- A100 one-step forward, batch 1, length 128: 0.209 seconds; approximately 155 MiB
  peak allocated model/forward memory; finite output.
- A100 pretrained unconditional sample, length 32, 500 steps: 10.13 seconds.
- Instantiated `cc89` parameter count: `25,169,439`.
- `cc89` CPU one-step forward, batch 1, length 32: 0.140 seconds; finite
  atom37 output.
- `cc89` A100 one-step forward, batch 1, length 128: 0.213 seconds;
  approximately 157 MiB peak allocated model/forward memory; finite atom37
  output with shape `(1, 128, 37, 3)`.
- No upstream source patch was required for these Linux smokes.

### Backbone Three-Residue Withheld Pilot

The pilot withheld residues 51-53 from chain B of the released 128-residue
`7eow_B_atom.pdb` example. Residues 1-50 and 54-128 were supplied through native
side-chain motif conditioning. One sample used seed `20260929`, 500 steps,
`step_scale=1.2`, and `s_churn=200`.

- Sampling time: 10.58 seconds.
- Output residues: 128.
- Conditioned backbone fit after rigid alignment: 0.028 A RMSD over 500 N/CA/C/O
  atoms; maximum displacement 0.103 A.
- Withheld three-residue backbone RMSD to hidden coordinates: 0.166 A.
- Gap-only splice junction C-N distances against original fixed atoms: 1.307 A
  and 1.281 A, both inside the 1.20-1.45 A gate.
- Generated residues are GLY placeholders because this checkpoint predicts
  backbone coordinates only and MPNN was deliberately disabled.

This is a positive runtime/backbone-conditioning signal, not an acceptance result
for the strict all-atom study. Native
conditioning does not preserve fixed coordinates exactly. A DVBFixer adapter must
rigidly map the proposal back to the source frame, take only generated gap
coordinates, retain every original fixed atom unchanged, restore authoritative
FASTA/SEQRES residue identities, and complete canonical side chains before normal
geometry, chirality, provenance, and atomic-publication validation.

An eight-candidate batch with the same seed and settings completed in 11.83
seconds. After rigid mapping and gap-only splicing, withheld-backbone RMSD ranged
from 0.134 to 0.339 A with median 0.228 A. Seven of eight candidates passed both
the 1.20-1.45 A junction C-N gate; the failing candidate had one 1.471 A
junction. This supports generate-many/validate/rank behavior and confirms that a
single native sample cannot be published without independent validation.

### All-Atom Prescribed-Sequence Insertion Pilot

The `cc89_epoch415` pilot used the same 128-residue `7eow_B_atom.pdb`
structure, seed `20260929`, and 500 diffusion steps. Residues 1-50 and 54-128
were supplied as the two fixed motifs. Residues 51-53 were inserted through
the low-level `Protpardelle.sample` API with the authoritative `ILE-SER-ARG`
sequence supplied through `gt_aatype`; MiniMPNN and full ProteinMPNN were both
disabled.

- Output residues: 128 with the prescribed sequence unchanged.
- Canonical heavy-atom completeness: 972/972 over the structure and 25/25 in
  the three generated residues; no hydrogens were emitted.
- Fixed-motif proposal fit after rigid alignment: 0.070 A RMSD, with 0.293 A
  maximum displacement. Production publication must still splice only the gap
  so source fixed atoms remain exact.
- Generated gap all-heavy-atom RMSD to hidden coordinates: 0.664 A.
- Gap-only splice junction C-N distances: 1.383 A and 1.342 A, both inside the
  1.20-1.45 A gate.
- `dvbfixer diagnose` found no chirality failure. Its only completeness finding
  was the expected absent terminal `OXT` from the upstream writer.

This passes the first native all-heavy-atom, prescribed-sequence insertion
gate. It does not establish general short/long-gap quality, leakage safety, or
production readiness. The upstream CLI does not expose this exact combination
of insertion motifs and a fixed sequence, so an adapter would need the
low-level API and must preserve atomic publication semantics.

### Frozen V5 Six-Case Pilot

The next pilot reused six requests from the completed v5 Protenix/Boltz/MODELLER
confirmatory materialization. The cases were frozen before Protpardelle outcomes
were inspected: the first three selected five-residue and first three selected
ten-residue cases in `screening_index` order within the released 512-residue
limit. This is adapter evidence, not an independent benchmark.

- Source aggregate SHA-256:
  `55e4cbc6c16ee6bf34a14aabf61ec71afcddfb436f15edb5195c3d8b3f717ad7`.
- Model/checkpoint/config/source pins are unchanged from the artifact record above.
- Runtime: Python 3.12.14, PyTorch 2.6.0+cu124, NVIDIA A100-SXM4-40GB.
- Adapter validation required MDAnalysis 2.10.0 to be added to the isolated
  Protpardelle environment; this did not change the DVBFixer core environment.
- Sampling: the request seed, 500 steps, `step_scale=1.2`, `s_churn=200`, supplied
  target sequence, and no MiniMPNN or ProteinMPNN.
- Final publication candidate: rigid synchronization followed by exact source
  fixed-atom reinjection; no generated-only boundary refinement in this ablation.
- Same-seed repeat on `36hb`: byte-identical candidate PDB with SHA-256
  `4ce30b48acf160ed9874561f193b64e6ff740f9af3a8f33c6b4f7827c98e5718`.

| Case | Gap | Protpardelle pass | Backbone RMSD | All-heavy RMSD | Protenix pass/RMSD | Boltz pass/RMSD | MODELLER pass/RMSD |
|---|---:|---:|---:|---:|---:|---:|---:|
| `36hb` | 5 | no: overlap | 0.244 A | 1.680 A | yes / 0.237 A | yes / 0.316 A | yes / 4.986 A |
| `9gtp` | 5 | no: overlap | 0.505 A | 1.121 A | yes / 0.491 A | yes / 0.562 A | no / 8.913 A |
| `9gae` | 5 | yes | 0.496 A | 1.101 A | yes / 0.204 A | yes / 0.268 A | yes / 3.133 A |
| `9ina` | 10 | no: overlap | 0.556 A | 1.279 A | yes / 0.338 A | yes / 0.352 A | yes / 9.945 A |
| `9eho` | 10 | no: overlap | 1.633 A | 2.046 A | yes / 0.949 A | yes / 1.533 A | yes / 4.502 A |
| `9gbg` | 10 | yes | 1.520 A | 1.890 A | yes / 1.426 A | yes / 0.957 A | yes / 7.204 A |

Protpardelle passed the complete scientific geometry contract in `2/6` cases:
`1/3` short and `1/3` long. All six preserved exact fixed coordinates, supplied
residue identities, and complete requested canonical generated heavy atoms.
The four failures each had only `severe-steric-overlap`; no case failed junction
connectivity, internal backbone continuity, bond geometry, peptide planarity,
Ramachandran, side-chain chi1/chi2, or C-alpha chirality. Median gap-backbone
RMSD was 0.531 A, median gap all-heavy RMSD was 1.480 A, median complete adapter
wall time was 11.05 seconds, and reported peak allocated VRAM was 302,253,056
bytes. The six candidate SHA-256 values, in table order, are:

- `4ce30b48acf160ed9874561f193b64e6ff740f9af3a8f33c6b4f7827c98e5718`
- `87d3906f6c30228f953d72abea5cde5734d419a8c961b28fa7a551953ce8377d`
- `a2ceddea80658c6c819a5c19ca9f0d929346fb3d40b8eb9d475e0859296929c6`
- `51ae3f22dd1b33937a4e2325085b6443dd3705118b7310fc18aab341db7776f8`
- `50b944e3846e95b3876485d35bdb5096ab6ed0f63d850fba4f6465b362e96cb8`
- `e0ed2f632971343db8da332ddcef0387a7675ed79d85a8e2189fa1fd8d90171e`

The pilot exposed one adapter defect before scientific aggregation: `9gtp`
contains a fixed terminal `OXT`, which is outside Protpardelle's atom37 state.
The adapter now aligns on represented fixed atoms and restores every source-only
fixed atom exactly, while still requiring the model to supply every generated
atom. A dependency-light regression test covers this behavior. The failed
pre-fix run is excluded because it never produced a candidate; the successful
rerun used the unchanged seed and sampling settings.

The predefined generated-only boundary-refinement follow-up reused these exact
six raw candidate digests without resampling. OpenMM `CPU` used eight threads,
eight deterministic input-seed restarts, and boundary-refinement revision
`dvbfixer-openmm-boundary-refinement-v4`. The CPU platform is not claimed to be
byte-deterministic.

| Case | Refined pass | Refined backbone RMSD | Refined all-heavy RMSD | Refinement time | Refined SHA-256 |
|---|---:|---:|---:|---:|---|
| `36hb` | yes | 0.274 A | 1.554 A | 17.09 s | `1fe6bb6a403df450665dd8d07679b6745ea9f4f21e2d28aa7fe44ddfa45b80da` |
| `9gtp` | yes | 0.560 A | 1.237 A | 32.13 s | `1948a998a232aeed1f7798022f1717380f1a74517c03d421eede53f6acba7823` |
| `9gae` | yes | 0.254 A | 1.045 A | 13.73 s | `ae772f547ba2cca660ee889a3b675d41297fce7e1b1c98e4506d4ba2770d4fd2` |
| `9ina` | yes | 0.401 A | 1.118 A | 73.20 s | `4069871c598bf853e132840cee4b431d03ec73d7e6d747b0636ed21a56841f51` |
| `9eho` | yes | 1.565 A | 1.933 A | 15.84 s | `799b07ae09285c88b055bb0c309d7e102f4ff9b56c9e62c68bf3405f774363a7` |
| `9gbg` | yes | 1.514 A | 1.892 A | 91.92 s | `a68c50cbf8e851d74fd0164f4a2520a1f08d48c05f1db641432fd65cf7b38804` |

Refinement raised the scientific hard-gate pass count from `2/6` to `6/6` while
keeping fixed-heavy RMSD exactly zero in every case. Median gap-backbone RMSD
improved from 0.531 A to 0.480 A and median all-heavy RMSD from 1.480 A to
1.396 A. Median added refinement time was 24.61 seconds. All four raw overlap
failures were removed without introducing a different hard-gate failure.

Two platform failures remain operational evidence. OpenMM `Reference` refined
all three five-residue cases to passing candidates in 47-126 seconds, but the
440-residue `9ina` run exceeded a 900-second timeout. OpenMM `CUDA` originally
failed before coordinate refinement with `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`.
The mismatch is now diagnosed as CUDA/NVRTC 12.9 running against a 535-series
driver with CUDA 12.2 compatibility. An isolated OpenMM 8.2/CUDA 12.2 Context
smoke passes; real-candidate tests had not yet run at this documentation
checkpoint. See [`openmm-cuda-refinement.md`](openmm-cuda-refinement.md). CPU
therefore remains the only completed same-platform six-case refinement result;
this does not establish Apple runtime feasibility.

### Expanded V5 Linux Pilot

The preregistered 24-case extension is recorded in
[`small-diffusion-v5-expanded-pilot.json`](small-diffusion-v5-expanded-pilot.json),
SHA-256 `fbeb95ea342ef706d123e3d22ef61a0e836ea65a50ed67af3355c21d6b79b3d3`.
It selected the next 12 five-residue and 12 ten-residue v5 cases by
`screening_index` after excluding the original six, without outcome-based
replacement. All 24 raw and CPU-refinement stages completed through digest-based
resumable orchestration.

| Stage | All | Gap 5 | Gap 10 | Median backbone RMSD | Median all-heavy RMSD |
|---|---:|---:|---:|---:|---:|
| Raw gap-only reinjection | 13/24 | 7/12 | 6/12 | 0.982 A | 1.868 A |
| Refined candidate alone | 21/24 | 12/12 | 9/12 | 0.796 A | 1.742 A |
| Validation-first selection | 22/24 | 12/12 | 10/12 | 0.684 A | 1.627 A |
| Protenix v1 comparator | 23/24 | 12/12 | 11/12 | 0.718 A | not reported |
| Boltz-2 proxy | 23/24 | 12/12 | 11/12 | 0.549 A | not reported |
| MODELLER 10.8 comparator | 21/24 | 12/12 | 9/12 | 5.637 A | not reported |

Validation-first selection means: publish the refined candidate when it passes;
otherwise retain the raw candidate only when the raw candidate passes; otherwise
count the case as failed. This selected 21 refined candidates and one raw
candidate. `9hab` is the important regression case: raw passed, refinement added
a severe overlap, and post-refinement validation correctly retained raw instead.
The two remaining failed cases are `9g9t` and `9e51`. Both retain peptide-
connectivity and local-geometry failures after refinement; they were not replaced.

Against Protenix, 22 cases passed both, Protenix alone passed one, and neither
passed one. Median selected-Protpardelle minus Protenix backbone RMSD was
`+0.080 A` among common passing cases. Against Boltz, 21 passed both,
Protpardelle alone passed one, Boltz alone passed two, and neither count was zero;
the median paired difference was `+0.185 A`. Against MODELLER, 19 passed both,
Protpardelle alone passed three, MODELLER alone passed two, and neither count was
zero; the median paired difference was `-4.173 A` for Protpardelle.

Native conditioning was not exact and therefore cannot replace reinjection.
After rigid alignment but before reinjection, median fixed-atom RMSD was
`0.549 A` and median maximum displacement was `4.654 A`. Final fixed-heavy RMSD
was zero in every raw and refined result. Median raw model time was 11.13 seconds,
median CPU refinement time was 55.05 seconds, median complete driver wall time
was 73.15 seconds, maximum complete wall time was 442.63 seconds, and peak model
VRAM was 302,253,056 bytes. The ignored aggregate is
`.artifacts/diffusion-confirmatory-cohort-500-v5/small-diffusion-v5-expanded-report.json`,
SHA-256 `fcde508ace08ecdc7576b1610e74b653772e0fdf28515a8d2c4916bca2c37aa0`.

The expanded pilot supports continued compact-model research but not backend
selection. The sample is small, training membership remains unresolved, Boltz is
proxy-only, no uncertainty interval is meaningful at this size, and the two long-
gap failures show that boundary refinement is not a universal repair.

### Per-Step Reinjection Ablation

The bounded same-seed `9dvi` experiment was frozen before callback inference in
[`small-diffusion-protpardelle-reinjection-ablation.json`](small-diffusion-protpardelle-reinjection-ablation.json),
SHA-256 `20f48f7959eb84ba6ac537157c850fea6e99ac5f279a10e0af0d129075847b71`.
It retained the existing request seed, 500 steps, `step_scale=1.2`, `s_churn=200`,
fixed sequence, and disabled MPNN paths.

The first preregistered implementation exposed a coordinate-frame defect rather
than a model result. It Kabsch-aligned every updated state to source PDB coordinates
while Protpardelle continued conditioning against its randomly rotated internal
motif. Although all 500 callbacks completed with zero projection error, the two
incompatible frames destroyed the generated region. That artifact is retained at
`protpardelle-reinjection-1` and is excluded from model comparison.

The corrected hook passes Protpardelle's internal `motif_all_atom` coordinates and
atom mask to the post-update callback. The callback overwrites fixed atoms in that
same frame before the next step; publication still performs the existing final
rigid synchronization and exact source-coordinate restoration. This is a disclosed
post-failure implementation correction, not a pristine preregistered result.

| Mode | Pass | Backbone RMSD | All-heavy RMSD | Fixed proposal fit | Candidate SHA-256 |
|---|---:|---:|---:|---:|---|
| Native conditioning plus final projection | no: junction, Ramachandran | 1.751 A | 2.250 A | 0.530 A RMSD, 4.807 A max | `b3af11aded576cdcc49fb9dff6a8e37979843286745076cff9a477f462902251` |
| Corrected per-step reinjection | no: junction, Ramachandran | 1.760 A | 2.271 A | 0.000008 A RMSD, 0.000014 A max | `cf9c2569fe763dc7bf19cbe3bda00b1d95a552718d6bb715c33d01bb40280cff` |
| Per-step reinjection plus CPU refinement, repeat 1 | yes | 1.397 A | 1.775 A | fixed-heavy RMSD 0.0 A | `0bf13ac5abafd173b3abb6004cd860122f953af2f5d5ff4de0fef1303a381a5c` |
| Per-step reinjection plus CPU refinement, repeat 2 | yes | 1.363 A | 1.676 A | fixed-heavy RMSD 0.0 A | `d627de4812f96241e7d9b96b8374e245a081816276865fbfb3496012eb6424c4` |

Two independent corrected sampler runs were byte-identical. Each invoked the
callback exactly 500 times, projected 569 represented fixed atoms with zero
post-projection error, emitted every requested generated heavy atom, and preserved
fixed-heavy RMSD at zero. Per-step reinjection therefore supplies exact trajectory
control but did not repair this case's raw geometry or materially improve withheld
RMSD.

Both generated-only CPU refinements passed all hard gates, but they were not
repeatable: generated-atom RMSD between repeats was `0.274 A` with `0.722 A`
maximum displacement; whole-candidate RMSD was `0.065 A`. A separate repeat of
the original native-conditioning candidate was also nondeterministic (`0.104 A`
whole-candidate RMSD). OpenMM CPU refinement must therefore be recorded as measured
nondeterministic on this host rather than described as deterministic from its input
restart seeds.

The corrected callback was then run as a separate follow-up over all 24 frozen
expanded-pilot cases. This reused the frozen case identities and comparator values,
but it is not an independent benchmark: the baseline outcomes and the one-case
callback ablation were already known. Outputs live under each workspace's
`small-diffusion-expanded-reinjection-v2` directory. The aggregate is
`.artifacts/diffusion-confirmatory-cohort-500-v5/small-diffusion-v5-expanded-reinjection-report.json`,
SHA-256 `2309ccce2d32aa67325661ace348dd2b57a63163b9e87d6240429651b4619ff0`.

| Stage | All | Gap 5 | Gap 10 | Median backbone RMSD | Median all-heavy RMSD |
|---|---:|---:|---:|---:|---:|
| Per-step raw | 13/24 | 6/12 | 7/12 | 0.868 A | 1.804 A |
| Per-step refined | 21/24 | 11/12 | 10/12 | 0.808 A | 1.733 A |
| Per-step validation-first selection | 21/24 | 11/12 | 10/12 | 0.718 A | 1.640 A |
| Earlier final-projection selection | 22/24 | 12/12 | 10/12 | 0.684 A | 1.627 A |
| Protenix v1 comparator | 23/24 | 12/12 | 11/12 | 0.718 A | not reported |
| Boltz-2 proxy | 23/24 | 12/12 | 11/12 | 0.549 A | not reported |
| MODELLER 10.8 comparator | 21/24 | 12/12 | 9/12 | 5.637 A | not reported |

Every raw run invoked exactly 500 callbacks, reported zero post-projection error,
preserved fixed-heavy RMSD at zero, and completed without an operational failure.
The corrected representable-atom mapping also passed all three cases with a
source-only terminal `OXT`; those atoms stayed on the final restoration path rather
than being misclassified as model-generated coordinates.

Against the earlier final-projection pipeline, 21 cases passed both selections.
The median paired backbone-RMSD difference was `-0.0006 A` for per-step
reinjection, with 11 lower and 10 higher RMSDs. `9ges` regressed from a passing
refined baseline to a generated-bond-length failure, reducing selected coverage by
one case; `9g9t` and `9e51` remained failed. Against Protenix, per-step selection
had 21 both-pass, zero Protpardelle-only, two Protenix-only, and one neither-pass
cases. The corresponding counts were 20/1/3/0 against Boltz and 18/3/3/0 against
MODELLER.

Per-step reinjection is therefore retained as a verified control mechanism, not
selected as the default compact-model sampling mode. It makes the trajectory's
fixed coordinates exact but provides no paired RMSD benefit and slightly worsens
validation-first pass coverage on this follow-up. The simpler native-conditioning
plus final exact projection path remains the compact exploratory baseline.

## Genie 3

### Artifact And Architecture Record

- Source: `https://github.com/aqlaboratory/genie3.git` at revision
  `d77ae5ac04212ff1e8b29b585859a3244c614804`.
- Model repository: `https://huggingface.co/yeqinglin/genie3` at revision
  `9ae31ebb8c56eebdc05ab282a8fd3f6a6d2a03a2`.
- Code and model license: Apache-2.0.
- Evaluated v1 checkpoint: `step=600000.ckpt`, 346,849,270 bytes.
- Checkpoint SHA-256:
  `57039afcd4e5a66e58c52ba1ece256bc311a4294a31c37fca7301dd3a4c4479e`,
  matching the Hugging Face LFS digest.
- The model state has 504 tensors and instantiates as exactly `23,269,484`
  trainable parameters. All checkpoint keys load strictly into the released v1
  architecture under PyTorch 2.6.0+cu124; model construction and CPU checkpoint
  loading took 8.97 seconds.
- The model card reports training on AFDB representatives and Pinder.

Genie 3 is excluded before a sampling benchmark. In the released feature
pipeline, every unconditioned scaffold residue has `residue_cond_group == 0`
and becomes one C-alpha token regardless of `prot_rep_mode=atom14`. Atomized
tokens represent conditioned motif atoms. The separate sidechain task expands
residues to atom14 tokens only while conditioning on an already supplied full
backbone and sequence. Consequently, an internal-gap workflow would require a
backbone generation/reconstruction stage followed by side-chain packing; the
gap's complete canonical heavy-atom coordinates are not one native generative
state. That fails this study's strict direct all-heavy-atom gate despite the
broader all-atom-equivariance description of the architecture.

## Current Decision

The tested `cc58` checkpoint and Genie 3 are excluded because their generated
gap scaffold is backbone/C-alpha-only. The all-atom `cc89` checkpoint is the sole
primary compact candidate. It has now demonstrated repeatable prescribed-sequence
sampling, complete generated heavy atoms, and exact fixed-coordinate publication
on frozen short/long pilots. The expanded 24-case validation-first result is
`22/24`, close to Protenix's `23/24` and above MODELLER's `21/24`, with much lower
RMSD than MODELLER but slightly higher paired RMSD than Protenix. This is a
positive Linux adaptation signal, not readiness for CLI or Apple Silicon work:
exact training membership remains unavailable and no independent confirmatory
study has been run. The next decision is whether this leakage-exploratory evidence
is sufficient to justify the separately planned Apple operator audit; it cannot
support backend selection. A future MPS audit must specifically inspect CUDA
autocast decorators, CUDA cache calls, unsupported operators, and silent CPU
fallback.
