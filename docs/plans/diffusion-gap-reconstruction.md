# Diffusion Gap-Reconstruction Implementation Plan

- Status: proposed.
- Branch: `diffusion`.
- Scope owner: Structure Preparation bounded context.
- Change group: `missing-atom-rebuild-chirality`.
- Coordination rule: assign one exclusive coordinating owner before implementation.
- Related research:
  - [`../research/reconstruction-and-modeling-backends.md`](../research/reconstruction-and-modeling-backends.md)
- Production baselines retained:
  - Salilab MODELLER for `model` and `homology`.
  - PDBFixer/OpenMM Modeller in the documented legacy preparation path.
  - `tleap-reduce` for its documented pure-protein preparation scope.

## Goals

- [ ] Add an experimental, backend-neutral diffusion path for sequence-guided protein gap reconstruction.
- [ ] Preserve deposited coordinates outside explicitly generated or movable regions.
- [ ] Implement the PATCHR-like method independently rather than importing, wrapping, vendoring, or copying PATCHR.
- [ ] Keep ML frameworks, CUDA libraries, and model-specific packages out of the core DVBFixer environment.
- [ ] Record enough provenance to reproduce or audit every generated candidate.
- [ ] Validate identity, connectivity, geometry, chirality, and clashes before publishing a candidate.
- [ ] Evaluate mosaic-first diffusion for comparative/homology modeling only after local gap reconstruction passes its acceptance gates.

## Non-goals

- [ ] Do not replace or deprecate MODELLER in this work.
- [ ] Do not change the default `model` or `homology` backend during the research phases.
- [ ] Do not add automatic fallback from diffusion to MODELLER.
- [ ] Do not use full-chain prediction as a substitute for fixed-coordinate gap repair.
- [ ] Do not treat missing heavy atoms inside an otherwise present residue as the same problem as missing whole residues.
- [ ] Do not move CIF parsing into the scientific modeling stage.
- [ ] Do not claim support for ligands, glycans, PTMs, cofactors, metals, or covalent non-peptide links until each class has separate evidence.
- [ ] Do not store model weights in Git.

## Initial Supported Slice

- [ ] Accept canonical L-protein inputs only.
- [ ] Require a complete, known target sequence from FASTA or SEQRES.
- [ ] Require unambiguous shared sequence placement.
- [ ] Support internal gaps with two observed anchors.
- [ ] Start with gap lengths of 3-12 residues.
- [ ] Reject terminal one-anchor gaps in the first slice.
- [ ] Reject alternate-location ambiguity in or adjacent to the generated region.
- [ ] Reject multiple `MODEL` blocks in the first slice.
- [ ] Reject retained unsupported heterogens or external covalent links near the generated region.
- [ ] Return an explicit `unsupported` result without creating a public PDB for every out-of-scope case.

## Mandatory Identity And Publication Rules

- [ ] Preserve case-sensitive chain IDs.
- [ ] Preserve residue insertion codes.
- [ ] Use `(chain, resid, icode, atom)` as atom identity.
- [ ] Use `(chain, resid, icode)` as residue identity.
- [ ] Preserve all atoms and explicit links outside the generated mask.
- [ ] Preserve `keepIds=True` on every OpenMM `PDBFile.writeFile` call.
- [ ] Normalize CIF input only at the CLI boundary.
- [ ] Use `dvbfixer.ffutils.dat.DatRecord` for the downstream `.dat` sidecar.
- [ ] Publish a candidate only after every hard validation gate passes.
- [ ] Write candidate PDB, `.dat`, and diffusion provenance atomically as one logical result.
- [ ] Leave no partial public output after runner, parsing, validation, or publication failure.

## Engine And License Decisions

### RFdiffusion v1 Baseline

- [ ] Use RFdiffusion v1 as the first external backbone-inpainting benchmark adapter.
- [ ] Pin the exact source revision.
- [ ] Pin and hash the exact checkpoint.
- [ ] Archive the applicable BSD-3-Clause license evidence for code and README-linked weights.
- [ ] Treat RFdiffusion output as a backbone-level benchmark until side-chain materialization and all-heavy-atom validation are complete.
- [ ] Do not expose raw RFdiffusion output as a production DVBFixer repair result.

### All-Atom Feasibility Spike

- [ ] Evaluate Protenix v1 first.
  - [ ] Pin a revision whose code and model-parameter terms are Apache-2.0.
  - [ ] Do not substitute Protenix v2 without a new license review.
  - [ ] Verify whether the denoising state can be intercepted before every next step.
- [ ] Evaluate Boltz-2 if Protenix v1 does not provide a maintainable hook.
  - [ ] Pin MIT code and weight revisions.
  - [ ] Treat ordinary template conditioning as insufficient for exact fixed-coordinate preservation.
  - [ ] Verify whether fixed-coordinate reinjection can be implemented without an unmaintainable upstream fork.
- [ ] Stop before user-facing integration if neither pinned engine provides stable denoising hooks.
- [ ] Record whether a maintained fork, upstream API change, or different sampler would be required.

### PATCHR Comparator

- [ ] Use PATCHR only as an external scientific comparator.
- [ ] Do not import PATCHR into DVBFixer.
- [ ] Do not vendor PATCHR source.
- [ ] Do not wrap the PATCHR CLI as the DVBFixer implementation.
- [ ] Require a separate ADR before changing this repository policy.

### Excluded Distribution Candidates

- [ ] Exclude FrameDiPT from the distributable backend because of CC BY-NC-SA terms.
- [ ] Exclude Chroma parameters from the distributable backend because of non-commercial parameter terms and API-key distribution.
- [ ] Exclude RFdiffusion2 until official checkpoint rights and optional Chai dependencies are unambiguous.
- [ ] Exclude engines with unclear checkpoint redistribution rights from the shipped adapter list.
- [ ] Keep unconditional-generation or motif-scaffolding engines as research comparators unless they demonstrate the required two-anchor repair contract.

### Per-Artifact License Inventory

- [x] Record code repository and revision.
- [x] Record source-code license and exact license-text revision.
- [ ] Record checkpoint URL and SHA-256. URLs are recorded where published; required hashes remain explicit blockers until artifact acquisition.
- [ ] Record checkpoint license and redistribution decision. Upstream claims and unresolved review status are recorded; no engine is approved for distribution.
- [ ] Record auxiliary model and cache URLs, hashes, and licenses. Required auxiliary artifacts remain unresolved.
- [ ] Record container base-image and package-lock provenance. Required digests and lock hashes remain unresolved.
- [x] Record relevant training-data cutoff claims for leakage analysis, including unresolved sequence-level membership.

## Backend-Neutral Contract

- [x] Add `src/dvbfixer/model/diffusion/contract.py`.
- [x] Define typed `AtomIdentity` and `ResidueIdentity` values.
- [x] Define a `GapRegion` with:
  - [x] target-chain identity;
  - [x] target-sequence interval;
  - [x] left and right anchors;
  - [x] generated residues;
  - [x] movable junction window.
- [x] Define a versioned `DiffusionRequest` with:
  - [x] normalized PDB input;
  - [x] target sequences;
  - [x] shared sequence placement;
  - [x] explicit gap ranges;
  - [x] fixed and generated atom masks;
  - [x] retained explicit links;
  - [x] candidate count;
  - [x] seed list;
  - [x] backend options.
- [x] Define a `DiffusionCandidate` with:
  - [x] candidate coordinate artifact;
  - [x] generated atom and residue identities;
  - [x] raw backend score;
  - [x] score provenance;
  - [x] warnings.
- [x] Define a `DiffusionResult` with:
  - [x] `success`, `unsupported`, or `failed` status;
  - [x] candidates;
  - [x] validation summaries;
  - [x] runner diagnostics;
  - [x] backend provenance.
- [x] Add strict JSON serialization and schema-version validation.
- [x] Reject unknown or incompatible schema versions explicitly.

## Mask And Sequence Placement Layer

- [x] Add `src/dvbfixer/model/diffusion/masks.py`.
- [x] Reuse the existing shared sequence-placement behavior instead of implementing a second alignment policy.
- [x] Construct fixed/generated masks from the authoritative alignment.
- [x] Preserve insertion-code-aware residue identity.
- [x] Preserve case-distinct chains such as `D` and `d`.
- [x] Define deterministic left and right anchor windows.
- [x] Define a separately bounded movable junction window.
- [x] Reject ambiguous placements instead of selecting a scientifically silent alternative.
- [x] Reject one-anchor and no-anchor regions in the initial slice.

## External Runner Boundary

- [x] Add `src/dvbfixer/model/diffusion/runner.py`.
- [x] Use a versioned subprocess protocol:
  - [x] DVBFixer writes `request.json` into an isolated workspace.
  - [x] The external runner writes a raw `RunnerResult` in `result.json` plus candidate artifacts; independent `DiffusionResult` validation summaries remain DVBFixer-owned.
  - [x] DVBFixer validates every returned path and digest.
- [x] Enforce workspace path containment.
- [x] Reject symlink and path-traversal escapes.
- [x] Bound runtime and captured output.
- [x] Capture exit code, stdout, and stderr without leaking credentials.
- [x] Treat malformed or missing result manifests as failures.
- [ ] Add preflight errors for missing runner, image, checkpoint, CUDA device, or incompatible protocol version. Missing executables and incompatible protocol versions fail now; image, checkpoint, and CUDA checks remain adapter-specific Phase 2 work.
- [x] Add a deterministic fake CPU runner for unit and integration tests.
- [x] Keep runner-specific Python packages outside the core DVBFixer dependency set.

## PATCHR-Like Constrained Sampling

- [ ] Build the complete target topology before sampling.
- [ ] Separate fixed observed atoms from generated missing-region atoms.
- [ ] Initialize or sample only the generated region.
- [ ] After every denoising update:
  - [ ] compute weighted Kabsch synchronization from reliable fixed anchors;
  - [ ] synchronize the generated frame to the experimental frame;
  - [ ] overwrite every fixed coordinate with its exact experimental coordinate;
  - [ ] verify that atom identities still match the request.
- [ ] Run a localized second pass around both peptide junctions.
- [ ] Keep atoms outside the declared movable junction window fixed during the localized pass.
- [ ] Materialize all expected canonical heavy atoms in generated residues.
- [ ] Leave hydrogen placement and protonation to the downstream preparation stage.
- [ ] Implement three explicit ablation modes:
  - [ ] template conditioning only;
  - [ ] template conditioning plus per-step reinjection;
  - [ ] reinjection plus local boundary refinement.

## Independent Validation And Ranking

- [x] Add `src/dvbfixer/model/diffusion/validate.py`.
- [x] Compare atom and residue identities outside the generated mask exactly.
- [x] Compare explicit links outside and across the generated region.
- [x] Measure fixed-heavy-atom RMSD and maximum displacement.
- [x] Validate peptide connectivity at both junctions.
- [ ] Validate generated bond lengths, angles, planarity, and atom completeness. Canonical heavy-atom completeness, generated C-N break gates, severe local bond-length outliers, broad canonical backbone/peptide-angle gates, and peptide-amide planarity are implemented through shared diagnose policies; finer residue-specific side-chain geometry remains pending.
- [x] Validate severe intra-region and region-context clashes.
- [ ] Validate Ramachandran and rotamer quality.
- [ ] Run `fix_ca_chirality` only as an explicitly recorded repair step. Validation currently rejects D geometry and performs no silent repair.
- [x] Run `assert_all_l` after the final heavy-atom coordinate change.
- [x] Reject any candidate that fails a hard gate.
- [x] Rank only passing candidates.
- [x] Keep backend confidence, closure metrics, geometry metrics, and physical energies as separately named score components. First-slice ranking keeps DVBFixer geometry metrics separate and does not use raw backend confidence as a shared-scale metric.
- [x] Do not compare backend confidence directly with MODELLER `molpdf` as though they share a scale.

## Provenance Manifest

- [ ] Write `<stem>.diffusion.json` beside each published diffusion candidate.
- [ ] Include a manifest schema version.
- [ ] Record DVBFixer version and commit.
- [ ] Record runner protocol version.
- [ ] Record engine repository and revision.
- [ ] Record container digest or environment-lock hash.
- [ ] Record checkpoint and auxiliary-artifact hashes and license identifiers.
- [ ] Record request hash, target sequence, gap masks, and fixed/generated identity sets.
- [ ] Record seed, device, precision, PyTorch, CUDA, driver, and deterministic flags.
- [ ] Record known nondeterministic kernels or operations.
- [ ] Record raw score and score provenance.
- [ ] Record every validation metric and pass/fail decision.
- [ ] Record final PDB and `.dat` SHA-256 values.
- [ ] Keep diffusion environment metadata out of `DatRecord` unless a later general sidecar decision explicitly changes that contract.

## Atomic Publication

- [ ] Add `src/dvbfixer/model/diffusion/pipeline.py`.
- [ ] Generate and validate candidates entirely inside a temporary workspace.
- [ ] Select the ranked passing candidates before touching public destinations.
- [ ] Build candidate-matched `.dat` records through `DatRecord`.
- [ ] Stage PDB, `.dat`, and provenance manifest in temporary destination files.
- [ ] Commit the complete candidate artifact set atomically.
- [ ] Remove staged artifacts after any failed commit.
- [ ] Preserve the source input on every failure path.

## Implementation Phases

### Phase 0: Policy, Corpus, And Thresholds

- [x] Update the research note with engine, licensing, and hardware findings.
- [x] Add or update the DDD implementation task without presenting the backend as shipped.
- [x] Accept an ADR covering:
  - [x] subprocess/container isolation;
  - [x] separate diffusion provenance;
  - [x] no automatic scientific fallback;
  - [x] MODELLER remaining the default;
  - [x] independent PATCHR-like implementation.
- [x] Select reviewed benchmark fixtures.
- [x] Record fixture provenance and checksums.
- [x] Record exact masks and expected supported/unsupported classification.
- [x] Record upstream training cutoffs and sequence-identity leakage metadata where available.
- [x] Freeze numerical thresholds before comparative inference runs.

### Phase 1: CPU-Testable Core

- [x] Implement the versioned request/result contract.
- [x] Implement deterministic mask construction.
- [x] Implement weighted Kabsch and coordinate-reinjection primitives.
- [x] Implement the isolated runner protocol.
- [x] Implement the deterministic fake runner.
- [ ] Implement independent validation and ranking. Identity, drift, completeness, connectivity, severe local bond-length, broad canonical backbone/peptide-angle, peptide-planarity, clash, chirality, and passing-candidate ranking gates are implemented; finer residue-specific geometry, Ramachandran, and rotamer validation remain pending.
- [ ] Implement atomic publication.
- [ ] Implement separate provenance manifests.
- [x] Keep public CLI behavior unchanged through the runner slice.
- [x] Keep MODELLER execution unchanged through the runner slice.

### Phase 2: RFdiffusion v1 GPU Baseline

- [ ] Build a pinned Linux/NVIDIA runner environment.
- [ ] Add checksum-verified checkpoint acquisition.
- [ ] Map `DiffusionRequest` to RFdiffusion contig/inpainting inputs.
- [ ] Map generated coordinates back to stable DVBFixer identities.
- [ ] Materialize canonical side chains through an explicitly recorded pure-protein route.
- [ ] Re-run all all-heavy-atom validation after side-chain materialization.
- [ ] Keep outputs inside benchmark workspaces until hard gates pass.
- [ ] Measure fixed-atom drift, closure, quality, runtime, RAM, VRAM, and seed variability.

### Phase 3: All-Atom Constrained Sampler

- [ ] Complete the Protenix v1 hook feasibility spike.
- [ ] Complete the Boltz-2 hook spike only if needed.
- [ ] Select a sampler only after per-step control is demonstrated.
- [ ] Implement weighted Kabsch synchronization.
- [ ] Implement exact per-step fixed-coordinate reinjection.
- [ ] Implement localized boundary refinement.
- [ ] Run the three-way ablation benchmark.
- [ ] Require complete canonical heavy atoms.
- [ ] Require final `assert_all_l` after every candidate's last heavy-coordinate change.

### Phase 4: Experimental `model` CLI

- [ ] Add `--backend {modeller,diffusion}` only after Phase 3 hard gates pass.
- [ ] Keep `modeller` as the default.
- [ ] Add diffusion runner/checkpoint/device options in a separate argparse group.
- [ ] Reject incompatible existing options before launching the runner.
- [ ] Preserve FASTA/SEQRES and content-selection semantics where the diffusion scope supports them.
- [ ] Do not fall back automatically to MODELLER.
- [ ] Add backend preflight to `doctor` without removing its stable report sections.
- [ ] Regenerate CLI reference and GUI command schema through their generators.
- [ ] Do not propagate diffusion options through `zbs` in this phase.

### Phase 5: Mosaic-First Homology

- [ ] Keep `selected_template_mosaic.pdb` as the authoritative coordinate frame.
- [ ] Export companion coverage metadata from `materialize_template_plan` instead of recomputing template ownership later.
- [ ] Preserve zero-based half-open template-plan masks.
- [ ] Mark covered template atoms as fixed.
- [ ] Generate only uncovered insertions, substitutions, and bounded junction windows.
- [ ] Start with one-chain internal insertions.
- [ ] Add multi-chain and multi-template cases only after one-chain acceptance.
- [ ] Preserve distinct antibody H/L chains and insertion codes.
- [ ] Compare against MODELLER using the same target and template plan.
- [ ] Use OpenFold or Boltz full-chain prediction only as an independent plausibility comparator.
- [ ] Do not expose a public homology diffusion backend until mosaic adherence and multi-chain gates pass.

### Phase 6: Production Decision

- [ ] Complete code, weight, auxiliary-artifact, and container license review.
- [ ] Complete representative GPU deployment and reliability evaluation.
- [ ] Compare against MODELLER across every claimed supported stratum.
- [ ] Document supported and unsupported inputs.
- [ ] Make production promotion a separate evidence-backed decision.
- [ ] Do not remove MODELLER or PDBFixer as a side effect of promotion.

## Benchmark Corpus

- [ ] Add withheld-coordinate internal gaps of 3-5 residues.
- [ ] Add withheld-coordinate internal gaps of 6-12 residues.
- [ ] Add a separate research stratum for gaps of 13-25 residues.
- [ ] Include regular loops.
- [ ] Include glycine-rich and proline-rich difficult loops.
- [ ] Include interface-adjacent gaps.
- [ ] Include antibody insertion codes.
- [ ] Include case-distinct chain IDs.
- [ ] Add terminal one-anchor gaps only as a separate later stratum.
- [ ] Add retained heterogen and covalent-link cases only as separate later strata.
- [ ] Keep reviewed structures under `tests/fixtures/`.
- [ ] Document fixture provenance in `tests/fixtures/README.md`.
- [ ] Regenerate `tests/fixtures/MANIFEST.sha256` after fixture changes.

## Initial Hard Gates

- [x] Preserve 100% of atom and residue identities outside the generated mask.
- [x] Preserve 100% of applicable explicit links.
- [x] Require fixed-heavy-atom RMSD no greater than `0.01 Å` after final serialization.
- [x] Require fixed-heavy-atom maximum displacement no greater than `0.03 Å`.
- [x] Require zero detectable D-C-alpha centers.
- [x] Require complete expected heavy atoms for every generated canonical residue.
- [x] Require both junction C-N distances to fall within `1.20-1.45 Å`.
- [x] Reject any backbone break greater than `1.8 Å` in the generated region.
- [x] Reject severe steric overlaps in the generated and junction neighborhoods.
- [ ] Require unsupported or ambiguous inputs to create no public PDB.
- [ ] Require same-seed repeatability within `0.01 Å` coordinate RMSD in the exact pinned environment, or explicitly classify and measure nondeterminism.
- [ ] Version threshold changes instead of adjusting them after viewing benchmark outcomes.

## Comparative Metrics

- [ ] Report gap backbone RMSD against withheld coordinates.
- [ ] Report gap all-heavy-atom RMSD.
- [ ] Report lDDT and GDT-HA or TM-score where meaningful.
- [ ] Report fixed-atom and anchor displacement.
- [ ] Report peptide closure and connectivity pass rates.
- [ ] Report bond, angle, planarity, Ramachandran, rotamer, chirality, and clash metrics.
- [ ] Report top-1 and oracle top-k quality.
- [ ] Report candidate diversity and ranking enrichment.
- [ ] Report success, unsupported, and failure rates by stratum.
- [ ] Report wall time, model-load time, peak RAM, and peak VRAM.
- [ ] Report external-process timeout and crash rates.
- [ ] Report conditioning/reinjection/boundary-refinement ablation results.
- [ ] Require experimental-CLI median gap-backbone RMSD to be no more than `0.25 Å` worse than MODELLER on the initial slice.
- [ ] Require junction pass rate to be no lower than MODELLER on the initial slice.
- [ ] Require stricter fixed-coordinate adherence than the MODELLER baseline.

## Hardware And CI Matrix

### CPU-Only Linux

- [ ] Run contract serialization tests.
- [ ] Run mask and identity tests.
- [ ] Run Kabsch and reinjection tests.
- [ ] Run fake-runner tests.
- [ ] Run provenance, checksum, and license-inventory tests.
- [ ] Run timeout, crash, containment, and atomic-publication tests.
- [ ] Run PDB validation, `.dat`, chirality, and benchmark-scoring tests.
- [ ] Treat any tiny CPU inference as an optional smoke test only when upstream officially supports it.
- [ ] Do not use CPU inference for representative ensembles or acceptance benchmarks.

### macOS Apple Silicon

- [ ] Run the same core, fake-runner, and validation suites as CPU Linux.
- [ ] Treat MPS inference as exploratory unless the selected pinned engine officially supports it.
- [ ] Do not use MPS output as a release acceptance gate.
- [ ] Document that Docker Desktop on macOS does not provide NVIDIA CUDA.

### Linux/NVIDIA

- [ ] Use Linux x86_64 with a pinned NVIDIA driver and CUDA/container runtime.
- [ ] Start short-gap evaluation on a machine with at least 24 GB VRAM.
- [ ] Prefer 48-80 GB VRAM for larger complexes, all-atom models, multi-sample ensembles, and profiling.
- [ ] Provide at least 64 GB system RAM and fast SSD/model-cache storage.
- [ ] Record actual peak memory instead of treating 24 GB as an upstream guarantee.
- [ ] Add a documented representative smoke profile before selecting long-term hardware.

### CI Separation

- [ ] Make the CPU contract lane mandatory.
- [ ] Keep the GPU inference lane opt-in or self-hosted until its environment is stable.
- [ ] Return a clear skip or preflight result when no GPU is available.
- [ ] Do not let missing CUDA fail the core test suite.

## Expected Files

- [x] Add `src/dvbfixer/model/diffusion/__init__.py`.
- [x] Add `src/dvbfixer/model/diffusion/contract.py`.
- [x] Add `src/dvbfixer/model/diffusion/masks.py`.
- [x] Add `src/dvbfixer/model/diffusion/geometry.py`.
- [x] Add `src/dvbfixer/model/diffusion/runner.py`.
- [x] Add `src/dvbfixer/model/diffusion/validate.py`.
- [ ] Add `src/dvbfixer/model/diffusion/pipeline.py`.
- [ ] Change `src/dvbfixer/model/cli.py` only in the experimental CLI phase.
- [ ] Change `src/dvbfixer/model/pipeline.py` only for explicit backend dispatch after gates pass.
- [ ] Change `src/dvbfixer/doctor.py` for runner, checkpoint, and device preflight.
- [ ] Change `src/dvbfixer/homology_plan.py` only in the mosaic-first phase.
- [ ] Change `src/dvbfixer/homology.py` only after the homology contract passes focused tests.
- [ ] Avoid changing `src/dvbfixer/ffutils/dat.py` unless a general sidecar requirement is demonstrated.
- [ ] Add isolated runner/container lock and license-inventory files.
- [x] Add focused diffusion unit and integration tests for the implemented contract, masks, geometry, research inventory, runner protocol, fake runner, and first-slice independent validation. Publication tests remain tied to the atomic-publication slice.
- [ ] Add reviewed fixtures and regenerate their manifest when benchmark structures are added.

## Focused Verification

- [x] Run the diffusion core tests:

  ```bash
  pytest -q tests/test_diffusion_contract.py tests/test_diffusion_geometry.py \
    tests/test_diffusion_masks.py tests/test_diffusion_research_inventory.py \
    tests/test_diffusion_runner.py tests/test_diffusion_validate.py
  ```

- [x] Run existing model and chirality regression tests:

  ```bash
  pytest -q tests/test_model_fasta_case.py tests/test_model_renumber.py \
    tests/test_model_strip_heterogens.py tests/test_ffutils_chirality.py
  ```

- [ ] Run homology tests when Phase 5 changes begin:

  ```bash
  pytest -q tests/test_homology_plan.py tests/test_homology.py \
    tests/test_homology_fasta.py
  ```

- [x] Run repository documentation and whitespace checks:

  ```bash
  python scripts/check_agent_docs.py
  git diff --check
  ```

- [ ] Run generated-reference updates and checks after argparse changes:

  ```bash
  python scripts/gen_cli_reference.py
  python scripts/gen_gui_spec.py
  python scripts/gen_cli_reference.py --check
  python scripts/gen_gui_spec.py --check
  ```

- [ ] Run a pinned GPU smoke test on one short internal gap before broader GPU benchmarks.
- [ ] Run a same-seed deterministic repeat on the pinned GPU environment.
- [ ] Run the representative benchmark subset before the full benchmark corpus.

## First Implementation Iteration Exit Criteria

- [ ] Complete Phase 0 and Phase 1 only.
- [ ] Keep the public CLI unchanged.
- [ ] Keep the MODELLER default path unchanged.
- [ ] Add no ML dependencies to the core environment.
- [ ] Provide a versioned backend-neutral contract.
- [ ] Provide deterministic mask construction.
- [ ] Provide weighted Kabsch and reinjection primitives.
- [ ] Provide a deterministic fake external runner.
- [ ] Provide validation and atomic publication.
- [ ] Provide a separate diffusion provenance manifest.
- [x] Provide a benchmark manifest with predeclared thresholds.
- [ ] Pass the CPU/Linux and macOS-compatible core suite before beginning the RFdiffusion GPU adapter.
