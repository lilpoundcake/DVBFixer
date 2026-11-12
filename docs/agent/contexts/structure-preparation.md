# Structure Preparation

Status: partial

Verified on: 2026-10-02

Verified at commit: `68d878f50d069268b60250cb162da2b9adf486cb`

## Purpose

Define the bounded context from sequence-guided residue rebuilding through
heavy-atom completion, protonation, hydrogen placement, and the `.dat` handoff
to minimization. It records ownership, ordering, and externally visible
contracts; it does not reproduce scientific algorithms.

## Scope

- `src/dvbfixer/model/` owns sequence/SEQRES placement, Modeller invocation,
  restoration of PDB identity/connectivity, and modeled-atom provenance.
- `src/dvbfixer/prepare/` owns the public preparation command and the default
  broad-input `legacy` path, including mutations, heterogens, and ligand SMILES.
- `src/dvbfixer/prep_backend.py` owns the opt-in `tleap-reduce` implementation.
- `src/dvbfixer/protonate.py` owns standalone post-hoc protonation. It is not a
  ZBS stage because legacy prepare already runs PROPKA and Reduce internally.
- `src/dvbfixer/ffutils/geometry.py` owns shared rebuild, hydrogen-geometry,
  disulfide-bond, clash, and C-alpha chirality guards.
- `src/dvbfixer/ffutils/dat.py` owns the sidecar schema and merge semantics.
- `src/dvbfixer/zbs.py` owns stage order, option propagation, artifact placement,
  final numbering, cleanup, and postflight, not scientific decisions.

## Capabilities

Alternate-state selection is implemented by `altloc_selection.py` with CLI
`altlocs.py`. All policies select a whole labeled state per residue and retain
shared atoms. The GUI Altlocs panel uses Python analysis and source-digest-bound
plans through `gui/server/altloc-api.ts`; Node never duplicates occupancy science.
Preview serial masks refer to the normalized PDB and explicit model ordinal.
Keep MODEL, chain case, insertion code, occurrence and connectivity identity
distinct. Apply writes a new structure; unresolved altlocs are rejected before
diffusion runner launch. Focused checks: `tests/test_altlocs.py`,
`gui/server/altloc-api.test.ts`, and `gui/src/lib/altloc-selection.test.ts`.

| Capability | Status | Owner | Boundary |
|---|---|---|---|
| Sequence/SEQRES gap rebuilding | implemented | `model/pipeline.py::main` | Modeller changes protein gaps; post-processing restores retained structure context |
| Broad-input preparation | implemented | `prepare/pipeline.py::run_pdbfixer` | Default `legacy`; supports glycans, ligands, PTMs, mutations, and covalent HETATM links |
| Pure-protein tleap/Reduce preparation | implemented | `prep_backend.py::run_prep` | Opt-in; unsupported chemistry fails rather than falling back automatically |
| Protonation decisions and H placement | implemented | `protonate.py::decide_protonation`, `prepare/pipeline.py::_run_propka_reduce_variants` | Decision evidence, topology H placement, and output naming remain separate concerns |
| Preparation provenance handoff | implemented | `ffutils/dat.py::DatRecord` | Model/homology produce; prepare merges; minimize consumes |
| Experimental diffusion gap reconstruction | partial | `model/diffusion_cli.py`, `model/diffusion/contract.py`, `scope.py`, `gap_capabilities.py`, `protein_context.py`, `adaptive_context.py`, `runner.py`, `validate.py`, `pipeline.py`, `provenance.py`, `deploy/protenix-v1/production_runner.py`, `deploy/protpardelle-1c/production_runner.py` | Explicit diffusion workflows publish independently validated atomic bundles with digest-linked sampler traces. Length evidence ranges remain profile-specific, but lengths outside them now warn and permit experimental attempts across CLI/mosaic/scope; no scientific acceptance is implied. Topology, chemistry, context, resource and output-integrity checks remain mandatory. Protenix supports explicit joint multichain/multi-gap context and mosaic edit windows. Richer cc89/adaptive context retains its separate research gates. MODELLER remains default. |
| Uniform hard L-chirality output gate | partial | `ffutils/geometry.py::assert_all_l` | Modeled candidates assert; prepare/protonate paths can only repair or warn before minimize |
| Private internal-coordinate sidechain kernel | research | `protein_completion/service.py::complete_sidechain`, `protein_completion/mixed_cohort.py::plan_mixed_cohort`, `protein_completion/mixed_evidence.py::read_mixed_evidence` | Standard/paired/anchored mixed joint/global-THR/rounded checks with distinct replay-verified manifest-last evidence for admission/raw/rounded failures. Old readers reject mixed kinds; production and empirical/native qualification stay gated |
| Matched own/seeded-template controls | research | `scripts/branched_template_baseline.py::template_proposal`, `docs/research/protein-completion-template-comparison.md` | 32 retrospective own-success-filtered controls: common passes 32 vs 15; template better on 5/6 mutually passing paired masks. Separate nearplanar gates/CIP inversions; not full legacy, broad superiority or production qualification |
| Rule-selected branched comparison | research | `scripts/unfiltered_branched_comparison.py::plan`, `docs/research/protein-completion-unfiltered.md` | 48 predeclared masks retain 16 source refusals; own/template pass 31/18 of 32. Native controls 29/32 expose VAL39 1–4 gate limitation; no threshold tuning. VAL swaps require both methyls generated; no fixed-partner relabelling |
| Rosetta repack mechanism | research | `scripts/research_rosetta_repack.py::probe`, `docs/research/rosetta-repack-mechanism.md` | Installed 2026.06: 16 cells show residue-level repack moves observed singleton anchors even with preserve-Cβ. Import packing is a separate stage. Count-pair cp_half=0.2, classic intra-rep differs from xover4. Not an atom-fixed backend or matched accuracy benchmark; no exported rotamers/dependency |
| Atom-fixed Rosetta scoring | research | `scripts/research_rosetta_fixed_scoring.py::run`, `docs/research/rosetta-fixed-scoring.md` | 128 retained inputs score without physical heavy motion after typed-virtual separation. Old gates unchanged: 31 own selections; three correlated THR masks improve retrospectively. No repack, no new sampling/public backend. Partial-graph diagnostics retain all old clashes; missing feasible scores block selection |
| Predeclared scoring follow-up | research | `scripts/prospective_branched_cohort.py::freeze`, `docs/research/rosetta-scoring-followup.md` | 1CRN/1AKI/1BTA: 48 masks, 16 H-scope refusals, own/template common passes 31/14 and 32 native passes. Explicit post-reconstruction CONECT-order amendment; same own proposals and scoring rule. All 128 score; one lysozyme VAL selection improves. Not an unamended/randomized/blind or broad qualification |
| Exact two-target scoring controls | research | `scripts/research_rosetta_joint.py::score`, `docs/research/rosetta-joint-controls.md` | 12 cases, 34 complete source-feasible combinations +12 native controls, 92 target-view imports and 34 canonical generated-only sidecars verified. Ten single-axis matrices, one empty and one numerically additive 3x3: no demonstrated coupled-packing benefit. No approximate graph/public integration; retain old guards |
| Input-backbone-selected close contacts | research | `scripts/research_rosetta_contacts.py::enroll`, `docs/research/rosetta-contact-controls.md` | CA-only rule retains zero 1CRN pairs and four 1AKI pairs/12 masks. 26 combinations, three real new–new joint-clash rejections, 23 canonical sidecars, 76 exact view imports. All matrices remain single-axis; independent winners avoid clashes. Native measured geometry 9/12, three correlated source-contact refusals; not generator admission or coupled-packing qualification |
| Private nine-chi contact controls | research | `scripts/research_rosetta_dense.py::prepare`, `docs/research/rosetta-dense-controls.md` | All frozen masks; 141 combinations, 12 joint rejections, 129 generated-only sidecars and 306 exact view imports. Seven two-axis matrices but ZERO changed coarse-grid score winners; 38 coarse/native poses/snapshots/scores repeat exactly. Two factorized disagreements use clashing references. Separate measured-chi/grid policy, not old/public admission or practical coupled-packing qualification |
| Exhaustive reference-sensitivity replay | research | `scripts/analyze_rosetta_references.py::analyze_matrix`, `docs/research/rosetta-reference-controls.md` | Retrospective saved scores only: all 141 anchors/381 rectangles, 255 fully feasible/nonadditive rectangles. All 104 feasible-term global comparisons and 129 supported-subset comparisons select the same minima; eight subsets exclude the global winner. Independent raw-cell/four-cell replay. No new sampling, favorable-reference choice, production route or demonstrated packing gain |
| Pre-download new-system controls | research | `scripts/prospective_independent_cohort.py::enroll`, `docs/research/rosetta-independent-controls.md` | 2CI2/1R69/1CTF, 44 masks (1CTF THR shortfall retained), 176 exact score snapshots, own/template common passes 41/16 and native 37. Six ranking gains are three correlated residues; matched template accuracy shows no own superiority. Two confirmed raw D-Cα and one CA-unassigned collapse after five seeded retries are rejected, not full-legacy outputs. No public/gate change |
| Template initialization causality | research | `scripts/research_template_causality.py::run_arm`, `docs/research/template-initialization-causality.md` | Four-case 2x2 factorial: 16 arms, 60 seeds, 300 full snapshots/75 target iterates; all 16 prior stock-normal seed states reproduce exactly. Primary D-Cα exists before minimize and without relaxation; proper whole-template L still yields fixed-backbone hybrid R. THR27 small-volume sign change is not confirmed D. No installed/production patch or qualified remedy |
| Proper-fit used-anchor ablation | research | `scripts/research_template_anchors.py::make_anchor_overlay`, `docs/research/template-anchor-controls.md` | Four-case mechanism factorial retained separately; actual supplied/used anchors differ explicitly. No production patch |
| Eleven whole-mask NCAC expansion | research | `scripts/audit_ncac_expansion.py::audit`, `docs/research/ncac-whole-mask.md` | All eleven retained, four seen/seven new only to variant. Stock/NCAC/own common passes 4/11/11; native 9/11. All 23 stock states/four seen NCAC traces repeat exactly. NCAC CA S throughout but accuracy mixed (median 0.323 vs own 0.294 Å); geometry published before native reporting, no blinded-family or public qualification |
| NCAC relative-rotamer selection | research | `scripts/audit_ncac_rotamers.py::audit`, `docs/research/ncac-relative-rotamers.md` | Frozen 0/+120/−120° distal-pair rotations, CB/observations fixed; 33 proposals, three geometry refusals, 30 audited fixed-heavy imports. Four changed THR choices all improve; median 0.323→0.226 Å. Retrospective, no empirical-library or independent-family qualification; next new enrollment, not tuning old angles |
| Prospective NCAC new-input controls | research | `scripts/audit_prospective_ncac.py::audit`, `docs/research/ncac-new-system.md` | Pre-download 1PGB/1BPI/1CSP: 11 masks, three shared refusals, four context/NCAC blockers; only four matched selected outputs. Thirty inspected scorer imports, six reject extra context atoms and observed PRO movement. NCAC scoring worsens median 0.439→1.223 Å; own baseline 1.172. Negative transfer result retained, not public readiness |
| Installed PyRosetta rotamer authorization | research only | `docs/plans/pdbfixer-reconstruction-handoff.md` | User authorizes installed-library local research; upstream archive is not a research blocker. No bundled-table redistribution or production dependency. Before production/public integration, explicitly request original tables and source/version/license permissions from the user, verify suitability and obtain production approval; scientific qualification remains separate |
| Installed-library prior-only ranking | research | `scripts/audit_dunbrack_prior.py::audit`, `docs/research/dunbrack-prior.md` | Twenty-two retained masks/fifteen queried pools, seven inherited nonexecutions. Active StpDwn_0-0-0 binary hashed; full API process repeat exact. Raw API weights are not normalized probabilities; explicit amended protocol retained. No new-input accuracy rescue, one old THR choice worsens; no table samples published or production qualification |
| Installed-library mean proposals | research | `scripts/audit_dunbrack_means.py::audit`, `docs/research/dunbrack-means.md` | Forty-five proposals/42 common passes and fixed-heavy imports, CB/observations fixed, seven inherited nonexecutions. Old median improves 0.226→0.198 Å, new-input scored median worsens 1.223→1.930; VAL29 score choice is wrong despite a good prior mean. No published table parameters or production qualification |
| Saved mean context/pair diagnosis | research | `scripts/audit_mean_comparability.py::audit`, `docs/research/mean-context-pairs.md` | All 84 saved mean/relative snapshots match context across all 15 executable pools; methyl symmetry does not remove VAL failures. All two prior/score disagreements have four exact-state replays and closed pair/onebody accounting; no changed non-target edges for the three audited pair terms. JSON-only tuple/list amendment retained, no numeric tolerance added. No weight refit or production qualification |
| Generated-only target-χ relaxation | research | `scripts/audit_generated_chi_relax.py::audit`, `docs/research/generated-chi-relax.md` | All 42 feasible mean starts refined, protected motion zero, CA S/THR CB R throughout; 41 common passes/one rejected clash. All 15 selected indices and rounded outputs unchanged. Fixed max200/tol0.01, no convergence claim; raw refined-H score semantics explicit, not fresh-import fixed-heavy scoring. No production qualification |
| χ stopping-policy controls | research | `scripts/audit_chi_stopping_controls.py::audit`, `docs/research/chi-stopping-controls.md` | Three predeclared sequential controls on all 42 starts each, 126 exact replays, zero protected motion, 40 common passes/mode. Stricter stopping changes one chosen index and up to all 15 rounded outputs; final two modes still differ. Native-free: no accuracy/profile choice or convergence certification, no new objective/weights/production route |

## Entry Points

| Change | Start symbol |
|---|---|
| Gap-model workflow | `src/dvbfixer/model/pipeline.py::main` |
| Experimental diffusion CLI | `src/dvbfixer/model/diffusion_cli.py::run_diffusion_model` |
| Diffusion contract | `src/dvbfixer/model/diffusion/contract.py::DiffusionRequest`, `DiffusionResult` |
| Partial-residue mask inventory (planning only; no sampling capability) | `src/dvbfixer/model/diffusion/atom_completion.py::plan_atom_completion`, `plan_request_atom_completion`, `audit_atom_completion_masks` |
| Private partial-residue protocol/insertion (not runner-enabled) | `src/dvbfixer/model/diffusion/atom_completion_contract.py::AtomCompletionRequest`, `require_completion_capability`; `src/dvbfixer/model/diffusion/atom_completion_materialize.py::materialize_atom_completion` |
| Private completion validation/refinement (not checkpoint qualification) | `src/dvbfixer/model/diffusion/atom_completion_validation.py::audit_completion_candidate`, `validate_completion_candidate`; `src/dvbfixer/model/diffusion/atom_completion_refinement.py::refine_completion_candidate` |
| Private cc94 completion sampling (not a public profile) | `src/dvbfixer/model/diffusion/atom_completion_sampling.py::assess_completion_sampling`, `completion_engine_view`; `deploy/protpardelle-1c/atom_completion_research.py` |
| Experimental completion trace / `.dat` / publication | `src/dvbfixer/model/diffusion/atom_completion_trace.py::CompletionSamplerTrace`; `atom_completion_evidence.py::audit_completion_evidence`; `atom_completion_publication.py::publish_completion_bundle`, `read_completion_bundle` |
| Private template comparison / held-out scoring | `src/dvbfixer/model/diffusion/atom_completion_benchmark.py::template_completion_proposal`, `heldout_completion_metrics` |
| Pre-diffusion observed-heavy preparation | `src/dvbfixer/model/diffusion/observed_completion.py::missing_observed_atoms`, `prepare_observed_completion`, `audit_completion_preparation` |
| Diffusion scope boundary | `src/dvbfixer/model/diffusion/scope.py::assess_diffusion_scope` |
| Diffusion pipeline | `src/dvbfixer/model/diffusion/pipeline.py::run_diffusion_pipeline` |
| Diffusion adapter preflight | `src/dvbfixer/model/diffusion/preflight.py::assess_adapter_preflight` |
| Diffusion provenance | `src/dvbfixer/model/diffusion/provenance.py::DiffusionProvenanceManifest` |
| Proposed diffusion CPU benchmarks | `src/dvbfixer/model/diffusion/benchmark.py::assess_same_seed_repeatability` |
| Proposed diffusion sampler conformance | `src/dvbfixer/model/diffusion/sampler.py::assess_sampler_conformance` |
| Diffusion runner boundary | `src/dvbfixer/model/diffusion/runner.py::run_diffusion_runner` |
| Internal RFdiffusion v1 benchmark adapter | `src/dvbfixer/model/diffusion/rfdiffusion_v1.py::run_adapter` |
| Diffusion validation boundary | `src/dvbfixer/model/diffusion/validate.py::validate_runner_result`, `build_validated_result` |
| Legacy atom completion | `src/dvbfixer/prepare/pipeline.py::run_pdbfixer` |
| Private own sidechain completion | `src/dvbfixer/protein_completion/service.py::complete_sidechain`, `src/dvbfixer/protein_completion/contracts.py::CompletionRequest` |
| Backend selection/final prepare writes | `src/dvbfixer/prepare/pipeline.py::main` |
| tleap/Reduce behavior | `src/dvbfixer/prep_backend.py::run_prep` |
| Standalone protonation | `src/dvbfixer/protonate.py::main` |
| Missing-atom rebuild policy | `src/dvbfixer/ffutils/geometry.py::rebuild_missing_atoms_with_retry` |
| Sidecar schema/merge policy | `src/dvbfixer/ffutils/dat.py::DatRecord` |
| Full pipeline order/propagation | `src/dvbfixer/zbs.py::_run_pipeline` |

## Contracts

- Scientific stages consume normalized PDB. CIF normalization and chain mapping
  belong at the CLI boundary; do not add stage-local CIF readers.
- PDB identity is case-sensitive. Atom identity is `(chain, resid, icode, atom)`;
  variant identity is `(chain, resid, icode)`. Preserve `keepIds=True` on every
  `PDBFile.writeFile` call.
- `DatRecord` is the only `.dat` reader/writer. `added_atoms` uses four-field
  atom identity; `variant_overrides` serializes as `chain:resid:icode`;
  `total_added` is derived. Merge deduplicates atoms and gives the downstream
  record precedence on variant collisions.
- Model marks every atom in a rebuilt residue as added. Prepare records its own
  additions and merges an adjacent upstream sidecar when found. Minimize
  re-resolves those identities against its final topology before restraint use.
- Diffusion is explicit opt-in and never falls back to MODELLER or another
  profile/device. Its `-o` destination is a new directory bundle; MODELLER keeps
  the existing file-output contract.
- Supplied ligand SMILES is authoritative for mapped, isolated, single-residue
  molecules. It is legacy-only, incompatible with heterogen stripping and
  disabled heterogen-H placement, and must fail rather than guess on an
  incompatible, ambiguous, or covalently attached target.

## Invariants

- Legacy PDBFixer order is fixed: `findMissingResidues` plus deletion scrub,
  `findNonstandardResidues`, `replaceNonstandardResidues`, `removeHeterogens`,
  `findMissingAtoms`, then rebuild. Topology-rebuilding calls invalidate
  residue-object-keyed missing-atom data.
- Never call `fixer.addMissingAtoms()` directly. Call
  `rebuild_missing_atoms_with_retry`; it owns explicit seeds, retries, topology
  snapshots, and post-attempt chirality/clash checks.
- Legacy PROPKA/Reduce runs only after heavy-atom rebuilding and chirality
  repair. Hydrogen placement uses `Modeller.addHydrogens(variants=...)`, never
  `PDBFixer.addMissingHydrogens`.
- Follow every legacy `Modeller.addHydrogens` call with
  `repair_misplaced_hydrogens`. Preserve the shared variant-parent rename and
  post-hydrogen restoration flow.
- Shared H postflight detects O–H at/above 1.2 Å without tightening C/N/S/P
  thresholds, and O-bound H within/at 1.0 Å of a same-residue nonparent heavy
  atom. Safe moderate stretches retain direction with nonparent clearance;
  valid H and heavy coordinates stay fixed. An unresolved O–H violation raises
  `HydrogenGeometryError`; prepare must not resample H via template fallback.
  Scope, qualification and PDB-rounding limitations are in the
  [implementation record](../../plans/hydrogen-postflight-contract.md).
- Explicit disulfide evidence forces CYX and genuine SG-SG bonds survive bond
  cleanup; extra SG-SG partners are reduced to one partner per sulfur.
- Any stage claiming a hard chirality gate must run `assert_all_l` after the
  final heavy-atom coordinate change. ZBS catches `ChiralityError`; minimize is
  the unconditional final safety net for warn-only preparation paths.
- Model performs CONECT inference before Modeller, uses shared affine sequence
  placement, restores chain/residue identity, makes atom serials unique, then
  remaps retained CONECT records. Bare `TER` records must remain newline-ended.
- `--number-from-1` applies to each model candidate and matching sidecar, but in
  ZBS only to the copied final output before postflight. Never renumber an
  intermediate sidecar independently of its PDB.

## Callers

- The command dispatcher invokes `model`, `prepare`, `protonate`, and `zbs`.
- ZBS orders `renumber -> model -> prepare -> minimize`; any stage may be
  skipped. It forwards sequence, heterogen, CONECT, force-field, backend,
  protonation, mutation, cap, SMILES, and minimization options to their owners.
- `minimize/pipeline.py` consumes the prepared PDB and merged sidecar. Homology
  can publish a model plus stub sidecar before invoking preparation/minimization.

## Adapters

- Modeller supplies sequence-guided coordinates. PDBFixer supplies legacy
  topology repair. OpenMM Modeller supplies legacy protein H placement.
- tleap supplies opt-in heavy-atom completion; Reduce supplies H placement and
  local HIS/ASN/GLN evidence; PROPKA supplies pKa evidence. Python owns evidence
  precedence and identity preservation.
- RDKit/Open Babel handle unmapped legacy heterogens. SMILES graph matching is a
  distinct authoritative adapter and does not inherit the requested protein pH.
- PDB fixed-column rewrites, subprocess execution, temporary files, and `.dat`
  JSON are adapters around the scientific stage contracts.

## Side Effects

- Model may materialize an inferred-CONECT input, creates a temporary workspace,
  and writes one or more candidate PDBs with candidate-matched sidecars when
  provenance exists. The no-SEQRES/no-FASTA shortcut only copies the PDB.
- The diffusion scope boundary verifies the normalized-input digest
  before parsing and returns explicit unsupported reasons for out-of-slice target
  sequence, gap, model, altloc, chemistry, heavy fixed-mask, placement, retained
  source-link, and detectable D-chirality cases. The runner creates an isolated
  workspace, copies and verifies the normalized input, writes `request.json`,
  launches one shell-free child process with bounded diagnostics, then verifies a
  raw `RunnerResult` manifest and candidate digests. Runner-provided data cannot
  claim independent validation. Validation reopens verified source/candidate
  artifacts, applies first-slice hard gates, calls final `assert_all_l`, rejects
  failed candidates, and deterministically ranks passing candidates without
  treating backend-native scores as commensurate. The pipeline launches
  no runner for unsupported scope, publishes nothing for failed results, builds
  each candidate `.dat` through `DatRecord`, writes sanitized separate provenance,
  fsyncs a hidden same-parent bundle, and commits the complete set with an
  atomic no-replace directory rename where available. On POSIX filesystems
  without rename flags, `model/diffusion/bundle_commit.py` exclusively reserves
  the destination with mkdir, transfers/fsyncs payloads, and atomically links
  `bundle.json` last without replacement. That index, not directory existence,
  establishes readiness; crashed reservations without it are rejected and never
  automatically adopted. Completed bundles retain schema 2 and single-linked
  regular artifacts. The public single-input `model` dispatcher
  builds the request and requires an explicit profile, runner, and local
  checkpoint. Local adapter preflight checks executable,
  immutable image identity, checkpoint presence/digest, CUDA availability, and
  protocol compatibility without downloading artifacts. Sampler conformance
  keeps ordinary template conditioning distinct from evidence-backed per-step
  reinjection and boundary refinement. Batch, GUI, ZBS, and Homology do not
  expose this path.
- Legacy prepare can create temporary inferred, renamed, deletion-cleaned,
  capped, GLYCAM, PROPKA, and canonical-CONECT PDBs. It writes/replaces the
  output PDB and sidecar and may rewrite the output for variants, heterogen H,
  CONECT, force-field naming, and terminal normalization.
- The tleap backend creates an isolated temporary directory, launches `tleap`
  and `reduce`, restores residue metadata, and includes subprocess diagnostics
  in failures. It writes a prepared PDB and sidecar.
- Protonate launches configured decision engines, writes one PDB, and creates
  sanitized/intermediate temporary PDBs. It does not publish a `.dat` sidecar.
- ZBS places named intermediates beside the final output, optionally aligns each
  stage to the original input, deletes tracked intermediates unless retained,
  and writes a diagnose JSON report unless postflight is disabled.

## Known Divergences

- `legacy` is the default and broad-input backend. `tleap-reduce` is opt-in,
  pure-protein-oriented, rejects `--mutate` and SMILES, and bypasses final
  `apply_variants_to_pdb_text` shifts that would break OpenMM terminal templates.
- Prepare's tleap wrapper always runs Reduce and variant assignment; it does not
  currently honor `--no-propka` or `--no-protassign`. Standalone tleap protonate
  uses `--no-propka` only to disable variant assignment and still runs Reduce.
- Legacy prepare degrades PROPKA/Reduce failures to reduced evidence. Standalone
  legacy protonate requires at least one engine and treats missing default
  Reduce, empty required PROPKA results, and H-placement failures as fatal.
- Model's no-gap sidecar emits legacy two-field variant keys, and model/protonate
  still contain local two-field variant restoration maps. These are known gaps,
  not alternatives to the three-field `DatRecord` identity contract.
- Model's rebuilt-candidate path hard-asserts L chirality. Legacy prepare repairs
  before H placement without a final assertion; tleap prepare and both protonate
  paths warn on residual D centers. Standalone outputs therefore have a weaker
  guarantee than a completed ZBS run through minimize.
- tleap metadata restoration is residue-ordinal-based and leaves unmatched extra
  output residues unchanged; supported pure-protein input is expected to keep
  residue order and count stable.

## Proposed Work

- Make protonation flags consistent across backends and migrate remaining local
  variant maps to insertion-code-aware identity.
- Add uniform final `assert_all_l` gates before changing the page status to
  `implemented`.
- No geometry-regularization backend is shipped. The
  [`whole-complex relaxation note`](../../research/whole-complex-relaxation.md)
  is research, not current preparation behavior.
- Alternative atom-reconstruction, loop-modeling, and homology-modeling
  backends remain research or proposed work. The accepted
  [`experimental diffusion boundary ADR`](../../adr/0010-experimental-diffusion-gap-reconstruction-boundary.md)
  and its [`public CLI successor`](../../adr/0011-experimental-diffusion-cli.md)
  fix the isolation, provenance, no-fallback, default-backend, narrow-scope, and
  profile-evidence policies. See
  the [`reconstruction and modeling backend note`](../../research/reconstruction-and-modeling-backends.md),
  its versioned
  [`diffusion research inventory`](../../research/diffusion-gap-reconstruction-inventory.toml),
  and the active
  [`diffusion gap-reconstruction implementation plan`](../../plans/diffusion-gap-reconstruction.md).
  PDBFixer and Salilab MODELLER remain supported production baselines; the
  experimental diffusion CLI does not deprecate either dependency.

## Focused Verification

```bash
pytest -q tests/test_ffutils_dat.py tests/test_ffutils_chirality.py \
  tests/test_seeded_atom_rebuild.py tests/test_model_renumber.py \
  tests/test_model_fasta_case.py tests/test_model_strip_heterogens.py \
  tests/test_prep_backend_variants.py tests/test_prepare_icode_variants.py \
  tests/test_prepare_input_file_safety.py tests/test_prepare_propka_integration.py \
  tests/test_prepare_smiles.py tests/test_ffutils_variants.py \
  tests/test_zbs_postflight.py
python scripts/check_agent_docs.py
git diff --check -- docs/agent/contexts/structure-preparation.md
```

External-tool and end-to-end coverage lives in
`tests/test_prepare_by_input_class.py`, `tests/test_zbs_e2e.py`, and
`tests/test_zbs_shit_inputs.py`. Use `--no-solvent` for development ZBS runs.
