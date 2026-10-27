# Geometry-gated anchor relaxation: headless handoff

Prepared: 2026-10-27. Status: core implementation in progress; production integration pending.
Working branch: `research/9cqd-anchor-relaxation`.

This is an execution handoff, not a claim that the feature is shipped. Track
completion in the [canonical checklist](implementation-checklist.md#geometry-gated-anchor-relaxation).
See the [B-factor study](../research/9cqd-bfactor-anchor-restraints.md) and
[earlier factorial experiment](../research/cc89-lrd-ablation-9cqd.md) for evidence.

## Execution progress

Started 2026-10-27 on the research branch. Implemented initial pure
`anchor_policy.py` (provenance/fallback, bounded chain-B weights, explicit
per-anchor movement bounds with no numeric defaults), `refinement_stages.py`
(typed geometry/admission outcomes and default-enabled raw-based escalation),
and an explicit `observed_anchor_restraints` argument in the scientific refiner.
The generated-only default of the existing low-level refiner is preserved.

The staged controller uses injected independent assessment/minimization functions;
it does not authorize publication of observed-atom movement. Its stage-2 results
explicitly require a separate anchor acceptance contract. Production finding
adapters, versioned acceptance, runner/CLI/GUI wiring, final movement notices and
calibrated limits remain unfinished. Do not treat these helpers as a shipped
default change.

Pure-policy/controller, force/units, and existing runner/validator/contract/pipeline
tests pass (108 in the updated narrow set). The full fast-lane snapshot passed
1506 tests with two optional-dependency skips; subsequent movement-bound tests
were covered in the updated narrow set. The lightweight CI selection passes
391 tests with 61 expected skips while OpenMM/PDBFixer imports are deliberately
blocked. The two new pure test files are included in that CI lane.
The experimental minimum of eight complete chain
backbone samples is an engineering pilot setting, not an accepted scientific
threshold. The new formula has been replayed with uniform controls and automatic
selection of problematic anchors using
`scripts/run_cc89_staged_anchor_pilot.py`. All four outputs pass tested geometry
but still fail the historical source-drift contract. Off selects W29 only;
on selects W25/W29. X anchors and every other unselected source record remain
unchanged. Every attempt is bound to the original raw digest. Detailed results
are in the B-factor study. All outputs remain research-only; the new B formula
has no consistent advantage over uniform restraints.

## Start here

1. Read repository `AGENTS.md`, inspect git status and current branch.
2. Preserve existing user/research changes. This session's research scripts,
   results and documentation may still be uncommitted. Do not assume a clean
   checkout or replace it with main and lose this work.
3. Continue on a dedicated branch based on main; never implement directly on main.
4. Inspect current code and deployed profile/source hashes before editing.
5. Implement and verify the requirements below. This document does not authorize
   a main merge, release/version bump or publication of ignored local directories.

## User-approved behavior

- Geometry determines whether any repair is needed. A geometrically acceptable
  candidate receives no coordinate-changing relaxation.
- Stage 1 relaxes generated atoms only; observed heavy atoms remain fixed.
- If relevant geometry remains defective, stage 2 relaxes generated atoms plus
  explicitly selected problematic immediate observed anchors.
- **Stage 2 starts from the original raw generation candidate, never stage-1
  coordinates.** Only diagnostics and mask selection may flow from stage 1.
- **Stage 2 is enabled by default.** The user can disable it; proposed CLI spelling
  is `--no-relax-anchors`. Do not accidentally implement an opt-in default.
- B factors, when suitable, determine positional restraint strength, not whether
  there is a geometric defect or which observed atoms count as generated.
- Emit visible notifications identifying every relaxed anchor and its movement.
- Do not weaken scientific geometry gates or silently change source identities.

## Relevant code and artifacts

Scientific entry points to inspect:

- `src/dvbfixer/model/diffusion/boundary_refinement.py`: existing v8 minimizer.
- `src/dvbfixer/model/diffusion/contract.py`, `validate.py`, `quality.py`,
  `scope.py`, `protein_context.py`: contracts, geometry and context identities.
- `src/dvbfixer/model/diffusion/pipeline.py`, `runner_refinement.py`,
  `provenance.py`, `publication.py`: orchestration, custody and publication.
- `src/dvbfixer/model/diffusion_cli.py`, model CLI, ZBS, Homology and GUI routing.
- `deploy/protpardelle-1c/` and other runners that delegate refinement: profile
  revisions, pinned hashes and replay contracts must remain internally consistent.

Existing research-only code, not production integration:

- [Restraint experiment](../../scripts/cc89_anchor_restraints.py).
- [Narrow tests](../../tests/test_cc89_anchor_restraints.py).
- [Visual exporter](../../scripts/export_cc89_anchor_review.py).
- [Candidate-bound summary](../research/9cqd-anchor-restraints-seed7-v1.tsv).

Ignored local evidence, if present:

```text
test/errors_diffusion/9cqd-lrd-factorial-seed7-fixed-frame/
test/errors_diffusion/9cqd-anchor-restraints-seed7-v1/
test/errors_diffusion/9cqd-anchor-restraints-visual-review-v1/
```

Tracked source: `tests/fixtures/protein_context/9cqd_split.pdb`, SHA-256
`6bcc9638c3d299ea12151c3894b7e6c3482ebd2fdfe1d8e54dcc1f6f6316a281`.
Ordinary CI must not depend on ignored local results. Do not add generated PDBs
to source fixtures unless explicitly needed as reviewed golden outputs.

## Evidence versus assumptions

9CQD is a 3.10 Å X-ray structure. Source anchor backbone B means are W25 88.725,
W29 89.9875, X139 97.325 and X142 94.0025 Å². Medians over residue-mean backbone
B are W 74.7925 and X 67.4175 Å².

Twelve paired seed-7 minimizations were completed: reinjection off/on crossed
with fixed, free, uniform 1000, B-weighted 1000, uniform 10000 and B-weighted
10000. Both fixed controls failed geometry; all ten mobile-anchor candidates
passed tested geometry but failed the existing exact-source contract.
Uniform k=10000 yielded combined W-aligned anchor backbone RMSD 0.2379/0.2324 Å,
maximum original-frame anchor heavy displacement 0.8240/0.8838 Å. B weighting
showed no consistent benefit. These are not missing-loop accuracy measurements.

The proposed chain-normalized B formula below is different from the experiment's
equal-mean-over-anchor-atoms normalization. **It is untested.** One source/seed
does not establish general B cutoffs, safe drift limits or experimental accuracy.

## Required state machine

Retain immutable raw bytes/digest and source bindings for every candidate.

1. Classify independent findings. Identity/completeness, unsupported chemistry,
   metadata and other nongeometric failures must not be disguised as repairable
   geometry errors. Always retain final full validation.
2. If no actionable generated/junction geometry error exists, preserve coordinates
   and proceed through the normal independent publication checks.
3. Run stage 1 from raw with generated-only mobility. Assess the resulting geometry.
4. If stage 1 passes acceptance, select it and stop.
5. If residual findings justify specific immediate anchors and stage 2 is enabled,
   run stage 2 independently from RAW with the enlarged explicit mask.
6. Publish only a candidate passing the applicable complete acceptance contract.
   If stage 2 fails geometry or movement bounds, retain evidence and reject.

For a stage-1 numerical exception, retain a typed failure. If raw diagnostics
justify a repair mask, stage 2 may use that evidence; otherwise fail closed.
Do not release all anchors just because stage 1 crashed. An execution failure,
a geometry failure and a validation failure are distinct report outcomes.

## Problematic anchor selection

Reuse existing geometry detectors and numerical thresholds. Findings must expose
affected atom/residue identities, including both partners of junctions/clashes.

- Junction connectivity/angle/planarity defects identify the corresponding flank.
- Relevant anchor backbone/Ramachandran defects or severe generated-to-immediate-
  anchor overlaps can justify releasing that anchor.
- Internal generated-only defects alone do not justify observed-atom motion.
- Do not automatically release distant observed clash partners or both ends of
  every gap. Inventory unrelated pre-existing source defects separately.
- Preserve the existing chirality invariant; never relax an observed D centre as
  a way to bypass a source-admission failure.
- Distinguish actionable errors from warnings, preserving amide-planarity's actual
  hard-gate semantics even if a diagnostic finding has warning severity.

Use case-sensitive `(chain, resid, icode)` identities and atom names. Terminal
gaps have only their existing flank. Deduplicate shared anchors. Explicitly resolve
overlapping repair masks/ownership before combining coordinates; no independent
last-writer-wins repairs of a shared anchor. Joint diffusion generation remains
deferred and is not part of this task.

## Scientific restraint API

Replace research monkeypatching with an explicit minimizer API distinguishing:

- original generated identities;
- authorized mobile observed anchors/atoms;
- restraint centres, coefficients and their source provenance.

Re-resolve atom indices against the final OpenMM topology after atom/H rebuilding.
Keep all unselected observed heavy atoms fixed. Do not silently release an entire
residue when only particular atoms were authorized. Preserve original generated
mask semantics in validation and provenance.

The restraint is harmonic, not a hard distance constraint:

`E_i = 0.5 * k_i * |x_i - x_source_i|^2`

Centres are the original experimental coordinates, never a stage-1 result or an
earlier relaxed position. OpenMM coordinates are nm; k is kJ mol⁻¹ nm⁻².

### Proposed B weighting to test

```text
B_r   = mean(original observed B_N, B_CA, B_C, B_O for anchor residue r)
B_ref = median(per-residue observed backbone mean B in the same source chain)
k_r   = k0 * clip(B_ref / B_r, 0.5, 2.0)
k0    = 10000 kJ mol^-1 nm^-2 as the initial candidate setting
```

Apply one residue coefficient to its authorized observed heavy atoms. Generated
atoms have no source-coordinate restraint. Uniform fallback is k0. Calculate chain
reference distributions before mask selection; do not renormalize only over
released anchors, since a single released anchor would lose all B dependence.

This is bounded regularization, not an ADP-derived coordinate uncertainty model.
No universal absolute B threshold, no direct B/(8*pi²)-to-drift-radius conversion,
and no thermodynamic force-constant estimate from diffraction temperature.
Harmonic restraints alone do not guarantee a maximum displacement.

### Provenance and fallback

A numeric B column is insufficient. Distinguish experimental ADPs from pLDDT,
other confidence scores and generated placeholders; use only original observed
atoms in reference statistics. Do not infer provenance from numerical magnitude.

Define and test deterministic handling of missing method/provenance, nonfinite
or nonpositive values, incomplete backbone, alternate states, constant placeholder
distributions, small chains and mixed valid/invalid anchors. Use uniform fallback
with an explicit reason where B evidence is unusable. Exact minimum-population
criteria are not yet approved. Provide an explicit way to disable/override B
interpretation without conflating it with disabling anchor motion.

Per-chain statistics can mix domains/TLS motion. Document that limitation rather
than inventing domain segmentation in this first implementation.

## Acceptance and default-change safety

Default-enabled anchor relaxation changes coordinate-preservation semantics.
Introduce an explicit versioned mobility policy and validation/provenance record;
do not silently weaken the historical strict contract for every candidate.

- Keep original generated atoms and selected observed anchors distinct.
- Unselected observed records remain exact, allowing only the existing internal
  platform tolerance before source-preserving PDB materialization.
- Add per-anchor backbone RMSD AND maximum heavy-displacement limits relative to
  source coordinates. Do not dilute movement over the entire complex.
- **Numeric drift bounds require calibration and are not approved yet.** Do not
  derive them from B or the largest successful 9CQD displacement alone.
- Check generated regions, released anchors, outer junctions and clashes against
  the rest of the complex; perform final cross-gap checks after merging.
- Preserve all-L, atom identity/completeness, retained-link and metadata gates.
- Kabsch metrics are explanatory; original-frame movement governs preservation.
  No independent chain fits that conceal interchain movement.
- Source B values remain inherited metadata, not newly calculated uncertainty.

If movement bounds cannot yet be justified, report that as a production
qualification blocker. Do not ship unbounded default source motion or silently
replace the user's default-enabled decision with an opt-in design.

## Notifications and workflow integration

Thread the default-enabled policy and negative option through diffusion model,
ZBS, batch/resume, Homology and GUI without changing unrelated MODELLER or prep
workflows. Include settings in request digests and resume compatibility.

On escalation emit a visible WARNING with the reason and exact selected anchor
IDs. After selecting stage 2, emit a visible WARNING listing relaxed anchors,
actual backbone RMSD/max heavy movement and effective restraint policy. Distinguish
attempted relaxation from the selected output's actual coordinate changes.
Preserve fd-level CLI capture/diagnostic summaries; GUI notice must not be hidden
only in raw logs. Clean unchanged runs should not emit an empty warning banner.

Structured reports must retain raw digest, selected stage, each attempt's input
binding, trigger findings, anchor identities, B/reference/fallback, effective k,
source bindings, movement metrics and independent gate outcomes. Keep bounded
runner reports and historical compatibility; do not bypass profile digest checks.

## Work packages

1. Inspect current implementation, working changes, runner boundaries and CI.
2. Add pure B-policy, finding-classification and identity-mask helpers/tests.
3. Extend the scientific refiner for authorized anchor mobility/restraints.
4. Add raw-based two-stage orchestration and independent acceptance policies.
5. Update versioned contracts, provenance, runner/profile custody and resume.
6. Add CLI/GUI routing, defaults, opt-out and visible notifications.
7. Compare the NEW chain-B formula with uniform on paired raw candidates.
8. Expand the cohort and calibrate movement bounds before production qualification.
9. Update user docs, scientific limitations and generated help/schema artifacts.

## Verification requirements

Unit/integration coverage must prove:

- Clean geometry means no refiner call or coordinate change.
- Stage-1 success prevents stage 2; default/negative-option routing is correct.
- Stage 2 gets byte/digest-identical raw input, not stage-1 coordinates.
- Only justified immediate anchors move; terminal/shared/overlapping masks,
  insertion codes and case-distinct chains retain correct identity.
- B means/reference/clipping, fallback/provenance and single-anchor weighting work.
- Correct units: k=1000 and displacement 1 Å gives 5 kJ/mol positional energy.
- Topology rebuilds cannot shift restraint selection; centres remain original.
- Outside-mask records are preserved, all-L holds and excessive drift or new outer
  defects reject publication. Failure/cancellation evidence is retained.
- Resume/digests reject incompatible mobility policies; old strict requests are
  not silently accepted under a new source-motion exception.

Start scientific replay with 9CQD off/on raw candidates, clean-geometry controls,
fixed/free/uniform and new chain-B arms. Preserve old evidence; use new directories.
Broaden to low-B/high-confidence, absent/invalid-B, different resolutions, terminal
gaps and coupled repair neighborhoods. Never select seeds by reference RMSD.
Record geometry, source motion and experimental agreement separately. If density
is not evaluated, say so. No cross-platform bit-identical coordinate claims.

Always use `--no-solvent` on dev minimize/ZBS commands. For direct refiner calls,
create no solvent. Run narrow tests first, then the appropriate fast/full lane
when available. Follow AGENTS.md for Ruff, strict mypy, agent docs, GUI tests,
fixtures/provenance/checksums and scientific invariants. Lightweight CI must not
import OpenMM merely to inspect pure policy/registry labels.

Regenerate CLI reference after argparse changes; never hand-edit generated
reference pages. Regenerate/check GUI schemas and routing outputs as required.

## Completion report

Tell the user what is implemented versus research-only, exact default/opt-out,
geometry triggers, B policy/fallback, calibrated bounds, test results and blockers.
Provide paired output/motif paths, original-frame and shared-frame aligned views,
per-anchor metrics and failed-attempt evidence. A geometry pass is neither
experimental accuracy nor exact source preservation.
