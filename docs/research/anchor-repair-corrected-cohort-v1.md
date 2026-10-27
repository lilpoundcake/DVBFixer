# Geometry-gated anchor repair: corrected cc89 replay cohort v1

Date: 2026-10-27. Status: protocol frozen before minimizations; all eight cases completed.
Research only; not a production qualification or a calibration of safe movement.

## Scope and frozen cases

Replay ALL eight raw candidates in the existing corrected-Apple pilot, without
new neural sampling or selecting candidates by reference accuracy:

| Case | Seeds | Gap length |
|---|---|---:|
| prospective-8ek6-H-cdr1-gap5 | 7, 19 | 5 |
| prospective-9ww8-gap10 | 7, 19 | 10 |
| prospective-9ww8-gap5 | 7, 19 | 5 |
| prospective-9z0u-H-cdr3-gap10 | 7, 19 | 10 |

Inputs reside under `.artifacts/apple-cc89/pilot-v1/` and are optional local
evidence, not CI inputs. Original structures are tracked under
`tests/fixtures/cc89_apple_handoff/`. This is a convenience regression cohort of
two Fab structures and the bacterial transcriptional regulator HosA (9WW8), not an independent general
protein population. Keep 9CQD and native 1UBQ controls separate in reporting.

## Purpose

Test whether geometry-triggered generated-only repair avoids unnecessary observed
motion on previously successful corrected cases; identify residual defects that
actually require anchor motion. Preserve all failures instead of choosing easier
cases. The historical success of this pilot is known; this is not a prospective
blind efficacy claim. Raw and previous refined evidence are never overwritten.

Three labelled B policies share the same raw candidate and generated-only stage:

1. `auto`: trust the normalized source method only. Missing experimental headers
   must trigger uniform fallback if stage 2 is necessary.
2. `disabled`: uniform positional restraints, independently of B.
3. `experimental`: explicit research provenance assertion, permitted only after
   verifying X-ray/neutron provenance and matching every normalized
   observed atom's name/resname/coordinate/occupancy/B record to that original.

Preflight found the tracked protein-only PDB copies ALSO lack experimental
headers. Before any minimization, method provenance is therefore supplied by
frozen RCSB entry JSON with matching entry IDs (8EK6, 9WW8, 9Z0U), fetched
2026-10-27 from `https://data.rcsb.org/rest/v1/core/entry/{id}`. All three report
X-ray diffraction; reported resolutions are 1.954, 2.03 and 1.87 Å respectively.
The cohort saves those exact JSON bytes/digests. No PDB header is fabricated and
no B value is reinterpreted without the separate observed-atom content check.

Normalized pilot PDBs can lack EXPDTA even though their input entries are
experimental. Do not infer method from B magnitude or fabricate a header to enable
weighting. Save the verified upstream header and digest in a provenance sidecar.
Production propagation of this information is a separate unfinished task.

## Frozen execution and evaluation

- Stage 2 enabled by default in the controller, from raw, only for named residual
  problematic immediate anchors. Generated masks are never redefined.
- CPU one thread, no solvent, request seed unchanged.
- Match existing runner settings: gap length <=5 uses 500 iterations/attempt,
  0.25 Å perturbation; longer gaps use 250 iterations and 0.75 Å. Up to eight
  seeded restarts, retaining the existing early-exit behavior.
- Base positional k=10000 kJ mol⁻¹ nm⁻², bounded chain-B weighting as in the
  [handoff](../plans/geometry-gated-anchor-relaxation-handoff.md).
- Cache the identical generated-only minimization across policy labels. If two
  policies produce identical atom/centre/k force parameters, reuse that exact
  minimization too. Record reuse and actual calls; labels are not independent
  stochastic replicates.
- Reconcile only dependent PDB metadata of archived raw candidates to current
  request semantics before testing; assert all atom records unchanged. Retain
  original and reconciled raw hashes. This is not a coordinate repair.
- Independently assess geometry and nonmobile records. Save raw, every attempt,
  diagnostics, effective coefficients, source/raw bindings and launch snapshots.
- Audit stage-2 candidates using explicit research sensitivity bounds
  backbone RMSD <=0.5 Å and maximum heavy displacement <=1.0 Å per anchor. These
  limits are NOT adopted defaults or justified by this replay's acceptance rate.
- Preserve historical strict results alongside any research mobility audit.
  All reports remain unqualified for production publication.
- Track clean-raw, stage-1 success, stage-2 attempts, geometry failures, source
  admission problems, fallback reasons and individual anchor displacements.

No thresholds, masks, seeds or settings will be tuned after observing outcomes.
If no additional cases reach stage 2, report that the cohort provides no new
anchor-drift calibration evidence rather than presenting generated-only passes
as validation of the B formula or movement bounds.

## Results

All eight cases complete with zero execution errors. Every selected output
passes historical strict validation and expanded geometry checks. All observed
atom records remain exact; no anchor is released or moved.

| Case | Seed | Selected stage | Initial blocking geometry | Native refiner calls |
|---|---:|---|---|---:|
| 8EK6 H CDR1 gap5 | 7 | Raw | None | 0 |
| 8EK6 H CDR1 gap5 | 19 | Generated-only | Severe overlap | 1 |
| 9WW8 gap10 | 7 | Raw | None | 0 |
| 9WW8 gap10 | 19 | Raw | None | 0 |
| 9WW8 gap5 | 7 | Raw | None | 0 |
| 9WW8 gap5 | 19 | Raw | None | 0 |
| 9Z0U H CDR3 gap10 | 7 | Generated-only | Severe overlap | 1 |
| 9Z0U H CDR3 gap10 | 19 | Generated-only | Severe overlap | 1 |

Thus 5/8 raw candidates avoid unnecessary minimization, 3/8 are repaired with
generated-only mobility, and 0/8 require stage 2. There are THREE native refiner
invocations, each including its configured internal seeded retries, not 24
independent trials. The three B-policy labels share the same generated-only
result, and B weighting is never consulted in these cases. The JSON field
`unique_native_minimizations` counts refiner invocations, not individual internal
`minimizeEnergy` calls.
No missing-B warning or anchor-motion warning is needed when stage 2 is unused.

This replay supports conditional repair/no-motion behavior on a bacterial
regulator and two Fab structures. It contributes **no new measurements of anchor
drift** and cannot qualify the B formula or the 0.5/1.0 Å sensitivity limits.
Additional problematic-anchor real-gap cases, rather than more already-correct
raw candidates, are needed to calibrate movement bounds. No reference-loop RMSD
was used for selection or claimed as an accuracy outcome.

The archived raw files remain unchanged. Their dependent metadata is reconciled
in separate copies; all atom records were checked unchanged before any geometry
decision. Raw selections preserve those reconciled coordinate records exactly.
Original and reconciled digests are retained in each report.

## Evidence and checks

[Replay script](../../scripts/run_corrected_anchor_cohort.py) and
[provenance/cache tests](../../tests/test_corrected_anchor_cohort.py). Evidence:
`test/errors_diffusion/corrected-anchor-cohort-v1/` (ignored), including raw and
stage outputs, full assessments, source/entry metadata, force parameters, actual
native-call counts and launch-time code/protocol snapshots. All reports retain
`publication_qualified=false`.

| File | SHA-256 |
|---|---|
| manifest.json | a4f204f59b200e0ae13f25ac2f06331fcb6d5a306dda9101bb1b47369cef1d5c |
| comparison.json | f02430d750d3c7c877136c20ccddef77ce463b25533f8108686c46d186504fa1 |
| completed.json | 8e14764c5a64560f3a393051d642f696a2fc4cd85de4494864c204bec67d8b3d |

RCSB metadata fetched before execution and saved unchanged:

| Entry | JSON SHA-256 |
|---|---|
| 8EK6 | 75c9bcfc14d4587907339b0bd0d53282bca8aa71e849d6e72b824901e5608b2d |
| 9WW8 | bdd04c0ba11b17f4a1ea9cd371a0a2e37230e7f20c901e1abb4a5d235885bfcd |
| 9Z0U | a506c6980d9201be284684ad0e985fd366e396407704172a2b25a263b47036c1 |

The tracked protein copies also lack EXPDTA. Each normalized observed atom was
matched to its tracked upstream copy by name, residue name, coordinates,
occupancy and B, independently of chain/residue renumbering. Matching entry IDs
in the saved RCSB records establish method provenance; B magnitude is never used
to infer experiment type. The production request/normalization path still needs
to propagate equivalent provenance or explicitly use uniform fallback.

The broader real-gap/anchor-motion cohort and production bound qualification
remain open. Later numerical-failure handling changes in the working code do
not retroactively change this launch-time snapshot or its recorded results.

Verification after the replay and typed-failure changes: full fast lane
1536 passed, two optional-dependency skips, 31 slow tests deselected. Lightweight
CI selection with OpenMM/PDBFixer deliberately unavailable: 394 passed and 61
expected skips. Ruff and local documentation checks pass. Numeric recovery is
covered by typed synthetic failure tests, not by failures observed in this cohort.
