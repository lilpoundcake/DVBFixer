# Private gap portfolios and anchor-movement calibration

Date: 2026-11-10. P2 has a private Python implementation and retained-coordinate
controls. Calibration has executed methodological pilots. Neither is a new public
CLI/GUI capability; no higher movement limit has been adopted. P3 is unchanged.

## P2: generated-only multi-parent contract

Continuation: [private publication/exact reuse and factorial symmetry analysis](portfolio-publication-and-anchor-symmetry.md)
now cover private bundle/`.dat`, preparation custody and readback. Public product
integration and scientific qualification remain gated.

`model/diffusion/gap_portfolio.py` defines a separate
`experimental-multiparent-gap-portfolio` kind, not a `RunnerCandidate` or a
schema-4 sampling result. Each parent records the exact common request/source,
candidate label and original seed, coordinate digest, staged sampler-trace digest
and refinement-evidence digest. The research loader verifies retained diagnostic
manifest hashes, trace-bound staged evidence and its generated-only stage before
admission. Moved-anchor donors are not supported in this first version.

All donors must share the exact request, including masks, template ownership,
source bytes and frame. Observations must match the normalized source. Generated
metadata and atom axes must agree. Selection copies only coordinate columns from
one complete donor gap into the immutable base projection; it never fabricates
a single sampling seed for the result.
The frozen original request remains a donor/mask reference; this explicitly
declared post-hoc operation does not rewrite it or claim that it sampled a
portfolio under its historical runner policy.

The deterministic bounded beam evaluates homogeneous parents as controls, then
tries whole-gap substitutions. Detailed local findings prioritize work; every
trial is assessed as a complete complex, so inter-gap clashes and junctions are
not assumed away. Budget exhaustion remains explicit. An issue-local refinement
step independently derives affected gaps, freezes all other generated gaps and
observations, and uses the P0/P1 full-score refiner. Unattributed residual findings
conservatively expand the scope instead of being ignored. Search scores and
runner-provided masks grant no publication authority.

`audit_portfolio_materialization` and a separate-process readback reconstruct the
expected selected coordinates, check lineage, whole-file non-coordinate records,
the independently derived local mobility mask, full geometry and global
L-chirality. Private outputs always retain `publication_qualified: false`.

The final fast suite passed **2,238 tests**, with the three existing optional
Torch/Biotite skips and 33 slow tests deselected. All 15 new focused tests passed;
Ruff passed. These Python-only changes did not rerun GUI or slow integrations.

### Executed controls

No new neural samples or replacement seeds were used. All original cohort
outcomes, including the original 1KTZ 0/5, remain unchanged.

| Case | Search budget | Private outcome |
|---|---|---|
| 1KTZ | 64 whole-complex evaluations, beam width 4 | Multi-parent choice leaves one D136 CA–C–O angle finding (166.2°) plus its aggregate gate; repair of only gap 7 (D131–136) yields a whole-complex geometry/custody pass |
| 1BRS | 16 evaluations maximum | Selects original seed-7 donor for both gaps; no unnecessary recombination or repair |
| 1AHW | 32 evaluations maximum | Selects the already-valid seed-103 donor for all gaps; no invented rescue |
| 1HE8 | 16 evaluations maximum | Still rejected with 37 actionable findings; not an exhaustive search or proof no portfolio exists |

1KTZ's gap donors, in request gap order, are seeds
`103, 104, 7, 103, 101, 7, 103, 103`. The final sixth-residue-gap refinement used
optimizer seed 7, at most 250 iterations per full phase and four seeded restarts.
That optimizer seed is not a new sampling provenance label. All seven other gaps
and every observed record remain unchanged by that operation.

The final PDB SHA-256 is
`b66042f5d2177ed8705a4e32ab2b063c18701ebafd2a3729937c6dc72c3d211b`.
Independent readback found no actionable geometry or admission findings.

Erratum: the first narrative described the pre-local-repair remainder as a clash.
The retained selected report always recorded a **bond-angle** finding instead.
Only this description is corrected; selected/final bytes, lineage, recorded
assessments and original cohort outcomes are unchanged.

Evidence is retained under `tmp/gap-portfolio-p2/`: common requests/source bytes,
five parent coordinate/trace/evidence artifacts, all beam trials and their
assessments, explicit choices, local repair protocol, before/after coordinates,
software snapshots and separate readback reports. This is case-specific
post-hoc workflow evidence, not structural accuracy or broad qualification.

### Reproduction interfaces

These scripts require the separately retained local cohort; those generated
artifacts are not distributed as tracked structural fixtures.

```bash
python scripts/research_gap_portfolio.py --cohort <retained-cohort> \
  --case 1KTZ --output <new-portfolio-directory>
python scripts/repair_gap_portfolio.py --portfolio <portfolio-directory> \
  --output <new-repair-directory>
python scripts/audit_gap_portfolio.py --portfolio <portfolio-directory> \
  --repair <repair-directory> --report <new-readback-report.json>
```

All output directories are new reservations. A malformed donor or changed source
fails closed; failed trials are retained. The public `dvbfixer model`, Homology,
GUI, schema-4 bundles and their readers are unchanged.

### Remaining P2 gates

- Private serialization/publication, generated-only `.dat`, preparation custody
  and exact reuse are implemented in the distinct kind; independent broader
  scientific qualification remains open. No relabelled single-seed bundles.
- A prospective multi-case evaluation with declared donor/repair budgets.
- Explicit public/API/CLI/Homology/GUI opt-in and resource controls, if approved
  for product integration after the scientific/custody gates pass.
- Eligibility for moved-anchor donors would need a new reconciled-frame/mobility
  contract; never mix them under this generated-only version.

## Anchor calibration: measurements, not a new ceiling

`model/diffusion/anchor_calibration.py` uses exact heavy-atom identities and
leave-anchor-out Kabsch fits. Measuring a residue that participates in its own
fit would underestimate its movement. Both a whole-backbone frame and a local
±4-residue frame are reported; degenerate/incomplete fits fail rather than pass.
Case-sensitive chains and insertion codes remain distinct. Contact-pair changes
use a declared 4 Å heavy-atom cutoff, excluding immediate sequence neighbors;
they are descriptive, not energetic or density acceptance.

The candidate grids were declared before coordinate comparisons:

| Grid | Backbone RMSD | Cα maximum | Heavy-atom maximum |
|---|---:|---:|---:|
| Current | 0.5 Å | 0.5 Å | 1.25 Å |
| Exploratory 0.75 | 0.75 Å | 0.75 Å | 1.5 Å |
| Exploratory 1.00 | 1.0 Å | 1.0 Å | 2.0 Å |
| Exploratory 1.50 | 1.5 Å | 1.5 Å | 2.5 Å |

The grid is sensitivity analysis, not a policy choice. No B factor is converted
into a displacement radius. Raw downloads, RCSB metadata, hashes, model choices,
occupancy/altloc/heavy-completeness exclusions, bond/chirality QC and all
per-residue results are retained under `tmp/anchor-calibration/`.

### Observations in the local leave-anchor-out frame

| Comparison | Matched residues | Backbone median / 95th percentile | Current-grid coverage | Exploratory-0.75 coverage |
|---|---:|---:|---:|---:|
| 1UBQ–1UBI, separate X-ray depositions | 72 | 0.090 / 0.149 Å | 61/72 | 61/72 |
| 1D3Z models 1–2, one NMR ensemble | 72 | 0.045 / 0.115 Å | 55/72 | 55/72 |
| 1UBQ–1D3Z model 1, cross-method | 72 | 0.215 / 0.433 Å | 33/72 | 39/72 |
| 1LYZ–193L, initial lysozyme pair | 125 | 0.389 / 0.964 Å | 72/125 | 90/125 |
| 1AKI–193L, QC follow-up / different crystal packing | 125 | 0.162 / 0.586 Å | 99/125 | 100/125 |

These are correlated residue measurements, not independent experiment counts.
NMR and cross-method strata are not pooled with independent crystal comparisons.
Whole-triplet coverage can remain unchanged when backbone limits increase:
sidechain-heavy displacement is often the limiting quantity.

Crucially, initial 1LYZ QC found one D-Cα and three severe bond findings. Its
measurements remain retained and labelled, not silently replaced or used to set
a production quantile. The separately declared 1AKI/193L follow-up had no bond
or chirality findings in these checks. That does not mean all aspects of those
deposited models or their density have been validated.

The evidence supports investigating backbone/Cα and sidechain bounds separately,
with local-context and contact checks. It does **not** support multiplying every
ceiling by a common factor or adopting a high heavy-atom bound: conformational
changes, crystal packing, exact-name symmetry ambiguities and model errors all
contribute. One lysozyme pair does not justify applying its tail to every complex.

### Calibration gates still open

- A larger predeclared experimental panel, including interfaces and more folds.
- Conservative FF-graph symmetry and factorial sensitivity are now implemented
  in the [continuation](portfolio-publication-and-anchor-symmetry.md); broader
  chemistry and prospective calibration remain open. Exact metrics are retained.
- Density/restraint support, covalent links and meaningful contact preservation;
  contact-pair churn alone is neither damage nor validation.
- Independent prospective validation of the separately varied backbone, Cα and
  sidechain grids; retrospective factorial analysis is not that validation.
- Any adopted ceiling needs a new request-bound policy and independent audits.
  Historical limits and frozen request serialization must remain unchanged.

Declarations and sources:
[ubiquitin pilot](anchor-calibration-pilot-v1.json),
[initial lysozyme pilot](anchor-calibration-lysozyme-v1.json),
[QC follow-up](anchor-calibration-lysozyme-qc-followup-v1.json);
RCSB entries [1UBQ](https://www.rcsb.org/structure/1UBQ),
[1UBI](https://www.rcsb.org/structure/1UBI),
[1D3Z](https://www.rcsb.org/structure/1D3Z),
[1LYZ](https://www.rcsb.org/structure/1LYZ),
[1AKI](https://www.rcsb.org/structure/1AKI),
[193L](https://www.rcsb.org/structure/193L).

See the [implementation plan](../plans/diffusion-geometry-repair.md) and
[canonical checklist](../plans/implementation-checklist.md). P3 reconstruction
and public portfolio integration are not implemented by this research route.
