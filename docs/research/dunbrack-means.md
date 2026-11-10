# Installed-library mean proposals: geometry retained, new-input ranking worse

Status: **research-only retrospective experiment**, not a production backend.
Continues [prior-only ranking](dunbrack-prior.md) with separately frozen coordinate
generation. The installed PyRosetta library remains a local research dependency;
the [original-table/license/user-approval gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
still applies before any production transition. No table rows or per-case library
weights/means/sigmas are included in tracked publications.

## Frozen intervention

For each of the same fifteen executable NCAC pools, use the three per-case χ1
means already queried from the installed StpDwn_0-0-0 library. Generate one proposal
per library well, in well order. Rotate both generated distal branches together
around the saved CA–CB axis to the desired measured χ. **CB, backbone and every
observation stay fixed.** No packing, H optimization, minimization, sigma sampling,
new library lookup, native-derived angle fit or geometry-threshold adjustment.

`scripts/research_dunbrack_means.py plan` froze source/query/implementation hashes
before generation. All 22 original masks remain in the denominator: eleven old
and four new-input pools execute; three 1BPI refusals and four 1CSP missing-context
blockers remain inherited nonexecutions. Both cohorts have already been evaluated;
this is not a new prospective accuracy test or broad qualification.

Every proposal receives the unchanged raw/rounded/global/target common gate.
Compare two separate selectors over the complete feasible mean pool: greatest
raw library-family weight, or smallest unchanged `ref2015` fixed-heavy score.
Neither reads native coordinates; stable original well index breaks ties. A
missing feasible score blocks the score choice. No prior/ref2015 weighted sum.

## Generation and independent evidence

- **45 proposals**, **42 common passes**, **three geometry refusals**; all retained.
- Fifteen score views, 42 independently audited fixed-heavy snapshots/imports.
- All selected outputs have L-Cα and correct THR CB stereochemistry under the
  common gate. Observed source records and raw CB coordinates remain unchanged.
- Rejected means: 2CI2 VAL32 wells 1/3 clash with I/33 OE2; THR55 well 2 clashes
  with I/69 CE2. They are never scored. Library well indices are not the parent
  relative-rotation indices; do not compare them as identical enumerations.

The [sanitized selection ledger](dunbrack-means-evidence-selection.json) is committed
before native reporting; [accuracy evidence](dunbrack-means-evidence.json) binds it
by digest. `scripts/audit_dunbrack_means.py` generates no coordinates and evaluates
no scores: it checks saved axial/radial transport, pair rigidity, declared/actual
mean χ, source/CB ownership, every common audit and full feasible score inventory,
physical/target/external graph import, score terms and selected-index replay.
Scorer extras/virtuals/context H are identical within every mean proposal pool.
Cross-arm equality of all scorer-owned context coordinates versus the older
relative proposal pools is not separately established here.

Subsequent [saved-context and residue-pair diagnosis](mean-context-pairs.md)
closes that cross-arm context check on all 84 retained feasible snapshots and
checks VAL symmetry. It does not alter this study's frozen results or selectors.

## Accuracy

Matched rounded generated-atom exact-name RMSD medians, Å:

| Original cohort | Planned / executed | Mean prior-only | Mean fixed-heavy scored | Parent relative scored |
|---|---:|---:|---:|---:|
| Old eleven | 11 / 11 | 0.198491 | 0.198491 | 0.226284 |
| New-input eleven | 11 / 4 | 1.210568 | **1.930166** | 1.222795 |

Equal old-cohort medians hide different choices: on 1R69 THR26 the mean prior
chooses a 2.022584 Å result, while score selection chooses 0.240256 Å. The new-
input cohort exposes the opposite scoring failure: **1PGB VAL29** has a good mean
proposal selected by the prior at **0.078348 Å**, but score chooses another at
**1.770814 Å**. This worsens the four-case median compared with the parent score
experiment. 1PGB THR11 remains poor at 2.089519 Å under both mean selectors.

The old-cohort improvement is bounded and retrospective; it does not outweigh
the negative transfer or establish score calibration. Family-prior selection is
also not a general answer, as THR26 illustrates. Native geometry caveats and all
seven inherited failures from the parent studies remain unchanged and explicit.
No native gate, refusal denominator or counterexample is removed.

## Finite saved-score diagnosis, no weight refit

`scripts/analyze_dunbrack_mean_terms.py` selects **all two prior-versus-score
disagreements**, without reading native data, and decomposes the recorded weighted
score difference (score-selected minus prior-selected). No new scores or weights:

- **1PGB VAL29:** total −0.404808 Rosetta score units. Attraction (`fa_atr`,
  −1.597526) and solvation (`fa_sol`, −0.867703) outweigh the unfavorable library
  term (`fa_dun`, +1.012951), repulsion (+0.790042) and other small changes. The
  false native preference is a balance of these recorded terms, not absence of
  an empirical prior from `ref2015`.
- **1R69 THR26:** total −9.906042, dominated by lower `fa_rep` (−9.448179), despite
  unfavorable `fa_dun` (+1.023393). Here context scoring overrides the family prior
  toward the more accurate candidate.

Both term-delta sums reproduce saved total differences. These are score-model
arithmetic explanations, not proof of physical/native causality or a prescription
to change weights. In particular, simply increasing the prior term to rescue
VAL29 could oppose the useful context correction in THR26. Private readback:
`tmp/pdbfixer-research/dunbrack-means-term-deltas-01.json`.

## Retained artifacts and checks

Private paths under `tmp/pdbfixer-research/`: `dunbrack-means-plan-01.json`,
`dunbrack-means-generation-01/` (private geometry and runner-facing score plan),
`dunbrack-means-score-run-01/`. Per-case library-derived metadata stays private;
publications contain selection indices, geometry findings and scientific outcome
metrics only, not raw library rows/parameters.

```bash
PYTHONPATH=src:. python scripts/audit_dunbrack_means.py \
  --plan tmp/pdbfixer-research/dunbrack-means-plan-01.json \
  --generation tmp/pdbfixer-research/dunbrack-means-generation-01 \
  --score-run tmp/pdbfixer-research/dunbrack-means-score-run-01 \
  --output <existing-new-directory>/dunbrack-means-evidence.json
```

Use fresh destination files. **720 focused tests pass in 26.88 s**, including the
saved-score arithmetic controls; this is not the full
repository/platform suite. No source fixture, installed library, production
consumer/default/dependency or public capability changed. No commit/push.

Further work must not fit weights to these native errors. Separate context
completion, a genuinely new
evaluation cohort, calibrated chemistry and a full legacy comparator remain open.
Do not add denser sigma/angle sampling to manufacture a positive cohort claim.
