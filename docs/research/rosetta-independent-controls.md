# Pre-download new-system controls: ranking gains, retained template stereo failures

Status: **prospectively enrolled follow-up on three deliberately selected new
systems**, not randomized/blinded qualification, a temporal holdout, independence
from Rosetta training/templates, or a replacement for the full legacy pipeline.
Follows the [bounded coupling arm's stopping result](rosetta-reference-controls.md).

**Result:** the unchanged fixed-heavy score improves six correlated masks on
three residues, one in each new structure. All source observations remain fixed.
Two raw template proposals have confirmed D-Cα; a third has small signed volume
and no CA CIP assignment. All three are rejected, not repaired or published.
Matched template accuracy does **not** establish own-method superiority.

## Enrollment, inputs and complete denominators

Before coordinate download/inspection or new proposals, the enrollment specified:

| Source | System | Selected chain | Actual mask cases |
|---|---|---|---:|
| 2CI2 | Barley chymotrypsin inhibitor 2 | I | 16 |
| 1R69 | Phage 434 repressor amino-terminal domain | A | 16 |
| 1CTF | E. coli L7/L12 carboxy-terminal domain | A | 12 |

Rule: first coordinate MODEL, first encountered ATOM chain (case-sensitive), up
to the first two VAL and two THR identities in file order, four masks each:
distal pair, whole sidechain including CB, first distal atom, second distal atom.
**1CTF has only one observed THR**, so there are 44 actual cases, not 48. Its
shortfall is retained; no replacement structure, chain, residue or mask was added.

Raw source PDBs are unchanged RCSB downloads from 2026-11-03 under
`tests/fixtures/protein_context/independent_followup/`, with
[provenance/checksums](../../tests/fixtures/README.md). Projection is explicit,
protein-chain-only, not intact-complex chemistry acceptance. It preserves ATOM
bytes including H/altloc/occupancy, chain TER and relevant original linkage
records. CONECT follows coordinates; excluded link endpoints remain for
fail-closed admission. No source repair, H stripping or occupancy rewrite.

All 44 cases pass shared admission and execute both arms. No baseline execution
or protocol blocker occurred. The PDBFixer resources/Reference preflight succeeded
in the existing environment (PDBFixer 1.12.0, OpenMM 8.5.2); no environment repair,
installation or production platform change was made.

## Unchanged proposals, scoring and gates

- Own arm: original three-staggered VAL/THR proposals or the exact fixed-partner
  anchored closure. Clearance selection is the old baseline; no nine-χ expansion.
- Template arm: unchanged shared seeded PDBFixer helper on Reference, seeds 1–5
  as needed, with every returned temporary state retained. No chirality reflection,
  relaxation, Modeller gap reconstruction or downstream production safety net.
- Same existing measured-χ common gate for both arms, including source records,
  raw/rounded geometry, global Cα/THR-Cβ and independent target RDKit CIP.
- Same full-pose `ref2015`, frozen Rosetta build/seed/import flags, fresh
  unoptimized H, typed virtual separation, bounded terminal OXT and complete
  geometry-feasible own-pool selection. Native controls run after selection is
  committed; truth never generates or chooses candidates.

The first score invocation used the wrong CLI subcommand `score`; argparse exited
before any science or run-directory creation. Its log is retained. The corrected
invocation uses `run` with unchanged code/policy. This is an invocation blocker,
not a scientific refusal, protocol amendment or failed sampling seed.

[Independent cohort evidence](rosetta-independent-controls-evidence.json) replays
saved proposals/requests, audits, clearance choices, template seed snapshots and
exact observation preservation without rebuilding coordinates. The unchanged
[score readback](rosetta-independent-controls-readback.json) independently checks
all **176 scored poses/snapshots**, physical-heavy identity/coordinates, canonical
target/external graphs, fixed context and score terms. No context-coordinate
comparison blocker occurs. Its historical benchmark flags remain unchanged;
prospective enrollment provenance is recorded separately, not retrofitted into
the generic reader's old scope labels.

## Coverage and native-control limitations

| System | Own feasible selections | Template common-gate passes | Native common-gate passes |
|---|---:|---:|---:|
| 2CI2 | 15/16 | 4/16 | 13/16 |
| 1R69 | 16/16 | 8/16 | 16/16 |
| 1CTF | 10/12 | 4/12 | 8/12 |
| Total | **41/44** | **16/44** | **37/44** |

Score changes ranking, not feasible coverage. Three own pools remain empty:
2CI2 VAL32 CG1 and 1CTF THR76 OG1/CG2 singleton masks. Their native controls do
not rescue them. Diagnostic rejected proposals remain retained and are never
made public bundles or automatically loaded GUI structures.

**This illustrative gate is not a calibrated native gold standard.** Seven native
failures are correlated masks of two source residues:

- 2CI2 VAL32: deposited CG1–C distance **2.758456 Å**, below the unchanged
  C–C cutoff of 2.8 Å; three masks fail. It is a known 1–4 contact, not a reason
  to introduce a blanket 1–4 exclusion or relax a frozen cutoff.
- 1CTF THR76: deposited CA–CB–CG2 angle **125.118487°**, outside the existing
  illustrative tetrahedral target ±12°; all four masks fail.

Do not label these deposited deviations reconstruction damage or proven chemical
impossibility. Idealized own construction can preferentially satisfy its geometry
gate; **41 versus 16 is filter-specific coverage, not proof of superior accuracy**.
The sixteen mutually passing masks have rounded exact-name RMSD medians
**0.198 Å own/scored versus 0.190 Å template**. Neither superiority nor full
legacy equivalence follows from this small, selected, correlated denominator.

## Prospective ranking result on the new structures

Six score choices change, all improving held-out rounded exact-name RMSD:

| Residue | Distal-pair RMSD, clearance → score | Whole-sidechain RMSD |
|---|---:|---:|
| 2CI2 VAL28 | 2.435 → 0.124 Å | 2.093 → 0.168 Å |
| 1R69 THR26 | 2.536 → 0.269 Å | 2.056 → 0.383 Å |
| 1CTF VAL66 | 2.609 → 0.277 Å | 2.182 → 0.352 Å |

These are **three residue-level events, not six independent successes**. Across
all 22 matched paired masks, median RMSD improves **0.524 → 0.249 Å**; nineteen
matched anchored masks remain **0.128 Å** because they have one algebraic option.
Generated-only VAL methyl symmetry is additionally reported, never used for
selection; fixed partners and THR branches cannot swap. The VAL ranking gains
persist under allowed generated-only symmetry. Three small classic systems do
not establish transfer to broader chemistry, repeat contexts or production inputs.

## Raw template chirality: two confirmed inversions, one unassigned collapse

The helper retries and then explicitly warns that it keeps the last attempt for
downstream safeguards. The comparator intentionally has no such repair. All five
states per failing case are retained, among **56 total template seed states**:

| Whole-sidechain case | Raw signed Cα volume (Å³) | Rounded target CIP | Interpretation |
|---|---:|---|---|
| 2CI2 VAL28 | −0.627386 | CA R | Confirmed D-Cα |
| 2CI2 THR55 | −0.555162 | CA R, CB S | Confirmed D-Cα and wrong THR-Cβ assignment |
| 1R69 THR27 | −0.082474 | CA unassigned, CB S | Small-volume collapse; **not** confirmed D-Cα |

All three fail the independent common gate. No reflection, observation restoration
or relabeling hides them. All 41 selected own candidates satisfy global/target
stereo checks. The observed backbone/other atoms remain unchanged in both arms.
This is evidence about **raw template reconstruction plus retry exhaustion**,
not unsafe final output from production legacy prepare/minimize, whose separate
chirality safeguards remain intact. No optimizer-stage root cause is established
by these retained endpoints alone.

## Reproduction, verification and next gate

Retained project-local roots: `independent-enrollment-01.json`,
`independent-download-01.json`, `independent-masks-01/`,
`independent-comparison-01/`, `independent-score-plan-01/`,
`independent-score-run-01/` below `tmp/pdbfixer-research/`. Only the three original
PDB inputs are tracked as fixtures; masked/generated candidates/logs stay private.

```bash
PYTHONPATH=src:. python scripts/prospective_independent_cohort.py enroll --output <enrollment.json>
# Download exactly the enrolled raw PDB URLs; preserve original bytes/hashes.
PYTHONPATH=src:. python scripts/prospective_independent_cohort.py freeze --enrollment <enrollment.json> --inputs <originals> --output <masks>
PYTHONPATH=src:. python scripts/prospective_independent_cohort.py run --plan <masks> --output <comparison>
PYTHONPATH=src:. python scripts/prospective_independent_cohort.py prepare-scores --enrollment <enrollment.json> --comparison <comparison>/evidence.json --output <score-plan>
# Existing user-authorized licensed PyRosetta interpreter:
PYTHONPATH=. python scripts/research_rosetta_fixed_scoring.py run --plan <score-plan> --output <scores>
PYTHONPATH=src:. python scripts/summarize_rosetta_fixed_scoring.py --plan <score-plan> --run <scores> --output <readback.json>
PYTHONPATH=src:. python scripts/audit_independent_cohort.py --masks <masks> --comparison <comparison> --score-plan <score-plan> --scores <scores> --score-readback <readback.json> --output <audit.json>
```

**658 focused tests passed in 26.67 s**; targeted Ruff passed. Tests cover
pre-download enrollment without science preview, first-model/first-chain identity,
case sensitivity, retained H/altloc/occupancy/link endpoints, template observation
movement/addition tampering, and the tracked 1CTF count shortfall. No full repository
suite or public capability qualification is claimed.

Next useful causal controls are the two new confirmed-D whole-sidechain cases:
separately specify an initial-overlay versus optimizer-stage probe with unchanged
observations and no installed-library/production edits. Keep the third unassigned
case distinct. More empirical candidate priors, calibrated context geometry and
full-legacy/independent-family qualification remain open; do not enable production
partial completion from this bounded ranking result.

Final cohort and scoring readbacks reproduce byte-identically without rebuilding
or rescoring. Package mypy, generated CLI/GUI checks, agent maps and all earlier
evidence verifiers pass. All 52 fixture hashes verify. Final read-only checks:
`tmp/pdbfixer-research/validate-independent-final.py`.
