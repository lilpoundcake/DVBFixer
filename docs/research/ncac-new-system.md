# New-input NCAC controls: accuracy gain does not transfer

Status: **prospectively enrolled inputs, negative bounded result**. Continues the
[retrospective relative-rotamer experiment](ncac-relative-rotamers.md), without
altering its algorithm, three relative angles, common gate or fixed-heavy scorer.
No production capability or general PDBFixer replacement is qualified.

## Enrollment and complete denominator

Before downloading coordinates, `scripts/prospective_ncac_cohort.py enroll` froze
1PGB (protein G B1), 1BPI (bovine pancreatic trypsin inhibitor) and 1CSP (cold-shock
protein). First coordinate model, first encountered ATOM chain and first two VAL/
first two THR identities in file order are selected. Only whole-sidechain masks
are considered; no replacement system or outcome-dependent filtering is allowed.
Original H, altloc/occupancy and touching linkage metadata are preserved by the
unchanged projection. This is deliberate enrollment, not randomized/blinded or
independent of Rosetta training/templates.

Original downloads and SHA-256 provenance are retained under
[`tests/fixtures/protein_context/ncac_followup/`](../../tests/fixtures/protein_context/ncac_followup/)
and described in the [fixture inventory](../../tests/fixtures/README.md).
The 12-mask maximum produces **11 actual masks**: 1BPI has only one selected VAL.
That shortfall remains explicit.

| System | Enrolled masks | Shared refusals | NCAC evaluated | NCAC protocol blockers | Selected NCAC / baseline |
|---|---:|---:|---:|---:|---:|
| 1PGB | 4 | 0 | 4 | 0 | 4 / 4 |
| 1BPI | 3 | 3 | 0 | 0 | 0 / 0 |
| 1CSP | 4 | 0 | 0 | 4 | 0 / 0 |

All three 1BPI requests refuse altloc/partial occupancy under the unchanged shared
contract. No alternate state is silently selected. All four 1CSP inputs have
missing **context** sidechain atoms beyond the target mask: stock refuses its
unrequested additions, and the strict NCAC adapter refuses an unexpected overlay
anchor set. These are retained execution/protocol blockers, not geometry passes,
not silently repaired sources and not successful NCAC trials.

## Independent readback and scorer refusal

The [input-only selection ledger](ncac-new-system-evidence-selection.json) is
published before native reference PDBs are opened for reporting. The
[accuracy evidence](ncac-new-system-evidence.json) binds it by digest. Generation
and scoring never use native target coordinates for candidate/seed/score choice.
Scorer plans transport separate complete NCAC/baseline feasible pools and no
native control cells. All case failures and partial-run files receive retained
digests; incomplete failed traces are not presented as complete audited runs.

Four complete NCAC traces contain 20 stage snapshots, four seed returns and eight
used-anchor fits. Independent stage/observation, source/template fit transport,
relative-rotation and common-gate audits pass for these executed cases. All twelve
1PGB NCAC relative proposals pass geometry. Stock common passes only two of four
1PGB outputs. This does not turn the four blocked 1CSP masks into executed successes.

Thirty score imports are independently inspected: **24 scored, six import-contract
blockers**. The six 1CSP baseline imports fabricate sixteen context GLU heavy atoms
(CG/CD/OE1/OE2 at A/3,21,36,66) and move observed A/PRO58 CD. They are correctly
ineligible, with no usable total score. Observations are not restored afterward to
hide the movement. Complete feasible-pool selection therefore blocks all four
1CSP baseline views. No scorer-owned missing-context exemption was added.

For the four selectable 1PGB masks, scorer-owned extras/virtual/context H are
identical across each pool; selected candidates retain source records and L-Cα.
That confirms the import/ownership protocol, not the accuracy of the score ranking.

## Accuracy: failure to reproduce the retrospective improvement

Rounded generated-atom exact-name RMSD, Å, on the four mutual selected cases:

| Mask | NCAC zero angle | NCAC scored | Prior-method scored |
|---|---:|---:|---:|
| 1PGB THR2 | 0.326266 | 0.326266 | 0.230854 |
| 1PGB THR11 | 0.551633 | 2.119324 | 2.107724 |
| 1PGB VAL21 | 2.360604 | 2.360604 | 2.214812 |
| 1PGB VAL29 | 0.076822 | 0.076822 | 0.236430 |
| **Median** | **0.438949** | **1.222795** | **1.172077** |

The only changed NCAC choice is THR11 +120°, which **worsens** native RMSD from
0.551633 to 2.119324 Å. It is chosen without native information and fails to
transfer the earlier four-improvement result. Both scored arms perform poorly
on THR11 and VAL21. No claim of NCAC superiority or calibrated ref2015 selection
is supported. These four masks are one effective system, not three independent
successful systems; denominators remain 11 enrolled / 8 admitted / 4 selectable.

Native common geometry passes 6/8 admitted masks: 1PGB THR11 has a generated
CG2–backbone-C clash and VAL29 fails CB–CG2 length. All four 1CSP native targets
pass their local/common audit despite the incomplete context blocking import.
The native geometry caveats and unsupported-context failures are distinct from
native RMSD; none are removed or used to loosen thresholds.

## Reproduction, checks and next boundary

Private retained paths under `tmp/pdbfixer-research/`:
`ncac-new-system-enrollment-01.json`, `ncac-new-system-download-01/`,
`ncac-new-system-masks-01/`, `ncac-new-system-run-01/`,
`ncac-new-system-score-run-01/`. Only original source PDBs are fixture assets.
All previous evidence and original source bytes remain intact.

```bash
PYTHONPATH=src:. python scripts/audit_prospective_ncac.py \
  --masks tmp/pdbfixer-research/ncac-new-system-masks-01 \
  --science-run tmp/pdbfixer-research/ncac-new-system-run-01 \
  --score-run tmp/pdbfixer-research/ncac-new-system-score-run-01 \
  --output <existing-new-directory>/ncac-new-system-evidence.json
```

Use fresh destinations; the selection-ledger filename is derived automatically.
**696 focused tests pass in 27.66 s**. This is not the full repository suite or
cross-platform acceptance. No installed-library, production default, dependency,
public bundle or GUI change; no commit/push.

Final independent selection/accuracy readback repeats byte-identically. All twelve
historical/current final verifiers, all 55 source-fixture hashes and 1,126 local
Markdown targets pass; targeted Ruff, package mypy, generated CLI/GUI checks,
agent maps and whitespace pass. Final verifier:
`tmp/pdbfixer-research/validate-ncac-new-system-final.py`.

Do not densify these same angles or fit weights to THR11. Further generalization
requires separately declared context completion and better calibrated empirical
priors/selection, retaining original sources and every new failure. Authentic
Dunbrack asset/license lineage remains unavailable locally: repeated official-site
fetches and a direct HTTPS connection timed out. Search snippets alone do not
establish a downloadable licensed asset; bundled PyRosetta tables were not extracted
in the experiments reported here.

### Subsequent research authorization (2026-11-03)

The user has explicitly authorized subsequent local research using the library
already installed with PyRosetta, **not production**. Original upstream archive
availability is therefore not a blocker for that research route. Keep installed
library/version provenance explicit and do not redistribute bundled tables or
add a production dependency. This authorization does not retroactively change
the experiment above or establish broader license rights.

The user will try to obtain original licensed tables. Before any production
integration, packaging or public capability, the agent **must ask the user** for
those original tables and documented source/version/license permissions, verify
their suitability for the intended use/distribution, and obtain explicit
production approval. Research permission is not production permission; scientific
qualification remains a separate requirement. Canonical policy:
[research authorization and mandatory production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03).

**Official-source retry, 2026-11-03 16:54 UTC:** the host is now reachable, but no
table archive has been obtained. The [official overview](https://dunbrack.fccc.edu/lab/bbdep2010)
states CC BY 4.0; its [linked license/application](https://dunbrack.fccc.edu/lab/bbdep2010_lic)
contains ODC Attribution terms and says an approved applicant receives a download
link by email. Example download links lead to that form. No application was
submitted on the user's behalf. Original HTML/provenance is retained under
`tmp/pdbfixer-research/dunbrack-official-retry-03/`. This does not change the
research authorization, production gate or reported scientific results.
