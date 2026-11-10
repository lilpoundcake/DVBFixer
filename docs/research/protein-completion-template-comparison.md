# Matched branched controls: own reconstruction versus seeded PDBFixer

Status: **retrospective, case-specific scientific comparison**, not production
qualification. Follows the [paired kernel](protein-completion-branched.md),
[singleton source controls](protein-completion-anchored-source.md) and
[mixed evidence work](protein-completion-mixed-publication.md).

Continuation: a [48-mask rule-selected, own-success-unfiltered study](protein-completion-unfiltered.md)
now includes source refusals and an own failure. Its native-reference control
exposes a coarse 1–4 clash-gate limitation. The historical results below remain
unchanged; neither experiment establishes calibrated packing or public readiness.

## Design and environment

All **32 previously accepted own-source controls** were declared before comparator
execution: 16 paired masks and 16 observed-partner singleton masks on the existing
explicit 1UBQ A1–72, observed 1FR2 chain-A and 1A3N A1–59 ATOM projections.
This selection is **retrospective and conditioned on own-method success**, not a
blind prospective cohort or a representative sample of missing-atom chemistry.

Both arms receive identical complete masked source/context; withheld coordinates
enter only later reporting, never proposals or selection. Own accepted PDB bytes
were checked against the earlier retained controls. No occupancy rewrite, observed
motion, source repair, hydrogen placement or downstream refinement/reflection.

`scripts/branched_template_baseline.py` uses PDBFixer 1.12.0/OpenMM 8.5.2 on the
explicitly declared **Reference** platform. Preflight actually loads `soft.xml`,
checks required templates and pins all 30 stock template hashes plus the upstream
source and shared helper. The active scientific environment passed; no missing
resources were patched and no failed imports were disguised as science failures.

PDBFixer gap insertion is disabled. The ordinary production seeded helper owns
rebuild calls/seeds 1–5; an instance-local observation hook retains every returned
temporary topology/coordinate state and the helper log. It does not change the
builder, retry budget or coordinates. All 32 inputs introduced exactly one extra
terminal OXT in temporary topology; it is recorded but never inserted into the
diagnostic masked-atom PDB. No other extra heavy atoms were accepted.

This is a **heavy-coordinate reconstruction comparison**, not the full legacy
MODELLER/PDBFixer/PROPKA/Reduce/minimize pipeline. PDBFixer retains its internal
local optimization; there is no claim of equal compute/optimizer budgets. Reference
is a numerical research control, not a new production platform requirement.

## Common gate, distinct from a generator's sampling policy

`scripts/compare_branched_completion.py` applies the same bond/angle/volume/clash,
fixed-source, PDB rounding, global Cα/THR-Cβ and independent RDKit target-CIP checks
to both arms. It uses **measured χ**, not the own method's staggered-grid membership.
The internal research-only `measured_chi` switch does not alter public candidate
audits or source/cohort/evidence policies; an automated non-grid test checks this.

The common geometry gate is centred on the own prototype's declared illustrative
bond/angle values. It is **not an independent calibrated gold standard**, a
MolProbity qualification or a force-field energy criterion. The RDKit CIP checks
provide separate chemical stereo evidence; naming and geometric deadbands remain
explicit. VAL CB naming orientation is not a molecular R/S stereocentre.

All original source records are reinserted unchanged; template temporary observed
motion was at most **1.435×10⁻¹⁴ Å**. Every proposal is retained as diagnostic data,
including failures, rather than silently promoted to accepted output/publication.

## Results: feasibility and accuracy must be read separately

[Full evidence](protein-completion-template-comparison-evidence.json),
[corrected aggregate](protein-completion-template-comparison-summary-02.json),
raw files: `tmp/pdbfixer-research/template-comparison-01/`.

| Mask class | Own common passes | Seeded PDBFixer common passes |
|---|---:|---:|
| Paired branches, with/without CB | 16/16 | 6/16 |
| Singleton, fixed CB/partner | 16/16 | 9/16 |
| Total | **32/32** | **15/32** |

The own arm was selected from its prior accepted controls, so **32/32 is not new
independent success-rate evidence**. Exact-name raw generated-coordinate RMSD in
the unchanged frame, including diagnostic rejected proposals:

| Class | Own median | PDBFixer median | Own lower raw error |
|---|---:|---:|---:|
| Paired | 0.179 Å | 0.992 Å | 10/16 |
| Singleton | 0.085 Å | 0.379 Å | 13/16 |

Do not interpret rejected-proposal RMSD as accepted structure quality. On the
**same subset where both arms pass**, the tradeoff is different:

| Class | Shared subset | Own median | PDBFixer median | Own lower error |
|---|---:|---:|---:|---:|
| Paired | 6 | 0.164 Å | **0.107 Å** | 1/6 |
| Singleton | 9 | **0.069 Å** | 0.133 Å | 6/9 |

Thus PDBFixer is more accurate on **5/6 mutually passing paired cases**. The own
whole-sidechain 1A3N VAL1 remains poor: raw RMSD **1.865 Å** versus template
**1.227 Å**, but the template fails the common geometry/stereo gate. Earlier
1.864 Å refers to rounded own coordinates, not a changed/repaired result.
Scores are exact-name, not symmetry-minimized for VAL; no truth-based reranking.

## Stereo and seeded retry evidence

- Four template proposals fail the common **non-L-or-nearplanar Cα** gate. This
  is **not four confirmed D-Cα**: three have small positive volumes, one negative.
- Independent rounded RDKit labels confirm **one CA=R** and **nine THR CB=S**.
  Ten proposals fail the broader THR-Cβ R-or-nondegenerate gate; the additional
  case has no independent CB CIP assignment, not a confirmed tenth inversion.
- Two proposals fail VAL named-branch orientation. This does not make VAL CB a
  physical stereocentre and is not evidence of two D amino acids.
- The helper accepted its first returned state for 31 requests; 1FR2 THR11 whole
  sidechain retried all five seeds and retained a final D-Cα state with its
  explicit exhaustion warning. All five states/logs are retained. They are
  internal retries, **not five independent scientific samples**.

The helper's deadband/distance gate differs from the stricter research audit.
These findings do not establish a new causal stage or invalidate the production
legacy chirality safety net. Downstream reflect/guarded-minimize safeguards remain
essential and untouched. No production backend switch or dependency removal.

The first aggregate JSON is preserved historically but used ambiguous field names
that conflated nearplanarity with inversion. The corrected `summary-02` separately
reports gate failures and independent CIP counts; raw evidence was not rewritten.

## Reproduction and remaining work

```bash
# Use the intact scientific environment; missing assets are blockers, not skips.
PYTHONPATH=src:. python scripts/compare_branched_completion.py --paired tmp/pdbfixer-research/branched-source-01 --anchored tmp/pdbfixer-research/anchored-source-01 --output tmp/pdbfixer-research/new-template-comparison --summary tmp/pdbfixer-research/new-template-comparison.json
python scripts/summarize_branched_comparison.py --evidence tmp/pdbfixer-research/new-template-comparison.json --output tmp/pdbfixer-research/new-template-aggregate.json
PYTHONPATH=src python -m pytest tests/test_branched_template_comparison.py -q
```

Verification: **551 focused tests passed in 21.69 s**; seventeen-file package mypy
and targeted Ruff passed. Tests cover actual resource usability, fixed-observation
preservation, no staggered-grid penalty, retained seeded states, exact-name scores,
matched accepted-subset statistics and nearplanar-versus-inversion reporting.

Next: declare broader masks/structures **without filtering on own success**, add
independent calibrated geometry/packing metrics and a separately defined full
legacy/refinement comparator if needed. Do not infer broad template superiority,
own-method superiority, root-cause resolution, MD readiness or public replacement
qualification from this deliberately narrow matched-control study.
