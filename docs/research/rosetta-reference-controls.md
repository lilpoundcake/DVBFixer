# All-reference replay: feasible nonadditivity without a demonstrated selection gain

Status: **retrospective descriptive analysis of saved matrices**, not new sampling,
new energy evaluation or independent prospective qualification. Continues the
[nine-χ controls](rosetta-dense-controls.md); original evidence and selected outputs
remain unchanged.

**Result:** nonadditivity survives on 255 completely geometry-feasible rectangles.
But every conditional factorized comparison using only feasible energy terms
chooses the same candidate as the full score. Earlier differing choices depend on
geometry-rejected reference-axis terms. This separates measurable interaction
effects from an unestablished practical advantage of coupled selection.

## Fixed descriptive protocol

`scripts/analyze_rosetta_references.py` records a separate analysis plan before
executing this new calculation. Existing matrices and their earlier outcomes
were already seen: the plan is explicitly **retrospective**, not a prospective
scientific enrollment. It binds the dense audit, original raw scores, coarse
scores, input plan and implementation hashes.

For every one of the twelve frozen matrices:

1. Retain every source-feasible cell, including joint rejections, and every
   possible reference pair `(a,b)`. Require a complete finite Cartesian table.
2. Recompute the conditional surrogate
   `S_ab(i,j) = E(i,b) + E(a,j) - E(a,b)` over the same global joint-feasible set.
3. Classify references by **all terms actually used**, not anchor feasibility
   alone: rejected anchor; feasible anchor with rejected terms; all terms feasible.
4. For each feasible anchor, also restrict both full and surrogate selection to
   the same supported subset where all four corners `(a,b)`, `(i,b)`, `(a,j)` and
   `(i,j)` are geometry-feasible. Record whether this subset retains the global
   full-score winner; never silently compare different coverage.
5. Enumerate all distinct 2×2 rectangles with four feasible corners and report
   `E(i,j)+E(k,l)-E(i,l)-E(k,j)`, independently of any distinguished reference.

Only saved indices, eligibility and total scores enter selection. Native
coordinates/RMSD are not consumed. Stable index breaks exact energy ties.
The fixed **1e-9 Rosetta-unit tolerance is reporting-only**, not a new geometry
cutoff, selector adjustment, significance test or accuracy qualification.
Single-axis matrices have no nondegenerate rectangle; their maximum is `null`,
not zero interaction. No favorable reference or rectangle is chosen after replay.

## Denominators and all-reference results

[Complete analysis](rosetta-reference-controls-evidence.json) and
[independent raw-cell audit](rosetta-reference-controls-readback.json) retain all
references and feasible rectangles, not just extrema or successful comparisons.

| Measure | Count |
|---|---:|
| Frozen cases / two-variable matrices | 12 / 7 |
| Source-feasible cells = reference choices | 141 |
| Joint-feasible cells = feasible anchors | 129 |
| Distinct possible rectangles | 381 |
| Fully feasible / rejected-corner rectangles | 255 / 126 |
| Fully feasible rectangles above reporting tolerance | 255 |

These are correlated masks/references/rectangles in one previously studied
lysozyme system, **not** independent biological replicates or 255 validation cases.

Global feasible-pool surrogate comparisons:

| Reference class | References | Different choice from full score |
|---|---:|---:|
| Anchor itself rejected | 12 | 8 |
| Anchor feasible, but some required axis terms rejected | 25 | 2 |
| Every used energy term feasible | 104 | **0** |

The last row includes eighteen single-axis references and **86 references in
genuinely two-variable matrices**. Thus the zero disagreement is not solely a
degenerate-axis effect. All ten differing choices use rejected energy terms;
anchor feasibility alone is insufficient protection.

Supported-pool comparisons cover all **129 feasible anchors**, with **zero**
full-versus-surrogate choice differences and zero positive score gaps. The global
winner remains in 121 supported pools; eight exclude it. Those eight compare
their own restricted-set full optimum, **not** the original global optimum and
not a claim of unchanged global coverage. Rejected anchors have no supported pool.

## Nonadditivity on fully feasible rectangles

| Pair/mask | Fully feasible rectangles | Maximum absolute rectangle Δ |
|---|---:|---:|
| THR89/VAL92 distal | 30 | 0.427269 |
| VAL2/THR40 distal | 108 | 0.006413 |
| VAL2/THR40 whole | 30 | 0.008327 |
| THR43/THR51 distal | 53 | 3.008091 |
| THR43/THR51 whole | 10 | 0.401633 |
| THR51/THR69 distal | 18 | 0.065705 |
| THR51/THR69 whole | 6 | 0.035662 |

All seven two-variable cases contain at least one fully feasible rectangle. The
remaining five have one-option axes, not a failed or negative interaction test.
The old first-reference maxima **48.693 / 5.441** for THR43/THR51 included rejected
corners. Their feasible-only maxima are **3.008 / 0.402**, not zero. This is a
descriptive change of comparison domain, not a retuned gate or replacement of the
historical numbers. Full-score Rosetta units are neither kcal/mol nor isolated
atom-pair interactions; target-H contributions remain part of the frozen protocol.

**Do not infer:** absent interactions, general equivalence of independent/coupled
packers, template/scorer superiority or public reconstruction qualification.
The specific conditional surrogates tested here select the same minima when their
terms are feasible. Different conditional surrogates, richer chemistry, other
structures and actual packer algorithms have not been qualified.

The earlier **zero changed coarse-to-dense winners** remains unchanged. No new
structure, `.dat`, selected candidate or GUI/public bundle is produced here.

## Independent verification and reproduction

Both new scripts use only the Python standard library. The independent
`scripts/audit_rosetta_references.py` reads the **original raw score cells**, not
the analyzer's audit-derived matrix, and enumerates **sets of four feasible cells**,
not axis-pair products. It checks every reference, supported pool, selector,
rectangle, class and aggregate against the analysis. No coordinate generator,
geometry relaxer, PyRosetta runtime or native-accuracy selector is imported.

```bash
python scripts/analyze_rosetta_references.py plan --audit docs/research/rosetta-dense-controls-evidence.json --output <new-plan.json>
python scripts/analyze_rosetta_references.py analyze --plan <new-plan.json> --output <new-analysis.json>
python scripts/audit_rosetta_references.py --plan <new-plan.json> --analysis <new-analysis.json> --output <new-audit.json>
```

All destinations are exclusive. Retained local files:
`tmp/pdbfixer-research/reference-plan-01.json`, `reference-analysis-01.json` and
`reference-audit-01.json`; published copies are byte-identical to those results.
**649 focused tests passed in 26.39 s** and targeted Ruff passed. Tests cover
additive and coupled controls, rejected reference terms, excluded global winners,
independent rectangle enumeration, omitted/tampered records, numeric reporting,
truth blindness, input immutability and unchanged single-axis interpretation.
This is not a new full repository or scientific-stack acceptance run.

## Decision and next gate

This bounded coupling arm supports a **stopping result**, not increasingly dense
sampling of the same four pairs to obtain a positive packing claim. Keep joint
geometry checks, but do not introduce an approximate graph or production route
based on these score tables. Any further accuracy qualification needs a separately
specified independent-family/reference-free comparison with all refusals retained;
context-chemistry calibration and empirical candidate priors remain separate open
questions. Do not relabel this saved-matrix analysis as new prospective evidence.

Final analysis and independent audit copies reproduce byte-identically. Package
mypy, generated CLI/GUI checks, agent maps, old evidence and all 49 fixture hashes
also pass. Final read-only checks are retained in
`tmp/pdbfixer-research/validate-reference-final.py`.
