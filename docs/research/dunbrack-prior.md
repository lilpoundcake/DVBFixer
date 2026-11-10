# Installed PyRosetta prior-only ranking: access verified, accuracy not rescued

Status: **research-only, retrospective ranking of saved proposals**. User permits
local installed-library use; bundled tables must not enter production or public
artifacts. See the [mandatory production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03).

This study asks whether source-backbone-dependent rotamer-family weights improve
selection among the unchanged [relative-rotation candidates](ncac-relative-rotamers.md).
It includes both that old cohort and the [new-input negative controls](ncac-new-system.md),
whose native outcomes are already known. It is not a newly prospective validation.

## API and active library provenance

PyRosetta revision `1a56185c2592611dec4c9c75ddc9468cd2227c1f` queries
`RotamerLibrary.get_library_by_aa(...).get_all_rotamer_samples(...)` with φ/ψ
computed only from masked-source observations and explicit peptide neighbors.
No pose is imported, no source atom is moved, and no packing, minimization,
hydrogen optimization, score evaluation or new coordinate generation occurs.
Canonical χ1 atom definitions are checked against the installed residue types:
N–CA–CB–CG1 for VAL and N–CA–CB–OG1 for THR. `bin_rotameric_chi` assigns the saved
candidate's measured χ to a family. Each three-candidate pool must cover each
library well exactly once; no probability cutoff or omitted feasible candidate.

The active installed binary is **`rotamer/shapovalov/StpDwn_0-0-0/Dunbrack10.lib.bin`**,
not an assumed copy of the original recommended 5%-stepdown upstream library.
Its SHA-256/size and software revision are bound in the
[selection evidence](dunbrack-prior-evidence-selection.json). No upstream-table
equivalence, standalone-library provenance or new redistribution right is inferred.
Raw per-case API samples, means, sigmas and weights remain private; neither table
bytes nor those query samples are included in the tracked evidence.

## Frozen protocol and explicit raw-weight amendment

`research_dunbrack_prior.py` froze 22 original masks, including parent refusals
and missing-context blockers. Eleven old and four new-input pools have saved
NCAC proposals and source-bound φ/ψ; seven original new-input masks remain
nonexecutions. No substituted inputs or silent context completion is introduced.

The first protocol required API `probability()` values to be normalized. Its first
run stopped before a complete selection ledger: off-grid API values did not sum
to one, and a value elsewhere exceeded one. On these 15 queries the raw sum range
is 0.960705–1.041386. This is an observed API property, **not an established cause
or a claim of a library defect**. Do not label the outputs calibrated probabilities.

The original script/plan/failed run remain unchanged. A separately frozen v2
adapter, `research_dunbrack_prior_weights.py`, uses finite nonnegative **raw family
weights** with positive total. It selects the greatest weight only among common-
geometry-feasible saved candidates, stable original-index tie. The optional
normalized display is diagnostic only and does not enter selection. No Gaussian
offset likelihood, context posterior, ref2015 sum or fitted score weighting is
claimed. No scientific geometry threshold, angle or source coordinate changed.

All 15 complete three-well queries pass φ+360° and ψ−360° controls within 1e−6.
A fresh process repeats the complete executed plan, environment, raw query and
weight-semantics files **byte-identically**. The independent auditor replays
input bindings, geometry feasibility and selection, commits a sanitized
selection ledger, and only then reads native coordinates for accuracy.
Cross-interpreter diagnostic sums tolerate 1e−14 arithmetic differences; raw
API bytes and geometry checks are not loosened.

## Results: family prior alone does not fix poor choices

[Independent accuracy evidence](dunbrack-prior-evidence.json): rounded generated-
atom exact-name RMSD medians, Å.

| Original cohort | Planned / queried | Zero-angle NCAC | Prior-only choice | Parent fixed-heavy scored choice |
|---|---:|---:|---:|---:|
| Old eleven | 11 / 11 | 0.322945 | 0.226284 | 0.226284 |
| New-input eleven | 11 / 4 | 0.438949 | 1.222795 | 1.222795 |

On the four available new-input pools, all prior choices equal the already
unsuccessful fixed-heavy choices, including 1PGB THR11: 0.551633 Å at zero angle
versus 2.119324 Å after selection. On the old cohort, only 1R69 THR26 differs
from the fixed-heavy choice and is worse: **2.083064 versus 0.337043 Å**.
Equal aggregate medians hide that deterioration. All selected candidates remain
geometry-feasible and observation-preserving; that is not native accuracy.

Seven inherited nonexecutions remain in the 22-mask denominator: three 1BPI
altloc/occupancy refusals and four 1CSP missing-context NCAC blockers. No library
query or successful completion is claimed for them. These outcomes support a
working source-bound installed-library access mechanism, not improved accuracy,
broader qualification or a production replacement.

## Retained evidence and continuation boundary

Private paths under `tmp/pdbfixer-research/`: `dunbrack-prior-plan-01.json`, failed
`dunbrack-prior-run-01/` and log; amended `dunbrack-prior-plan-02.json`,
`dunbrack-prior-run-02/` and exact repeat `dunbrack-prior-run-03/`.
The intermediate auditor's strict cross-Python scalar-sum comparison failed
before publication; that diagnostic failure is retained separately.

```bash
PYTHONPATH=src:. python scripts/audit_dunbrack_prior.py \
  --plan tmp/pdbfixer-research/dunbrack-prior-plan-02.json \
  --run tmp/pdbfixer-research/dunbrack-prior-run-02 \
  --repeat tmp/pdbfixer-research/dunbrack-prior-run-03 \
  --output <existing-new-directory>/dunbrack-prior-evidence.json
```

**712 focused tests passed in 26.91 s**; not full repository/platform acceptance.
Production, installed libraries, source fixtures, original scientific evidence
and thresholds remain unchanged. No commit/push or table redistribution.

Potential next experiment must separately freeze actual library-mean proposals
while preserving CB/observations, with unchanged gates and explicit parent
controls. This ranking study generated no such proposals and cannot qualify them.
Calibrated context chemistry, missing-context admission and full legacy comparison
remain open. Original upstream data/license provenance and explicit user approval
must be requested and verified before any production transition.
