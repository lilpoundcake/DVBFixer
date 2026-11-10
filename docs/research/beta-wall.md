# Private THR Cβ wall: correct handedness alone does not ensure acceptable geometry

Status: **bounded mechanism/stress research**, not a production fix or a causal
replay of the original unforced THR2 failure. Continues the
[complete-start controls](legacy-start-controls.md).

## Frozen experiment, not a new seed search

`research_beta_wall.py` binds all **12 existing complete starts**: six source
masks, both previously declared seed labels. No new H placement, seeds, force
strength grid, native reference coordinates or accuracy-based choices.

Four arms were fixed before execution:

1. **unchanged** — original serialized physical system;
2. **wall** — original system plus a target-THR signed-volume wall;
3. **drive** — original system plus an explicitly artificial inversion drive;
4. **drive + wall** — identical artificial drive, with the wall added.

All use **Reference as a private numerical control**, identical saved starting
coordinates and two native minimization phases, tolerance 10 kJ mol⁻¹ nm⁻¹ and
maximum 1,000 iterations per phase. Phase 2 reduces only the original positional
restraint strengths tenfold, without changing their anchors, just as production
does. No emergency reflection, observation restoration or branch swapping occurs.
This is not the complete production CLI or a replacement production platform.

The wall and drive use the dimensionless target volume
`q = ((CA−CB) × (OG1−CB)) · (CG2−CB) / (2.5 Å³)`:

- wall: `0.5 × 10000 × max(0, 0.4−q)²` kJ/mol;
- artificial drive: `0.5 × 1000 × (q+1)²` kJ/mol.

Thus the wall's requested floor is **1.0 Å³**. Its strength uses the existing Cα
guard's scale as one predeclared research setting; neither parameter is fitted
to native accuracy or subsequently tuned. It is a **finite soft penalty**, not a
hard constraint or mathematical guarantee of that floor or chirality.
Reference-only primitive tests check units/energies and **36 Cartesian derivative
comparisons** across active/inactive wall and drive conditions.

## Admission and original denominators

Initial target THR must have positive β volume and CB R. A positive but small
volume is deliberately allowed for this repair/stress experiment; this is **not
new production input admission**. An already wrong-handed source is refused,
not preserved by assigning it a negative reference volume or silently reflected.
VAL has no stereogenic Cβ and is explicitly not applicable.

Of 12 retained starts per arm:

- **six execute**: THR27, THR1 and THR2 under both old labels;
- **two are refused**: THR6 already has CB S;
- **four are not applicable**: the two VAL masks under both labels.

The two THR27 labels share an identical complete start. Six executions per arm
are therefore **not six independent structures**; the executed targets come from
three masks in two proteins, with five distinct coordinate starts.

Across all arms: **48 declared cells, 24 executed cells / 48 native minimizations,
eight initial refusals and 16 not-applicable cells**, zero execution blockers.
Every executed arm starts with exactly the saved heavy and H coordinates.

## Results

[Independent evidence](beta-wall-evidence.json), with all initial refusals and
not-applicable cells retained in the denominators:

| Arm | Executed | Final target CB R | Final raw-and-rounded geometry passes |
|---|---:|---:|---:|
| Unchanged | 6 | 6 | 6 |
| Wall | 6 | 6 | 6 |
| Artificial drive | 6 | 0 | 0 |
| Artificial drive + wall | 6 | 6 | **2** |

Without a β wall, the artificial drive produces CB S in all six executed cells,
while target Cα remains S. Adding the wall keeps target CB R in all six tested
stress endpoints. This directly demonstrates a **separate Cβ failure mode not
excluded by Cα protection**, under an explicitly imposed drive. It does **not**
reconstruct why the original ordinary THR2 run crossed.

The stress-plus-wall β volumes are approximately **0.707–0.712 Å³**, below the
requested 1.0 Å³ soft floor. The wall therefore does not enforce its floor exactly.
Correct handedness is also insufficient for the other geometry checks:

- **THR2, both old seed labels:** raw and rounded geometry pass;
- **THR27 and THR1, both labels:** fail the existing **CA–CB–CG2 bond-angle gate**,
  despite CB R and positive volumes above the near-planarity threshold.

Those four endpoints are explicitly **refused for geometry**, not called repaired
or accepted merely because CB R was preserved. No strength/budget/angle threshold
is adjusted after seeing that result. The unchanged and non-stress wall arms both
pass all six geometry checks, but that is not evidence of improved native accuracy
or broad qualification.

## Independent audit and ownership

`audit_beta_wall.py` verifies all frozen inputs, complete arm inventory, initial
coordinates and every retained artifact. Whole original systems are compared after
removing only the declared research forces; the only phase-2 physical-system
change is the prescribed original-restraint reduction. Masses, constraints,
original force-field terms and coordinate identities are not altered.

Each phase is checked using an **actual-frame abstract heavy-atom geometry view**,
with captured external bonds retained. Raw and PDB-precision-rounded *whole*
coordinates are audited separately. The existing measured-χ bond/angle/clash
checks, independent RDKit target CIP and global Cα/THR-Cβ signed-volume checks
remain separate from the experimental wall. No source-file materialization or
observation-preserving-output claim is made by this diagnostic view.

**Original observations are free to move under the original restraints**, and
their full-coordinate movement is retained; none are restored afterward. H is
present in every physical minimization, but this common heavy-geometry validator
does **not** qualify H geometry. Finite endpoint passes do not certify every
intermediate iteration, general convergence, or immunity to a stronger drive.

## Retention and implications

Private evidence is under `tmp/pdbfixer-research/`: `beta-wall-plan-01.json`,
`beta-wall-run-01/` and original `beta-wall-audit-01.json`. All failed raw endpoint
coordinates remain private and retained. The strengthened full-system auditor
preserves the same numerical results; it runs no new minimization. Reproduction:

```bash
PYTHONPATH=src:. python scripts/audit_beta_wall.py \
  --plan tmp/pdbfixer-research/beta-wall-plan-01.json \
  --run tmp/pdbfixer-research/beta-wall-run-01 \
  --output <fresh-evidence.json>
```

The bounded conclusion is architectural: a stereochemical restraint needs an
independent fail-closed final **geometry** audit, not merely a sign check. Already
wrong-handed inputs also require explicit refusal or separately qualified repair;
this test does not repair them. It does not justify shipping this wall, applying
it to VAL, or generalizing to ILE/PTMs/other chemistry.

No production code/default/dependency/public capability, installed libraries,
scientific acceptance thresholds, commit or push changed. The original unforced
THR2 inversion and all previous negative/nonexact replays remain unresolved and
retained, rather than replaced with the easier artificial-drive result.

Verification: **815 focused tests passed in 27.66 s**. Independent saved-only
evidence regeneration is byte-identical. Targeted Ruff and the 17-file completion-
package mypy check pass; this is not full repository/platform qualification.

Subsequent [fixed-χ source-geometry controls](start-geometry.md) test upstream
reconstruction of only provenance-owned generated atoms, rather than a stronger
wall. Both rebuild arms pass final geometry12/12 in that bounded controlled
procedure; original THR2 causality and initial geometry failures remain distinct.
