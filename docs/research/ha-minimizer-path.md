# Native H minimizer: the retained HA crosses inward by the first reported step

Status: **original-stage optimizer-path observation and static force readback**,
not a production repair, physical trajectory or original 1CRN THR2 explanation.
Follows the [original H-stage boundary trace](ha-placement-trace.md). Production
backends, thresholds, force parameters and public completion capabilities remain
unchanged.

## Complete, unchanged native replay

[Complete evidence](ha-minimizer-path-evidence.json) contains **24/24 completed
cells, 1,200 callbacks and 1,248 measured nodes** (initial + every callback + final).
These are four masks, two original H seeds and three source-geometry arms on
**one protein, 1SHG**, not 24 independent proteins. No retries, execution blockers,
timeouts or budget nonexecutions occurred.

A read-only native `MinimizationReporter` is attached at the original H-only
minimizer boundary, against the same preserved historical source and original
H-placement seeds. It copies coordinates, objective gradients and statistics and
always returns `False`. The original Reference call remains tolerance
**1.0 kJ/mol/nm**, maximum **50 iterations**. No alternative H placement, new seed,
force, restraint, step limit, integration or whole-system minimization is used.
All captured Systems have zero constraints, no PBC and no virtual sites. A
512-callback storage bound would retain an execution failure, not report early
convergence; every observed cell produced exactly 50 callbacks.

The frozen plan is `tmp/pdbfixer-research/ha-minimizer-path-plan-01.json`, SHA256
`4d159a5be0823c18577e0cb821c9a534361d6a41c28cb572617369137c0cfb2e`.
Raw path coordinates/gradients are private under
`tmp/pdbfixer-research/ha-minimizer-path-run-01/`, with 900-second per-cell and
3,600-second total budgets. Those paths do not transfer through Git.

The native Context is observed directly, rather than reconstructing a start by
round-tripping saved Å coordinates into nm. The raw start/final nm arrays match
their audited Å snapshots exactly. For all 24 cells, adding the reporter leaves
the original start and final H endpoint **coordinate-exact**, and all existing
heavy coordinates remain exact at **every callback**. The last callback equals
the returned final coordinates. Original H-stage System parameters also match.

## Force readback and interpretation

At every measured node, an independent saved-only Reference Context measures each
original force component and their combined force/energy. Measurement groups are
the only System metadata changed; restoring them must reproduce the serialized
original System. No minimization, integration or coordinate perturbation occurs
in this audit. Force/energy components must sum to the combined values.

Let **u** be the fixed outward direction opposite the sum of the three unit
CA-to-N/C/CB heavy-bond vectors. We report HA's signed outward distance,
its normalized outward cosine and the Cartesian force component along **u**.
A separate **sphere-tangent** projection points toward increasing outward cosine
at fixed current CA–HA radius. A singular tangent is explicitly undefined, not
assigned an arbitrary direction. These projections are not interchangeable:
bond stretching can dominate the first while contributing zero to the second.

At all **1,200 native callbacks**, measured physical HA force equals the negative
recorded objective gradient exactly: maximum component difference **0.0
kJ/mol/nm**. Native energies match the readback, and constraint-restraint energy
and constraint error are zero. Thus the force data are checked against the
actual callback state, not evaluated at a different geometry or conflated with a
constraint penalty.

## The known THR24 / full-generated / seed-1 case

The following is an explanatory view; all nodes and components remain published.
Callback ordinals/native indices are zero-based here and retained separately.

| Node | Outward cosine | Signed outward distance | CA–HA radius |
|---|---:|---:|---:|
| Initial | +0.962637 | +0.962637 Å | 1.000000 Å |
| Callback 0 | **−0.912389** | **−2.286774 Å** | **2.506358 Å** |
| Callback 1 | −0.830338 | −3.679187 Å | **4.430951 Å** |
| Callback 17 | −0.971254 | −1.043797 Å | 1.074691 Å |
| Callback 49 / final | −0.994893 | −1.070726 Å | 1.076222 Å |

The first negative outward sample is already callback **0**. The first callback
where **all three** HA angles are acute is **17**; a negative outward cosine
alone is not the same descriptor. Later callbacks briefly return to positive
outward cosine before the final inward basin, so the path is not described as a
single monotonic flip. No intermediate line-search evaluation is claimed.

At the **initial state**, outward-axis force projections (kJ/mol/nm) are:

| Component | Outward-axis projection | Sphere-tangent projection |
|---|---:|---:|
| Harmonic bonds | +2464.937 | ≈0 |
| Harmonic angles | +464.623 | +1715.780 |
| Periodic torsions | +0.036 | +0.132 |
| Nonbonded | **−56105.474** | **+3808.931** |
| **Combined** | **−53175.879** | **+5524.842** |

The enormous initial inward-axis force is dominated by the **nonbonded**
component, not a flattening heavy-atom angle force. Yet its initial tangent
projection is outward: this is not evidence of an initially inward *rotational*
force. By the first reported optimizer step HA lies on the opposite side with a
strongly stretched parent bond. This is a large unconstrained Cartesian optimizer
displacement, not ordinary bond-preserving rotation or molecular dynamics.
Force projections alone do not identify the full L-BFGS displacement or prove a
specific line-search/bond-crossing sequence.

The initial whole-H-System energy is 145,674.319 kJ/mol and falls to 5,464.180
kJ/mol at the final callback. Reaching the declared 50-iteration endpoint is
**not a convergence certificate**. At that endpoint the sphere-tangent components
include angles −741.258, nonbonded +697.866 and torsions +47.328 kJ/mol/nm,
leaving +3.936 combined; this cancellation is not a stable-basin qualification.

## Complete cohort and remaining boundary

Only the two full-generated THR24 starts have a negative-outward callback, both
first observed at ordinal 0. Seed 2 also makes a large early displacement
(CA–HA radius 5.373800 Å at callback 0), but returns to positive outward cosine
and never has an all-acute callback. The known seed-1 endpoint therefore must not
be generalized to every large early step or every original H seed. The heavy Cα
remains fixed throughout this stage; its later flattening is a different event
established by the [paired one-H intervention](ha-start-control.md).

The largest recorded CA–HA radii are 4.430951/5.697322 Å for those two seeds;
the other 22 cells never exceed 1.210774 Å at the measured nodes. Some raw THR24
starts also have negative initial outward-axis force but never cross inward, so
its sign alone is not a predictive gate. Recorded energies decrease throughout
each of the 24 paths; that does not establish local force convergence or acceptable
geometry at the fixed iteration limit.

This localizes a large early nonbonded-dominated force and observed optimizer
displacement, but does **not** identify which interacting particle/pair creates
it. A next static readback should separately freeze and validate exact HA-touched
nonbonded pair sums against the unmodified force before naming a clash partner.
Do not infer a particular pair from spatial proximity alone, suppress failures,
increase iterations, or promote an outward-H construction to production.
Original THR2 causality, remaining THR/VAL failures, broader protein/platform
qualification and the empirical-table production/license gate remain open.

Validation: five native-reporter tests and six independent frame/force/gradient
tests pass; the combined angle/H-stage set passes 59 tests. Targeted Ruff passes.
Independent publication regenerates byte-identically. The full fast suite passes
**2,849 tests**, three optional Torch/Biotite skips and 33 slow tests deselected;
agent-map and generated CLI/GUI checks pass. This is not slow-suite or broader
scientific/platform qualification.

Code: `scripts/research_ha_minimizer_path.py`,
`scripts/audit_ha_minimizer_path.py`.
