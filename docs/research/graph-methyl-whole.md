# Source-diverse downstream control: native ASP49 flattening persists in selected arm

Status: **case-specific native-phase research**, not production repair or complete
structure qualification. Follows the [graph-methyl H-stage control](methyl-diversity-control.md).
All source refusals and earlier negative results remain. No production code,
default, threshold, force or chirality safeguard changed.

## Independent baseline and exact transport

[Complete centre-level evidence](graph-methyl-whole-centres-evidence.json) reports **48 continuation
cells: 30 completed and 18 refused-source**, with no retries, execution blockers,
timeouts or budget nonexecutions. Ten original saved-State native baselines were
executed separately before those continuations; six source/seed refusals remain
in that baseline enrollment. Five protein views × two seeds are ten contexts,
not ten independent proteins.

The unchanged production H repair is applied once under the actual internal H
labels, then positions are mapped through the captured topology-object index
bijection into the full labels. All ten original starts preflight **exactly in
nm**. The original full-System/State/integrator and forces/anchors are reused.
Independent baselines call native `Simulation.minimizeEnergy`; continuations use
the corresponding direct minimizer with the same **Reference / no-solvent /
tolerance 10 kJ/mol/nm / maximum 1,000 iterations per phase** settings. Phase 2
divides only the original positional restraint weights by ten, retaining anchors.
All ten unchanged starts and **both endpoints reproduce their native baselines
exactly**, every atom in nm. Every continuation's serialized System agrees with
its corresponding baseline phase.

These are **two native phases before production's emergency chirality recovery**,
not the complete production pipeline. No emergency reflection, observation
restoration, candidate bundle or PDB output is introduced. Production's final
zero-D safeguard remains unchanged. Budgets are 1,800 seconds per cell and 14,400
seconds per baseline/continuation batch, without retry or favorable-case replacement.

Plan SHA256:
`cb35b44b30346854d3b4a216271120d0fddec81959085f74fa5abb0ba1a680eb`.
Private plan, baseline and run are `graph-methyl-whole-plan-01.json`,
`graph-methyl-whole-baseline-01/` and `graph-methyl-whole-run-01/` below
`tmp/pdbfixer-research/`; those paths do not transfer through Git.

## All 90 continuation states: continuous geometry, not a new acceptance gate

The auditor preserves raw and PDB-rounded **continuous signed volumes for every
Cα and THR-Cβ centre**, plus the unchanged 0.5 Å³ shape descriptor. It also reports
all-acute HA, H-bond/angle strain, short H contacts and original heavy-atom motion.
An output-reporting extension adds all continuous volumes without changing inputs,
coordinates, thresholds or the prior numerical evidence; the
[initial compact audit](graph-methyl-whole-evidence.json) is also preserved.
There is no masked-target or native-RMSD score in this
cohort. Do not compare its global-only descriptor to the earlier 1SHG target-heavy
gate as if they were the same scientific acceptance test.

| Final descriptor | Unchanged | Graph zero | Graph selected |
|---|---:|---:|---:|
| Executed contexts | 10 | 10 | 10 |
| Negative Cα centres | 0 | 0 | 0 |
| Cα centres at/below 0.5 Å³ | 0 | 0 | **1** |
| THR-Cβ centres at/below 0.5 Å³ | 0 | 0 | 0 |
| All-acute HA | 0 | 0 | 0 |

The selected-arm shape finding is **1TIM/A/ASP49, seed 1**, raw **0.498283 Å³**,
rounded **0.497121 Å³**. However, unchanged ends at **0.500677 Å³** (rounded
0.500967), just above the pre-existing threshold. This is **similar substantial
flattening in both arms**, not a large effect established by their opposite binary
classification. No threshold is retuned to hide or magnify that distinction.

## 1TIM ASP49: H-stage inward HA is followed by observed-heavy flattening

All arms start from exactly the same heavy volume. Only generated methyl-H
initialization differs upstream:

| ASP49 Cα signed volume | Unchanged | Zero | Selected |
|---|---:|---:|---:|
| Initial whole-system start | 2.575429 Å³ | 2.575429 Å³ | 2.575429 Å³ |
| Phase 1 endpoint | 1.065215 Å³ | 2.577428 Å³ | 1.060820 Å³ |
| Phase 2 endpoint | **0.500677 Å³** | **2.609578 Å³** | **0.498283 Å³** |

The original inward H-only endpoints in unchanged/selected do not remain formally
all-acute after full relaxation, but the final HA angles are still distorted:

| Final ASP49 HA measurement | Unchanged | Zero | Selected |
|---|---:|---:|---:|
| N–CA–HA | 80.871° | 109.209° | 80.973° |
| C–CA–HA | 93.586° | 108.454° | 93.574° |
| CB–CA–HA | 86.077° | 109.152° | 86.044° |

Thus a zero final all-acute count alone misses this failure mode. The fixed-zero
graph arm avoids this specific outcome; selected initial global energy does not.
This is an independent protein-context example of native H-placement-associated
heavy flattening, not yet a pair-specific causal attribution or an explanation of
the original unforced 1CRN THR2 Cβ event. Earlier fixed-phase 1SHG TYR15/negative-Cα
regressions still rule out calling zero phase a safe general repair.

## Original observations move and H descriptors are incomplete

For 1TIM seed 1, all **1,870 original heavy observations** move during whole-system
relaxation; their RMS/max movement from the prepared-conditioning start is:

| Movement | Unchanged | Zero | Selected |
|---|---:|---:|---:|
| RMS | 0.288697 Å | 0.285958 Å | 0.289334 Å |
| Maximum | 1.892546 Å | 1.880936 Å | 1.932975 Å |

This is not observation-immutable reconstruction. Across the ten contexts,
median maximum final H-bond deviations are **0.03114 / 0.03218 / 0.03166 Å**, and
maximum H-angle deviations are **28.629 / 11.852 / 28.527°**. No final <1 Å
H-contact descriptor is triggered, but that limited distance check is not complete
steric/chemical acceptance. No full-complex, convergence, template-superiority or
general repair qualification follows from these results.

## Validation and next boundary

Five focused native-baseline/transport/audit tests pass, including end-to-end
saved-only audit and parameter-drift refusal with retained partial evidence.
Ruff, generated CLI/GUI and agent-map checks pass. The fast suite through native
baseline integration passes **2,969 tests**, three optional Torch/Biotite skips and
33 slow deselected; later audit tests are separately exercised. Independent
evidence publications regenerate byte-identically.

The subsequent [saved-only force/pair attribution](graph-methyl-ASP49-attribution.md)
now covers all six matched paths and 1,172,184 pair observations. Seed-1 inward
crossings are delayed until callbacks 11/12 and are preceded by strong LEU48 HB2
methylene-H repulsion, not the initial methyl-H impulse seen in 1SHG. Selected
seed 2 also has a retained transient inward path despite an outward endpoint.
This is exact attribution, not a controlled pair-causality experiment.

The next step declared here was a separately frozen **saved-only** attribution
of all three arms and both existing seeds without more minimizer/phase/seed search.
Retain adverse
cases and source-bound physics; do not integrate this initializer into production.
Original THR2 causality and empirical-library production provenance/license/
approval remain separate unresolved gates.

Code: `scripts/research_graph_methyl_whole.py`,
`scripts/audit_graph_methyl_whole.py`.
