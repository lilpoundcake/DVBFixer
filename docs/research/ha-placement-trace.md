# The retained inward GLN16 HA appears during OpenMM's H-only minimization

Status: **original-stage provenance**, not a new repair or original THR2 replay.
Follows the [one-HA controlled replay](ha-start-control.md): changing only the
generated HA start prevents the later GLN16 heavy-Cα flattening for one retained
1SHG start. This study asks **where that inward HA came from**, without changing
the original H-placement algorithm, seed or minimizer settings.

## Faithful capture of all original H-placement calls

[Complete evidence](ha-placement-trace-evidence.json) retains all **24 enrolled
cells**, all completed with no retries, timeouts or budget nonexecutions. These
remain four masks × two original H seeds × three source-geometry arms of **one
protein**, not 24 independent protein examples.

Workers execute the preserved historical source/worker from the detached baseline,
with byte-identical copies of the corresponding original prepared PDB and `.dat`.
They observe the native OpenMM **8.5.2** `Modeller.addHydrogens` internal Context
immediately before/after its existing minimizer call, forwarding that call
unchanged: **tolerance 1.0 kJ/mol/nm, maximum 50 iterations, Reference**. Existing
atoms have zero mass in that H-placement System; newly added H are mobile.
An expected boundary stop prevents the first whole-system minimization from
executing. No H method, force, strength, seed, iteration budget or production
source is changed. Raw worker logs retain this intentional early-stop boundary.

The frozen plan is private at `tmp/pdbfixer-research/ha-placement-plan-01.json`,
SHA256 `92e057f80deb827070721c2832f4de0c96aca52ee214dedd0ef6227342be8919`.
It binds the preserved scientific source and research scripts, installed Modeller
source, all saved inputs and the tracing driver. Budgets were 900 seconds per
cell and 3,600 seconds total. Neither the baseline tree nor the saved historical
inputs are edited. Raw Context/System records remain private below
`tmp/pdbfixer-research/ha-placement-run-01/`; those paths do not transfer via Git.

## Independent reproduction and immobility checks

For **all 24** cells, the independent saved-only auditor confirms:

- all existing input atoms remain coordinate-exact throughout H placement;
- original particles have zero mass, added particles are mobile H;
- every H-stage endpoint exactly matches the original saved H endpoint;
- the subsequent first whole-system start exactly matches its original saved
  coordinates **and** serialized System parameters;
- full-system minimization has not executed, and hook scope is restored.

The maximum endpoint/start difference is **0.0 Å** for both heavy atoms and H.
Thus the observed stage is a faithful replay, not an alternative H construction
or an approximate attempt that happens to have similar angles. Separate synthetic
tests check forwarding of the real H minimizer, the intentional stop, restoration
on failure and independent ownership/geometry refusals. Eight focused tests and
targeted Ruff pass. The audit/publication regenerates byte-identically. The prior
HA-control full fast suite passed 2,830 tests; this additional eight-test increment
has focused acceptance, not an inferred new full-suite result.

## The problematic HA crosses into an inward endpoint during this stage

For THR24 / full-generated / original H seed 1:

| GLN16 measurement | Before native H-only minimization | After |
|---|---:|---:|
| N–CA–HA | 117.734° | **64.207°** |
| HA–CA–CB | 91.407° | **72.748°** |
| HA–CA–C | 114.548° | **78.228°** |
| CA–HA distance | 1.000000 Å | 1.076222 Å |
| Heavy Cα signed volume | 2.528406 Å³ | **2.528406 Å³** |
| Existing heavy-atom maximum movement | — | **0.0 Å** |

All three angles are initially obtuse, including one near the hemisphere boundary;
all three are acute afterward. The **472 existing heavy atoms are immobile**;
483 H are added. Therefore the retained inward HA does **not** begin as an
all-acute random-placement endpoint: it appears during this native H-only
minimization. The heavy Cα itself does **not** flatten during this stage, because
the heavy atoms cannot move. Its later flattening occurs during the separate
whole-system phases documented in the paired replay.

Across all 24 cells, initial triads are 12 all-obtuse and 12 mixed/boundary, with
zero all-acute. After H placement they are 18 all-obtuse, five mixed/boundary and
one all-acute. These are geometric descriptors, not new gates, population rates
or successful general H sampling. All cells and descriptors remain published.

## Interpretation and next boundary

The evidence now links three distinct observations for one retained 1SHG case:
the native H-only relaxation produces the inward HA; H-containing angle forces
oppose restoration at the later flattened endpoint; and a matched one-HA start
change prevents that later failure under an exactly reproduced full-system
baseline. This localizes a mechanism, but does not identify a particular internal
line-search step, prove convergence of the 50-iteration H stage or isolate which
H-stage force/interaction drives HA inward.

The subsequent [native minimizer path readback](ha-minimizer-path.md) now records
all 1,200 callbacks with exact original endpoints: the known HA is already inward
by the first reported step, after a large parent-bond stretch. Its enormous
initial inward-axis force is nonbonded-dominated, while the initial tangent force
is outward. This is not ordinary bond-preserving rotation or a physical path.
A particular nonbonded partner remains unidentified and needs separate validated
pair decomposition before making that causal attribution.

The next experiment declared at this boundary was a separately frozen intermediate-step
or force-component readback under the same H-stage System, retaining all cases.
Do not silently increase iterations or install the outward construction as a fix.
The remaining THR target flattening and VAL naming failures, original unforced
1CRN THR2 inversion, broader protein/platform qualification and empirical-library
production/license gate remain separate and unresolved. Production is unchanged.

Code: `scripts/research_ha_placement_trace.py`,
`scripts/audit_ha_placement_trace.py`.
