# Private public-API post-H adapter and independent 1VCU continuation

## Scope and frozen enrollment

Continuation after `1c8140a` under the user's autonomous-work instruction.
The [prospective contract](../plans/post-hydrogen-public-api.md) and
[two-arm 1VCU plan](post-hydrogen-1vcu-plan.json) are frozen at **`94b8098`**,
before either scientific worker executes. Work is isolated on
`research/post-hydrogen-public-api`; production, installed OpenMM and the earlier
variant-1/variant-2 sources/evidence remain unchanged.

The independently activated source is retained 1VCU A/SER0/HA after forcefield-
free prepare fallback. Its final prepared PDB and DatRecord contain 13 protonation
variants, so ordinary minimize strips protein H and adds them again. Simply
patching the prepared HA does not establish that a correction survives this
boundary. The new pair therefore compares an observational baseline with a
guarded prepared HA plus a guard after each successful public H-add call.
Both retain all ordinary variant/ligand/chirality policies and --no-solvent.
This constructs new matched whole-complex experiments; it is not an exact replay
of a previously captured whole-system 1VCU Context (none was available).

## Transactional adapter and engineering evidence

`scripts/post_hydrogen_public_api.py` delegates to the original public
`Modeller.addHydrogens` on an independent working model, resolves exact before/
after identities and graph, runs the unchanged geometry rule and commits last.
An API error, typed guard refusal or diagnostic-write failure leaves the caller's
model untouched. Actual native FF template errors still propagate to the normal
pipeline fallback; guard/ownership refusals use the geometry-error escape, not
ValueError/KeyError. No global LocalEnergyMinimizer interception is introduced.
Scoped public-method instrumentation exists only inside isolated research
workers and is restored on exit, not installed as production behavior.

The caller supplies original observed-H lineage explicitly. Missing/stripped
protected H refuse before the API call; they are never silently classified new.
In this declared 1VCU source, the original tracked PDB has **zero H/D**. All
pre-call existing H remain immutable to the guard regardless of their upstream
origin. Atomic identity preserves chain case, insertion codes and component
occurrence, with shared variant aliases. Existing bond graph **and bond metadata**
are checked. Nonempty object-keyed residueTemplates are explicitly unsupported.

API rounding and guard motion are separate: before→raw API uses the predeclared
native transport bound (exact when Reference explicitly requested); raw→guarded
preserves every atom outside the new-H move mask exactly. No restoration to
pre-API coordinates hides native single-precision rounding. The API's return
value, selected variants and caller RNG restoration remain reported.

Twenty-one pre-execution engineering controls pass: transactional rollback,
observations/protected H, deletion/alias/graph/bond-order failures, malformed
returns, exact Reference transport, return value, RNG neutrality, PDB coordinate-
only handoff and finite-difference constraint derivatives. With the two existing
multi-path accounting controls, the focused set passed 23 tests.

An engineering test caught a real adapter-construction hazard **before any native
qualification**: naive deepcopy(Topology) on the installed runtime left copied
Bond endpoints referring to original atoms. The adapter now reconstructs the
working topology through public constructors, explicitly rebinds every bond and
preserves chain/residue/atom metadata and periodic box vectors. A regression
checks both endpoint ownership and snapshot equality. This is a bounded local
finding, not a claim about every OpenMM/Python version or the earlier trajectories.

## Native and audit boundaries

Two workers run once, seed label 1 before each public H call, ordinary minimize
defaults, owned ligand caches, complete actual FF/fallback and API traces.
Baseline parameters are copied to an isolated candidate cache; actual whole
Systems/anchors and raw API returns are compared, not assumed identical. GPU
sampling can remain numerically non-bit-identical despite the same saved RNG.
Never relabel a strip-and-splice output as a full-system ligand qualification.

Kept waters can add distance constraints. The independent audit reconstructs
their reported harmonic penalty/gradient in addition to physical System forces,
using the prospectively fixed single-precision numerical envelopes. This does
not add or retune native forces. Relevant exact-revision source references are
[CommonMinimizeKernel](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/common/src/CommonMinimizeKernel.cpp)
and [ReferenceMinimize](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/platforms/reference/src/SimTKReference/ReferenceMinimize.cpp).
Source/version inspection is not independent attestation of installed GPU
binaries. Constraint projection and bounded numerical agreement must not be
reported as exact terminal equality or certified constrained stationarity.

Private PDBs, native States/checkpoints, callback arrays, RNG states and parameter
caches belong under the ignored project run directory only. Public records are
limited to identities, hashes, scalar findings and actual outcome statuses.
The adapter is not public/production-ready: broader observations/chemistries,
formal-charge-rich API inputs, multi-model/altloc transport and custom
residueTemplates remain outside this bounded live qualification.

## Actual paired result: guarded integration refused

The pair executed once under `94b8098` with no retries. The
[complete first audit](post-hydrogen-1vcu-audit.json) retains **baseline completed /
guarded execution-error**. Baseline ran both whole phases with DAN and EPE present
and wrote 6,205 atoms with zero detected D-Cα. Guarded stopped after the native
H call returned, before any whole phase or output publication. Its exception
was `PostHGuardError`, not a template ValueError/KeyError, and triggered no FF
fallback. Both actual H placements used provided AMBER19/GAFF-enabled physics,
not the earlier prepare-time forcefield-free route. This new native boundary is
not interchangeable with the old fallback snapshot.

| Arm | Underlying native H call | Caller adapter commit | Strict guard refusal |
|---|---|---|---|
| baseline | returned | raw return committed | GLU228/HA radius **2.137 Å** |
| guarded | returned | **no commit** | PRO379/HA radius **1.651 Å** |

The original worker `h_calls.returned` flag records whether its **wrapper**
returned; it is false for guarded even though the underlying API returned.
Likewise the first audit's legacy `api_returned` field mirrors that wrapper flag.
Interpret completion from the retained raw-return journal/stage separately; do
not call the typed refusal a failed native H placement or invent an API retry.
The private raw journals retain selected variants and the entire returned model.

Raw SER0/HA remains sign-incompatible in both arms, with C–H lengths about
1.133/1.123 Å. The current strict rule refuses atomically because another
supported center is outside the 1.5 Å radius admission. No threshold was changed
and no partial projection was committed. Source/DatRecord handoff and the
prepared HA-only patch were retained, but prepare correction alone is lost by
the ordinary variant-driven strip/readd; the proposed early public hook did not
reach a validated whole-system guarded state.

### Baseline completion is not clean qualification

Baseline printed an astronomical starting energy of **2.345410118750321e13
kJ/mol**. Existing production postflight did not eliminate all starting defects.
Both whole calls returned, with 2,282/2,058 callbacks (constraint-restraint outer
iterations can exceed the per-inner-run maxIterations label). There were **213
distance constraints**. Final constraint projection changed last-callback
coordinates by ~2.93e-6/3.44e-6 nm, kept explicitly, not labelled exact return.
Final relative constraint errors were ~2.38e-6/3.13e-6; constrained stationarity
is not certified by this audit.

The direct same-GPU physical audit failed its fixed envelopes at **26 phase-1
callbacks**, with zero failed phase-2 callbacks. Keep that first audit failed:
early reported forces exceed the GPU fixed-point range, and the exact-revision
source documents a CPU overflow fallback. The separately
[frozen read-only attribution protocol](../plans/post-hydrogen-overflow-readback.md)
tests CPU and Reference for all 26 callbacks with the same original bounds; it
does not rerun or repair the native trajectory or turn the original test green.

The [secondary readback plan](post-hydrogen-overflow-readback-plan.json) was
frozen at `9a0aeba`; [all 52 CPU/Reference readbacks](post-hydrogen-overflow-readback-evidence.json)
executed once, without budget nonexecutions or new optimization. CPU system
energies agree with the reporter within **7.63e-6 kJ/mol**, consistent with the
documented overflow-fallback hypothesis. The maximum gradient residual drops
from ~4.03e12 on direct GPU readback to ~2.68e5 on CPU (~5.69e6 on Reference).
Nevertheless **0/26 callbacks pass the complete original component-wise
envelopes on either secondary backend**. No bound was enlarged. Native gradients
exceed fixed-point force capacity in these selected callbacks, but energy
agreement/source inspection is not full gradient attribution or proof of an
internal event. Keep physical attestation **incomplete**, not repaired/green;
further numerical localization would be a separate read-only research stage.

Even after baseline whole relaxation, SER0's HA is spatially incompatible; its
heavy Cα volume stays positive but falls from ~0.925 Å³ to ~0.180 Å³. Zero D-Cα
and a normal ~1.091 Å HA bond are therefore not sufficient to call its local
tetrahedral geometry correct. Preserve this distinction from chemical CIP or a
claim of successful comprehensive geometry regularization.

## Offline ordering control: repair first does not qualify both states

Follow-up [contract](../plans/post-hydrogen-after-repair.md) and
[plan](post-hydrogen-after-repair-plan.json) were frozen at **`8c458ea`** before
applying ordinary repair to either retained raw return. Both states were enrolled,
not just baseline. The [complete outcomes](post-hydrogen-after-repair-evidence.json)
used zero Contexts, minimizations, H samples or changed thresholds.

- Ordinary repair returned counts 362/343; heavy coordinates and both Python/
  NumPy RNG states were unchanged. These are the existing method's counts,
  not a new sampling cohort.
- Shadow baseline ordinary repair matches the actual captured whole-System
  start within the predeclared float32 transport bound (maximum component error
  ~2.074e-7 nm). It fixes GLU228/HA to 1.09 Å; the unchanged guard can then propose
  SER0/HA and the HAs at residues 283/286.
- Shadow guarded ordinary repair leaves PRO379/HA at **1.650683 Å** and the
  unchanged strict guard still refuses. Retrospective inspection of that retained
  center finds positive heavy/four-ligand volumes (~1.372/~6.558 Å³): it is
  **spatially compatible**, but outside radius admission. This sign inspection is
  descriptive, not a retrofitted pass or proof of acceptable C–H length.

The predeclared **both-states-admitted native gate is false**. No new native
continuation was launched from only the easier baseline state. Moving the guard
after ordinary repair is not a qualified live solution to this failed pair.
Variant 1, variant 2 and their earlier positive/negative results remain intact.

## Decision boundary

The current conservative policy validates distance admission at every supported
new-H center, including compatible centers it would never move. That policy
blocks this 1VCU path on the compatible but elongated PRO379/HA. Continuing to
native testing now requires an explicit policy choice, not a silent wider cutoff:
keep global strict refusal, or prospectively study a new local policy in which
compatible centers are exact no-ops with visible out-of-scope warnings while all
actually moved H retain the original strict distance/sign/contact bounds.
No such policy change has been implemented or enabled here. Earlier tests,
thresholds, native forces and failed outcomes are unchanged; public readiness
remains false.

Final engineering validation is **3,269 fast tests passed**, three optional
Torch/Biotite skips and 33 slow deselections. All 24 new focused controls,
Ruff/strict-target mypy, generated CLI/GUI specs and agent maps pass. The
[validation/quality record](post-hydrogen-public-api-validation.json) deliberately
marks scientific qualification **failed/blocked** despite those engineering
passes. Original native failures, direct-GPU residuals and unsuccessful secondary
gradient comparisons remain bound alongside the successful baseline output.
