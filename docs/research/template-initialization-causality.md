# Template initialization controls: D-Cα precedes optimization, proper rotation is insufficient

Status: **retrospective case-specific mechanism controls**, not an accuracy cohort,
a production repair or a qualification of a new reconstruction backend. Extends
the [earlier PDBFixer mechanism study](pdbfixer-reconstruction.md) to endpoints from
the [pre-download new-system comparison](rosetta-independent-controls.md).

**Main result:** both confirmed-D whole-sidechain cases already place the new CB
on the D side during construction, before minimization. Disabling all relaxation
does not remove it. A determinant-corrected proper fit preserves the complete
template's L geometry but still produces a D-side **hybrid** when its sidechain is
joined to unchanged observed backbone coordinates. Correcting `det(R)` alone is
therefore not a sufficient remedy for these inputs.

## Frozen factorial and observational noninterference

Before the new stage runs, `research_template_causality.py plan` bound the original
masked bytes, untraced seeded proposals, scientific sources/resources and this
four-case × two-overlay × two-relaxation factorial:

| Case | Deliberate role |
|---|---|
| 2CI2 I/VAL28 whole sidechain | Previously confirmed D-Cα endpoint |
| 2CI2 I/THR55 whole sidechain | Previously confirmed D-Cα endpoint |
| 1R69 A/THR27 whole sidechain | Previously small-volume, CA-CIP-unassigned endpoint; not confirmed D |
| 2CI2 I/THR22 whole sidechain | First whole-sidechain case in original enrollment order; not selected by success |

Overlay modes are upstream stock fit and determinant-corrected proper fit.
Relaxation modes are normal upstream behavior and **none**, which bypasses both
local minimization and clash-escape integrator steps. The latter is an ablation,
not a geometry-valid output or a suggested user workflow. No anchor set, force
field, geometry cutoff, source atom, protonation or native coordinate is changed.
Both overlay modes retain the original observed-anchor set; no backbone-only
selection intervention is silently introduced.

The unchanged seeded helper owns every rebuild call, with its original maximum
of five attempts. Hooks observe actual heavy/all-atom construction boundaries,
before/after minimization, reporter iterates, dynamics-call boundaries and returned
seed states. No extra construction preview is called. Hooks are process-local,
restored after each arm; installed source files and production code are untouched.

**All four stock/normal arms reproduce every untraced seed state exactly**:
sixteen complete atom/element/coordinate inventories match the earlier comparison.
Thus the stock observations are not reinterpreted after a tracer-induced change.
Such a mismatch would have been an instrumentation blocker, not a scientific
rejection or permission to adjust the old expected output.

## Executed trace and independent audit

[Independent stage/template readback](template-initialization-causality-evidence.json)
validates:

- Four cases, sixteen factorial arms, **60 returned helper seed states**.
- **300 complete stage snapshots** with exact identity and observed-coordinate
  checks, and **75 target-coordinate minimizer reporter iterates** with independently
  measured volumes and observed target-backbone checks.
- Every overlay determinant and canonical template anchor ordering/resource hash.
- Complete-template transformed geometry versus actual fixed-backbone hybrid.
- No-relaxation construction/final equality, restoration of process-local hooks,
  original untraced baseline replay, and unchanged installed scientific sources.

No clash-escape dynamics step is called in any arm: **zero dynamics calls**.
Returned payloads are identical across retries within each arm on this pinned
Reference control. The five seeds are therefore not five independent conformational
samples; these particular endpoints do not arise from stochastic Langevin escape.
Do not generalize that result to other inputs/platforms or all PDBFixer rebuilds.

All observed source atoms remain fixed. The missing CB has no observed input
coordinate: this is wrong-side placement of a **new** atom, not movement/reflection
of an observed sidechain. The tracer does not restore observations post hoc.

## Construction versus minimization

First-seed signed Cα volumes, Å³; every later returned seed in the same arm is
identical. Final states here are raw diagnostic endpoints, not accepted outputs.

| Case | Stock construction → normal return | Proper construction → normal return |
|---|---:|---:|
| 2CI2 VAL28 | −0.38863 → −0.62739 | −1.55533 → −1.53251 |
| 2CI2 THR55 | −0.62022 → −0.55516 | −1.27528 → −0.88048 |
| 1R69 THR27 | +0.06179 → −0.08247 | −1.08344 → −0.48755 |
| 2CI2 THR22 | +3.09377 → +3.09377 | +3.09377 → +3.09377 |

Both primary cases have CA **R already at construction** under stock fitting.
Minimization changes their geometry but is not necessary to create these wrong-side
placements: stock/no-relaxation returns exactly the constructed target coordinates.
The original final-D observations therefore cannot be attributed solely to
optimization in these two cases.

The 1R69 stock case is different: its tiny volume changes sign during minimization,
but CA CIP remains unassigned at both boundaries. This is an optimizer-stage
change in a small-volume collapse, **not demonstrated L→D inversion**. The proper
intervention instead starts with CA R on this case; it is not a rescue. Keep
signed-volume sign, degeneracy and independent CIP assignment separate.

The identity-chosen THR22 control stays CA S/CB R with positive volume in all four
arms. The tracer/proper intervention does not simply label every case inverted.

## Why a proper rigid transform does not protect the hybrid

The stock fit permits `det(R)=-1` in the two primary cases and the small-volume
case; it uses `det(R)=+1` for THR22. A proper fit enforces a rotation, not reflection.
Independent canonical resource measurements show L templates: VAL volume
**2.407856 Å³**, THR **2.457112 Å³**.

For each recorded proper transformation, the auditor transforms **every template
atom** and verifies both the determinant volume law and CA S (THR CB R). It then
replaces only the observed target atom coordinates with the unchanged actual
input, reproducing the exact constructed hybrid:

| Proper construction | Fully transformed template volume | Actual hybrid volume | Template / hybrid CA CIP |
|---|---:|---:|---|
| VAL28 | +2.407856 | −1.555334 | S / R |
| THR55 | +2.457112 | −1.275285 | S / R |
| THR27 | +2.457112 | −1.083445 | S / R |
| THR22 | +2.457112 | +3.093773 | S / S |

The determinant guarantee applies to one rigidly transformed set of coordinates.
It does **not** apply after mixing transformed sidechain coordinates with a
different fixed observed backbone. This is directly verified reconstruction of
the saved stages, not a native-RMSD selection or an assertion that all proper-fit
methods must fail. Which specific anchors drive the incompatible placement is
not isolated by this factorial; removing O or changing anchor weights requires a
separate declared experiment.

## Scope, reproduction and next gate

No final PDB, `.dat`, public completion bundle or GUI artifact is published by
this experiment. Installed libraries are read but never edited. Raw diagnostic
coordinate arrays/logs stay under `tmp/pdbfixer-research/causality-run-01/`; the
frozen plan is `causality-plan-01.json`. The independent audit never instantiates
PDBFixer or invokes reconstruction/minimization/dynamics.

```bash
PYTHONPATH=src:. python scripts/research_template_causality.py plan --output <new-plan.json>
PYTHONPATH=src:. python scripts/research_template_causality.py run --plan <new-plan.json> --output <new-run>
PYTHONPATH=src:. python scripts/audit_template_causality.py --plan <new-plan.json> --run <new-run> --output <new-audit.json>
```

The four cases are retrospectively selected mechanisms, not a new family-level
accuracy comparison. These raw pre-minimize failures do not establish unsafe
final production legacy output: downstream chirality defenses remain unchanged.
No successful geometry remedy is claimed from suppressing relaxation or correcting
the rotation determinant.

**664 focused tests passed in 25.80 s**; targeted Ruff passed. Synthetic tests
cover proper-fit rigid geometry/input immutability, a proper-template/hybrid
counterexample, independent determinant/volume/CIP tamper rejection and complete
factorial enrollment without running candidates. This is not a full repository
or cross-platform acceptance run.

Next gate: separately declare an N–CA–C-only proper-fit control to test the
fourth-anchor hypothesis on **all four cases**, not just successful rescues. Keep
the source observations, all raw failures and independent hybrid/stereo/geometry
checks. Do not promote such an untested anchor change or a determinant-only patch
into production; empirical priors, calibrated chemistry and broader replacement
qualification remain separate.

Final saved-stage readback reproduces byte-identically without rerunning controls.
Package mypy, generated CLI/GUI checks, agent maps, prior evidence and all 52
fixture checksums pass. Final read-only checks are retained in
`tmp/pdbfixer-research/validate-causality-final.py`.
