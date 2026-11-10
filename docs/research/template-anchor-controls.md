# Proper-fit anchor controls: omitting O rescues four mechanism cases, not a public backend

Status: **separate retrospective four-case causal ablation**, not native-accuracy
qualification, a complete PDBFixer replacement or a production patch. Continues
the [initialization-stage study](template-initialization-causality.md), where a
properly transformed L template still made a D-side fixed-backbone hybrid.

**Result:** using only N–CA–C to steer the proper fit produces L-Cα at construction
and after normal relaxation in all four enrolled cases. All four pass the same
independent raw/rounded/global/target geometry gate. O remains an unchanged
observed atom; only its role as a fitting anchor is removed. This supports the
anchor-set explanation **on these cases**, not a general safety/accuracy guarantee.

## Frozen intervention and complete factorial

Before new reconstruction, `scripts/research_template_anchors.py plan` froze the
same four cases and roles, without selecting successful rescues:

- 2CI2 I/VAL28 and I/THR55: previously confirmed raw D-Cα whole-sidechain endpoints.
- 1R69 A/THR27: previously stock small-volume/CA-CIP-unassigned case, not confirmed D.
- 2CI2 I/THR22: first whole-sidechain case in original enrollment order, not a
  success-selected or already geometry-valid control.

The full 4×2×2 experiment compares proper fitting with all supplied N/CA/C/O
anchors versus exactly N/CA/C, each with normal relaxation and with both
minimization/dynamics suppressed. The same source bytes, template resources,
proper SVD routine, observed atom order, helper seeds 1–5, force field and numerical
platform remain fixed. No native sidechain coordinates or score-based selector
enter this experiment. No source atom is moved, deleted or restored post hoc.

The adapter requires exactly the four source-bound N/CA/C/O observations at the
target and refuses unexpected sets/order instead of falling back. It does not
handle retained distal atoms, arbitrary partial sidechains, loops, ligands or PTMs.
This restriction is part of the research contract, not a new public capability.

### Supplied anchors are not necessarily used anchors

The unchanged parent tracer continues to record four **supplied** anchors and
labels the rotation class `proper`. A separate `anchor-policy.json` records, for
every construction, the exact **used** atom names/indices, coordinates, centers and
transform. The NCAC subset preserves upstream relative atom order; it is not an
assumed first-three array slice. The new run/plan schemas identify this distinct
experiment; do not reinterpret its old trace fields as four-anchor fitting.

Only a process-local research fit hook implements the new subset; the parent
tracer file and installed scientific libraries are unchanged. All hooks restore.
Every all-anchor arm must reproduce its corresponding saved parent proper arm,
including complete trace JSON, reporter iterates and compressed stage coordinates.
Any mismatch would block interpretation rather than replace the old expectations.

## Execution and independent audit

[Independent used-anchor evidence](template-anchor-controls-evidence.json) verifies:

| Measure | Count |
|---|---:|
| Cases / factorial arms | 4 / 16 |
| Helper-returned seed states | 40 |
| Complete stage snapshots / target reporter iterates | 200 / 43 |
| Independently checked used-anchor fits | 80 |
| Exact parent all-anchor arms / seed states | 8 / 32 |
| Escape-dynamics calls | 0 |
| NCAC normal-relaxation initial / final CA S | 4 / 4 |
| NCAC normal-relaxation common-gate passes | **4/4** |

All eight all-anchor controls reproduce the earlier proper traces/payloads exactly.
Each new NCAC arm returns on seed 1; repeating the same cases with relaxation
suppressed adds controls, not independent accuracy samples. No output is chosen
by a successful seed, native RMSD or a favorable reference.

`scripts/audit_template_anchors.py` does not instantiate PDBFixer or reconstruct,
reflect, minimize or run dynamics. It independently checks complete-stage source
identities/observations, target iterate volumes/backbone, actual supplied-versus-used
transport, proper orthogonality/determinant, used-anchor centroids and the SVD
singular-value bound for the optimal fit residual. From original canonical
resources and saved transforms it reproduces both the full template and actual
hybrid, including O's unchanged observed position. Every returned seed then receives
the **unchanged comparator's** common geometry/CIP audit; no thresholds are relaxed.

The helper's notion of a clean rebuild is not a replacement for that common gate.
For example, the all-anchor THR22 control is L and returns on seed 1 but fails the
common CA–CB length/N–CA–CB angle checks. It was not a geometry-passing control.

## Results with normal relaxation

First-seed construction and final signed Cα volumes, Å³:

| Case | Four-anchor construction → return | N–CA–C construction → return | NCAC returned CIP / common gate |
|---|---:|---:|---|
| 2CI2 VAL28 | −1.555334 → −1.532514 | +2.423129 → +2.423129 | CA S / pass |
| 2CI2 THR55 | −1.275285 → −0.880485 | +2.459139 → +2.459139 | CA S, CB R / pass |
| 1R69 THR27 | −1.083445 → −0.487546 | +2.450350 → +2.306998 | CA S, CB R / pass |
| 2CI2 THR22 | +3.093773 → +3.093773 | +2.365314 → +2.365314 | CA S, CB R / pass |

The all-anchor normal arms pass **0/4** common gates; three have wrong-side CA
placement and the fourth has the length/angle failures above. NCAC passes **4/4**
with O and every other observation still fixed. The no-relaxation NCAC controls
also preserve constructed L geometry and pass their common gates; that is not a
recommendation to skip relaxation on real inputs. THR27's normal run moves its
generated atoms while retaining L geometry, so the result is not merely an
unchanged frozen-coordinate artifact.

Both modes use a proper rotation. The only declared fit intervention is removing
O's orientation/centering influence; the force field still receives observed O.
This isolates a sufficient case-specific anchor-set change, not proof that O's
experimental coordinate is wrong or that all Cα failures have this mechanism.
Optimization-induced inversion on other inputs remains a separate known risk.

## Scope and next gate

No native RMSD, rotamer accuracy, broad-family success or full-legacy comparison
was evaluated here. Four selected mechanism cases cannot qualify a general
anchor-only repair, guarantee future L outputs, or eliminate downstream stereo
checks. No final PDB/`.dat`, public bundle or GUI artifact is published; raw stages
remain private under `tmp/pdbfixer-research/anchor-run-01/`. The frozen protocol is
`anchor-plan-01.json`; the parent source/evidence remains unchanged.

```bash
PYTHONPATH=src:. python scripts/research_template_anchors.py plan --output <new-plan.json>
PYTHONPATH=src:. python scripts/research_template_anchors.py run --plan <new-plan.json> --output <new-run>
PYTHONPATH=src:. python scripts/audit_template_anchors.py --plan <new-plan.json> --run <new-run> --output <new-audit.json>
```

**673 focused tests passed in 25.95 s**, targeted Ruff passed. Synthetic tests
prove NCAC fit independence from observed O, exact full-anchor control equivalence,
source/atom-order binding, immutable inputs and independent rejection of hidden
used-anchor changes, wrong centering or a nonoptimal proper rotation. This is not
a full repository or cross-platform acceptance run.

Next gate: separately freeze evaluation on **all eleven whole-sidechain masks**
from the existing three-system cohort, not only these four rescues. Keep the four
already tested cases distinguishable from the seven not yet executed under this
variant. Assess native accuracy only after proposals and selection, retain all
raw stereo/geometry failures and native-control gate deviations, and do not call
this a new blinded/independent-family cohort. A public or general partial-sidechain
route remains disabled; installed libraries and production defenses stay intact.

Final independent used-anchor/common-gate readback reproduces byte-identically
without new reconstruction. Package mypy, generated CLI/GUI checks, agent maps,
prior evidence and all 52 fixture hashes pass. Final read-only checks:
`tmp/pdbfixer-research/validate-anchor-final.py`.
