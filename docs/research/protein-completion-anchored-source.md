# Explicit singleton source policy: fixed partner, global stereo and PDB readback

Status: **separate private source experiment**, following the
[typed observed-partner kernel](protein-completion-anchored.md). The paired VAL/THR
policy, default ALA/SER/CYS source/cohort/evidence APIs and production callers keep
their previous capability and semantics. This is not a new public backend.

Continuation: a [distinct mixed-method joint source contract](protein-completion-mixed.md)
now composes standard, paired and anchored masks with cross-target readback. It
does not change this single-target policy or enroll old cohort/evidence formats.

## New explicit entry points

`protein_completion/anchored_source.py` provides:

- `plan_anchored_source(source, exact_target)`: rederive the complete request and
  links from source bytes, requiring exactly one missing VAL/THR distal atom with
  observed CB and the other distal partner.
- `AnchoredSourcePlan`: the full request, **exact fixed-partner observation**,
  source THR-Cβ checks and distinct policy
  `singleton-val-thr-source-fixed-partner-global-thr-beta-v1`, all digest-bound.
- `materialize_anchored_source(...)`: rebind source/partner/context/policy, audit
  the raw one-atom candidate, insert it in memory and independently audit readback.
- `audit_anchored_materialized_source(...)`: exact original records/identities,
  generated ownership, bounded coordinate rounding, anchored geometry/stereo/chi/
  clashes, global Cα and post-output THR-Cβ checks.
- `MaterializedAnchoredSource.dat_record()`: a fresh canonical `DatRecord` with
  exactly the generated singleton. The observed partner and CB never enter it.

Source assembly/layout/link/graph/insertion/readback helpers are reused. The
internal admission helper accepts an explicitly selected validator, but all old
entry points keep their old default; no serialized input or ordinary caller can
implicitly opt an older policy into this new mask contract. Paired plans cannot
be passed as anchored plans or vice versa.

Changing both a claimed partner and its duplicate observation is still detected:
readback rederives them from the original source, not just from each other. Case,
insertion code, coordinates, element, occupancy, context and complete source SHA
remain part of binding. No observed atom is moved, renamed or reflected.

## Scientific and filesystem boundaries

The new policy shares the explicit canonical-heavy protein/water input scope,
unit occupancy/no-altloc rules, link/symmetry refusals and global protein Cα gate.
Every unselected THR Cβ must already be complete and R; only the exact target may
be pending. **No pending centre is allowed after readback.** Negative/nearplanar
or undefined source THR centres refuse even when Cα remains correctly L.

VAL CB remains a naming convention, not a physical stereocentre. ILE/ligand/other
secondary centres are not globally qualified. Existing hydrogens, unsupported
chemistry and linked target sidechains do not gain repair support.

Raw geometry is independently audited before writing. Rounded coordinates must
differ by <=0.000500001 Å per component, keep existing physical thresholds and
pass the separate 0.15° rounded measured-chi tolerance. Boundary-water tests show
both an O-tip and C-tip raw pass can become a rounded clash; the service raises
without returning accepted bytes. It never switches hemisphere, moves the partner
or retries to disguise that failure.

Original headers, partner ANISOU, valid CONECT and bare/sparse TER remain exact.
Serial exhaustion returns no partial output. The core service writes **no files**;
only the research harness persists explicitly requested results. There is no new
evidence-bundle kind, CLI/GUI registration, cohort enrollment or production route.

## Executed readback controls

[Portable evidence](protein-completion-anchored-source-evidence.json), complete
retained artifacts under `tmp/pdbfixer-research/anchored-source-01/`.

All **18 cases were declared before execution**:

| Declared projection | Singleton masks | Outcome |
|---|---:|---|
| 1UBQ ATOM A1–72: VAL5/26, THR7/55, either distal atom withheld | 8 | Accepted after independent rounded/source/sidecar checks |
| All observed 1FR2 chain-A ATOM records: VAL19, THR11, either atom | 4 | Accepted |
| 1A3N ATOM A1–59: VAL1, THR8, either atom | 4 | Accepted |
| Full observed 1A3N chain-A ATOM records: VAL1/THR8 singleton controls | 2 | Admission refused because partial occupancy remains |

These are explicit projections with original ATOM bytes, withheld masks and bare
TER/END, **not whole deposited-complex acceptance**. Original structures and hashes
are retained. No occupancy was normalized or partial atom silently selected.

All 16 accepted files passed generated-only `.dat` readback and **108 independent
RDKit checks** on rounded THR centres and the selected VAL/THR target. Those calls
repeat fixed centres across masks; they are not 108 independent scientific samples.
The fixed distal partner remains physically byte-identical in every result.

Held-out singleton exact-name error was computed only after selection/readback:
**0.012–0.243 Å**. This easier fixed-CB/fixed-partner task does not fix the earlier
whole-sidechain 1A3N VAL1 1.864 Å result. There is no empirical rotamer data,
native-quality guarantee, template comparator or MD evaluation in this step.

The existing source-control harness gained an explicit research-only
`--mode anchored`; default `--mode paired` is unchanged. Modes select different
planners, generators, materializers, policy IDs and result schemas; they do not
alter application command registries or historical records.

Verification: **473 focused tests passed in 17.65 s**; fifteen-file package mypy
and targeted Ruff passed. Tests include coherent partner/observation forgery,
policy/context substitution, old API isolation, global source D-Cα, wrong/undefined
THR centres, case/icode exemptions, source/partner edits, coordinate tampering,
rounding-created clashes, ANISOU/CONECT and sparse/exhausted serials. Independent
CIP tests cover rounded source results with either branch absent.

```bash
# Use the scientific environment for independent RDKit controls.
PYTHONPATH=src python -m pytest tests/test_protein_anchored_source.py tests/test_protein_anchored_stereo.py -q
PYTHONPATH=src python scripts/probe_branched_source.py --mode anchored --output tmp/pdbfixer-research/new-anchored-source --summary tmp/pdbfixer-research/new-anchored-source-summary.json
```

Next gates: versioned joint paired/singleton composition and replay-verified
evidence integration without changing old formats, broader chemistry and harder
mask comparisons. Genuine empirical statistics, ILE/multi-chi support and all
production PDBFixer consumer adapters remain open.
