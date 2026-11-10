# New mixed-method joint source contract

Status: **private, explicitly versioned joint experiment**. This composes the
default ALA/SER/free-CYS controls, [paired VAL/THR](protein-completion-branched-source.md)
and [fixed-partner singleton](protein-completion-anchored-source.md) methods.
It does not expand old cohort/evidence formats or production callers.

Continuation: [a separate mixed evidence schema](protein-completion-mixed-publication.md)
now preserves admission/raw/rounded refusals and accepted files with complete
deterministic replay. The old cohort/evidence formats and historical hashes below
remain unchanged; this does not establish empirical or production qualification.

## Source-derived modes, not a larger default whitelist

`protein_completion/mixed_cohort.py` introduces `MixedTarget`, `MixedCohortPlan`,
joint completion/materialization/readback and policy
`mixed-canonical-paired-anchored-source-global-thr-beta-v1`.

One to four exact targets are canonicalized into **source encounter order**.
Capability is explicitly frozen to ALA/SER/CYS/VAL/THR, not derived from future
kernel additions. Mode is rederived from actual source/missing atoms:

| Mode | Admitted target | Generator/audit |
|---|---|---|
| `standard` | ALA/SER/free CYS existing masks | Existing staggered controls |
| `paired` | VAL/THR paired branches, with or without CB | Existing paired control |
| `anchored` | One VAL/THR distal atom, CB/partner observed | Existing two-angle singleton solver |

The immutable plan binds the complete source/context/link requests, mode, exact
observed partner where applicable, source THR-Cβ checks and policy. A caller
cannot select a more permissive method by editing mode or duplicating a forged
partner in its observations. Whole gaps, no-op targets, unsupported graph types,
linked target chemistry and other prior source/layout refusals remain excluded.

Only exact selected incomplete centres may be pending **before** reconstruction.
Every target must independently pass its admitted mask. Unselected missing/wrong
THR Cβ or missing/non-L protein Cα still refuses. After output there are **no**
pending exemptions: all source THR Cβ must be complete and R, and all assessed
protein Cα must be L. VAL CB is not a physical stereocentre; ILE/ligand secondary
centres and full source chemistry remain outside this global stereo claim.

## Bounded joint checks and one accepted output

Per-target proposals and all failed individual attempts are retained. Only
source-feasible proposals enter the Cartesian product, capped at **81** jointly
audited combinations. Every target is re-audited using its own method with all
other targets' new atoms included as context. No observed partner is regenerated
or treated as another target's output.

The original cohort now shares its enumeration helper with this new policy;
its plan, result semantics and default validation remain unchanged. Selection is
still capped worst-clearance feasibility with stable tie-breaking, **not** coupled
energetic packing, a rotamer posterior or native-conformation prediction. No
empirical statistics or fallback authority is accepted in this control path.

Selected candidates are rebound and individually/jointly checked before a single
insertion pass. Readback verifies exact original records, combined generated
ownership, serial/metadata constraints, per-coordinate rounding bounds and joint
rounded geometry/clashes, followed by global Cα and THR-Cβ checks. Any failure
returns no accepted PDB or sidecar. Other raw combinations are not retried to
hide a selected combination's rounding failure.

The fresh canonical `DatRecord` owns only combined generated atoms, preserving
case-sensitive chain and insertion-code identity. CB/partners already observed
never become generated. There is no upstream sidecar merge or claim that all
unselected missing sidechains/context-context strain are repaired.

## Old contracts remain isolated

The retained result container and search engine are reused, not the old plan or
evidence semantics. Old cohort entry points reject `MixedCohortPlan`. The old
evidence preflight now rejects wrong plan classes explicitly before staging;
it also rejects a mixed result paired with an otherwise valid all-standard old
plan because deterministic policy/plan replay differs.

Old artifact schemas and valid replay semantics are unchanged. **No mixed evidence
publisher/reader or public API/CLI/GUI/backend was added.** The new service returns
bytes in memory; the research harness explicitly retains its files.

## Executed controls

[Current portable record](protein-completion-mixed-evidence-02.json), raw artifacts
under `tmp/pdbfixer-research/mixed-cohort-02/`. Input masks/outcomes were declared
before that execution; some 1UBQ compositions were already exercised in unit tests.
This is not blind prospective native-model qualification.

| Declared control | Outcome |
|---|---|
| 1UBQ A1–72: THR7 OG1 + VAL26 CG1/CG2 + ALA28 CB + SER57 OG | Accepted; all three modes composed |
| 1UBQ A1–72: THR7 full sidechain + THR55 CG2 | Accepted; multiple source THR centres initially pending |
| Observed 1FR2 chain A: THR11 OG1 + VAL19 paired branches | Accepted |
| 1A3N A1–59: VAL1 full sidechain + THR8 OG1 | Accepted; poorer native-fit VAL1 remains visible |
| Full observed 1A3N chain A, partial occupancy retained | Admission refused |
| Synthetic paired THR + anchored THR | Three source-feasible combinations: one joint pass, two retained clash refusals |
| Synthetic two anchored THRs | Each passes alone; the sole joint combination clashes and refuses |
| Synthetic anchored THR + ALA + boundary water | Raw joint pass; entire output refused after rounding |

Thus five records have accepted PDB/`.dat`, one refuses admission, one refuses
joint geometry and one refuses rounded output. Rejected records have no accepted
structure/sidecar. The accepted records passed **30 RDKit checks** on rounded
source THR centres/selected VAL targets; repeated fixed centres are not independent
samples. Cases use explicit ATOM projections and retained originals, **not intact
complex admission**, occupancy rewriting or force-field preparation.

The 1A3N VAL1 whole-sidechain exact-name RMSD remains **1.864 Å**, not corrected by
joint feasibility. Held-out scoring is reporting-only, after proposal selection
and readback. Other masks/native-fit values and all joint rejections are retained.

The [initial record](protein-completion-mixed-evidence.json) and
`mixed-cohort-01/` remain historical. After adding explicit old-evidence wrong-kind
preflight, `mixed-cohort-02/` repeated the **same deterministic controls**, pinning
current source hashes; these are not new masks or extra scientific seeds.

Verification: **492 focused tests passed in 19.86 s**, sixteen-file package mypy
and targeted Ruff passed. Tests include all **81/81** feasible four-paired controls,
free CYS composition, source/mode/partner/policy forgery, exact case/icode ownership,
new–new clashes, no-joint-solution refusal, whole-output rounding refusal, bounded
contexts/search, candidate/record edits and old cohort/evidence rejection.
No full-suite, empirical or production acceptance is claimed.

```bash
PYTHONPATH=src python -m pytest tests/test_protein_mixed_cohort.py -q
# Use the scientific environment for independent RDKit controls.
PYTHONPATH=src python scripts/probe_mixed_cohort.py --output tmp/pdbfixer-research/new-mixed-controls --summary tmp/pdbfixer-research/new-mixed-summary.json
```

Next: a separate mixed evidence/publication schema with deterministic mode/partner/
global-stereo replay and retained admission/raw/joint/rounded failures. Do not
silently repurpose old schemas. Broader chemistry, harder masks, ILE/multi-chi,
real licensed statistics, comparative science and production adapters remain open.
