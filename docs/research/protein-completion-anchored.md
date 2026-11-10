# Explicit observed-partner singleton VAL/THR kernel

Status: **separate private typed-kernel experiment**. This adds a narrower mask
than the [paired-branch kernel](protein-completion-branched.md), not a new mode
inside existing source/cohort/evidence APIs. The
[paired VAL/THR source adapter](protein-completion-branched-source.md) still
explicitly refuses these singleton masks.

Continuation: [an explicit anchored-source policy](protein-completion-anchored-source.md)
now binds the fixed partner, global THR-Cβ and rounded PDB/.dat readback. It is a
separate entry point: earlier typed-only measurements/hashes and paired-source
refusals remain historical and unchanged.

## New construction, not chi sampling

`protein_completion/anchored_branch.py` provides `complete_anchored_branch` and
`audit_anchored_branch`. The only admitted masks are one of VAL CG1/CG2 or THR
OG1/CG2, with CB and the other distal branch **already observed and immutable**.
Missing CB, both missing branches, missing backbone, linked target chemistry and
the existing occupancy/altloc/context defects do not gain an implicit fallback.

`geometry.py::two_angle_tip` solves the intersection of two angular cones about
CB, using the fixed CA and partner directions. It enforces the declared new-atom
bond length and two ideal 109.4712° angles. The algebraic hemisphere is selected
by natural THR CB=R or the existing VAL named-branch convention. When the missing
atom is the first rather than second branch, the volume sign is reversed in the
construction's anchor ordering; the independent named-volume check is unchanged.

No retained atom is moved, reflected, optimized or returned as generated. There
is exactly **one** proposal, not three chi samples, a rotamer prior or a global
energy optimization. The chi in candidate metadata is measured from the resulting
geometry (possibly determined by an already observed primary branch), not forced
onto a staggered grid. The independent audit recomputes it and rejects nonfinite,
out-of-range or false metadata, wrong ownership/stereo and generated clashes.

Collinear, incompatible or nearplanar cone intersections refuse. A solvable cone
can still produce a clashing proposal: that failed attempt is retained, with no
optimization, resampling or threshold relaxation. Numerical construction tests
explicitly retain such refusals instead of assuming every synthetic angle is a
physically acceptable conformation.

Result revision: `canonical-val-thr-observed-branch-anchor-v1`; policy:
`observed-partner-two-angle-closure-v1; no empirical rotamer probabilities`.
There is no empirical-library argument. The default `complete_sidechain`,
`audit_candidate` and all previous materializers keep their original singleton
refusal. Shared validation internals select this explicitly isolated mask/chi
contract; old callers do not opt in merely because the implementation exists.

## Independent stereo and immutable-anchor controls

- Synthetic non-grid chi examples recover the known algebraic tip to <1e-12 Å
  without ever passing that withheld tip to construction. Other source geometry
  is unchanged; synthetic clashes remain explicit rejections.
- Reflecting only the proposed tip across its CA–CB–fixed-partner plane preserves
  both constraint angles. Even with corrected measured chi metadata, independent
  stereo rejects the other solution. For THR, RDKit confirms CB R→S while CA
  remains S; VAL still has no physical CB stereocentre.
- Real 1UBQ controls exercise both possible missing branches of VAL5 and THR7
  using the independent RDKit oracle. Core code imports no RDKit.
- Tests cover immutable anchors, single-attempt retention, impossible/collinear
  anchors, invalid hemisphere selection, forged masks/digests/chi and explicit
  refusal by the existing paired source API.

## Executed 16-mask control

[Portable evidence](protein-completion-anchored-evidence.json); originals, declared
plan, masked inputs, full typed requests, attempts and independent audits under
`tmp/pdbfixer-research/anchored-branch-01/`. Python 3.12.13, RDKit 2025.09.6.

All masks were declared before execution: both singleton possibilities for VAL5,
VAL26, THR7 and THR55 on the explicit 1UBQ A1–72 projection; VAL19/THR11 on all
observed 1FR2 chain-A ATOM records; and VAL1/THR8 on the 1A3N A1–59 projection.
The research harness reuses the existing PDB parser/layout and observed-peptide
continuity helper, supplying the full declared projection minus one atom. It is
**not** a new public/source-bound request admission or global stereo policy.

All **16/16** one-atom proposals passed the scoped geometry/clash audit and
independent CIP checks. Withheld exact-name coordinate error was consulted only
after generation/selection:

| Residue | First branch missing | CG2 missing |
|---|---:|---:|
| 1UBQ VAL5 | 0.026 Å | 0.032 Å |
| 1UBQ VAL26 | 0.090 Å | 0.138 Å |
| 1UBQ THR7 | 0.231 Å | 0.118 Å |
| 1UBQ THR55 | 0.243 Å | 0.215 Å |
| 1FR2 VAL19 | 0.068 Å | 0.038 Å |
| 1FR2 THR11 | 0.011 Å | 0.024 Å |
| 1A3N VAL1 | 0.099 Å | 0.069 Å |
| 1A3N THR8 | 0.081 Å | 0.104 Å |

This is a **one-atom, observed-partner task**, not evidence that the earlier
whole-sidechain 1A3N VAL1 result (1.864 Å) is fixed. It uses more observed
information and cannot be compared as an improvement on that harder mask.
There was no template comparator, empirical dataset, truth-based reranking or
native-quality qualification. Requests remain immutable; only generated tuples
are returned. No PDB or `.dat` was materialized for these new singleton cases.

Verification: **443 focused tests passed in 15.40 s**; fourteen-file package mypy
and targeted Ruff passed. Historical paired/default/evidence records are not
rewritten to match changes in shared validator source hashes.

```bash
# Independent CIP tests and the probe require the scientific environment.
PYTHONPATH=src python -m pytest tests/test_protein_anchored_branch.py tests/test_protein_anchored_stereo.py -q
PYTHONPATH=src python scripts/probe_anchored_branch.py --output tmp/pdbfixer-research/new-anchored-controls --summary tmp/pdbfixer-research/new-anchored-summary.json
```

Next: a separate explicit anchored-source policy must bind the observed partner,
global source stereo, rounding, immutable physical records and generated-only
ownership before any materialization/publication enrollment. Multi-chi/ILE,
broader source chemistry, genuine empirical data and production adapters remain
open. Never silently broaden the paired-source v1 policy or legacy consumers.
