# Graph-methyl diversity study: source admission and preparation

Status: **source/preparation checkpoint only**. No new graph-H control outcomes
have been evaluated. Follows the limited one-protein
[whole-system controls](methyl-whole-control.md), including their adverse findings.
Production preparation, safeguards and public completion capabilities are unchanged.

## Fixed convenience enrollment

The [immutable enrollment](methyl-diversity-enrollment.json) was committed before
source admission or graph-initializer implementation (`d43d4f5`), SHA256
`716060b015b1e553ef127318e573f1b60ff7cf601d3e9fd117502e16a9a896b5`.
These eight repository fixtures have appeared in other research: they are not
blinded proteins or a temporal holdout. No replacements are permitted.

Projection uses the existing first-coordinate-model/first-encountered-ATOM-chain
policy. Exact ATOM records, including H, altloc, occupancy and numbering, and
touching LINK/SSBOND/CONECT endpoints are retained. This is an explicit protein
view, not full-assembly/heterogen acceptance or an occupancy cleanup.

## Independently audited denominator

[Complete admission/preparation evidence](methyl-diversity-admission-evidence.json)
binds original/projection bytes, the frozen plan, runtime resources, canonical
`.dat` decoding, execution controls and all retained artifact hashes.

| Input | Source admission / preparation | Prepared added heavy / all atoms |
|---|---|---:|
| 1CRN | admitted / prepared | 0 / 315 |
| 1CTF | admitted / prepared | 0 / 518 |
| 1TIM | admitted / prepared | 1 / 1888 |
| 1UBQ | refused source | — |
| 1EMV | admitted / prepared | 1 / 618 |
| 1FR2 | admitted / prepared | 1 / 617 |
| 2VLN | refused source | — |
| 2VLQ | refused source | — |

All three refusals reproduce the shared layout reason:
`altloc/partial occupancy or noncanonical residue number`. That combined reason
is not an assertion that every listed defect occurs in each input. Nothing was
stripped, cropped or replaced to gain entry; refused views never ran preparation.

All five admitted views ran exactly one unchanged legacy preparation with
standard naming and default PROPKA/Reduce/pH. Captured controls confirm Reference
PDBFixer and Reference H placement with seed 1, restored random state and restored
hooks. Borrowed projected input bytes remain exact. There were no execution
blockers, timeouts or budget nonexecutions. Runtime: PDBFixer 1.12.0, OpenMM 8.5.2,
Reference, 31 bound soft.xml/template resources. Scope checks are canonical source
layout/global Cα/THR-Cβ and declared-link endpoints, **not exhaustive chemistry or
prepared-geometry qualification**. A successful prepare is not an H-control pass.

Private records live under `tmp/pdbfixer-research/methyl-diversity-*`: original
source plan `methyl-diversity-source-plan-01.json`, projections, logs and
`methyl-diversity-prepared-01/`. These runtime paths do not transfer through Git.

## Frozen scientific boundary and next work

Graph methyls must be carbon with exactly three newly mobile H and one heavy
bonded neighbor in canonical protein context. Observed/kept H are never rewritten.
The heavy second-shell frame prefers the same residue, then full case-sensitive
identity; undefined frames are refused. Canonical tetrahedral directions retain
each original C–H radius and every nonmethyl initial coordinate.

Before any arm evaluation, implement/test and freeze the generic algorithm and
capture original H Systems/starts/endpoints for exactly seeds 1/2, plus their first
full-System starts (stop before whole minimization). Preserve actual topology
index/object/bond lineage and both label inventories through variant restoration;
do not guess name aliases. Three arms are unchanged, phase-zero graph tetrahedral
and initial-energy-selected graph tetrahedral. The last starts with every group at
zero, then makes **one** lexicographic parent-identity pass over 0–345° in 15° steps,
no revisits, fixed 1e−8 kJ/mol/shortest-signed-phase ties, no final/native scoring.
Exact unchanged endpoints and all adverse geometry are mandatory.

No production repair or general methyl qualification follows from these admissions.
Original unforced 1CRN THR2 causality remains unresolved; empirical-library
production provenance/license/approval gates remain unchanged.

Code: `scripts/freeze_methyl_diversity.py`, `scripts/prepare_methyl_diversity.py`,
`scripts/audit_methyl_diversity_admission.py`. Focused enrollment/admission/
preparation/audit tests: **19 passed**. Independent audit publications agree
byte-for-byte; targeted Ruff passes.
