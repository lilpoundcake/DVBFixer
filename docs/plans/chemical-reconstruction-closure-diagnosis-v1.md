# WP2Q: read-only diagnosis of the two retained CD-only rejections

Date: 2026-11-22; baseline `6810f76`; persistent
`feature/chemical-reconstruction`. Freeze BEFORE new diagnostic measurements.

## Scope and documentation plan

User approved diagnosing 1TIM-CD and 2CI2-CD after WP2P. This is retrospective
analysis of saved coordinates and construction declarations, NOT a new proposal,
repair, benchmark rerun, solver, H completion or permission to enable a backend.
Keep the complete nine-cell WP2P context, including the passing 1CRN-CD and the
six larger-mask cases. Do not turn replay/diagnostic counts into new trials.

Documentation: this protocol, frozen saved-artifact ledger, diagnostic report,
validation JSON, active handoff/task discovery and local quality note. Pending
links: none until report/validation are created. Version stays 0.10.2. No source,
fixture, chemical state, mask, metric bound, constructor/compiler, parent admission,
dependency/default/public/DatRecord or production/data-license change.

## Frozen inputs and independent source qualification

The adjacent `chemical-reconstruction-closure-diagnosis-inputs-v1.json` pins ALL
nine saved WP2P evidence files plus summary `8501897ac1b3ccb841f17e7ef8bca387c69206bfdebd7b94385beca2e4d60a59`.
Reuse the existing WP2P constructor-free source/authority/parent/template/H/native
readback before any diagnostic interpretation. Do not silently re-create missing
artifacts, run a new baseline, replace coordinates, alter the old report, choose
another site/mask, or update its pass/failure denominator. Missing/corrupt files
are evidence blockers before new output. Retain reports inside the project.

## Measurements, not new coordinate construction

1. Enumerate each actual target heavy bond/angle and explicit constraint, using
   existing authority_geometry, explicit metric policy and measurement owners.
   Report exact identities, authority/allowed values, measured values, pass status,
   raw-vs-rounded phase and whether the term is fixed-only or generated-involving.
   Measure ORIGINAL own-request semantics separately from COMMON v2-port semantics.
   Every explained failure must map to the existing independent candidate audit;
   no new diagnosis may weaken or replace that audit.
2. Inspect the saved deterministic chain plan and recompile its declaration only
   (no propose/extend/rotation/solver). For each placement report the actual anchor
   order, intended length/angle/ideal torsion and measured saved output values.
   Enumerate generated-involving graph edges which are NOT enforced as placement
   length edges. Identify ring membership by generic graph connectivity after edge
   removal, not a PRO atom-name rule. A missing edge in the construction frame is
   not a missing authoritative bond: the independent validator must still enforce it.
3. Compare saved raw/rounded original coordinates with retained PDBFixer and native
   reference under identical scalar terms, and saved projections where available.
   References were already extracted after WP2P proposal phases; here they are
   diagnostic ONLY, never a source of a new candidate or an alternative anchor/chi.
   Do not construct counterfactual poses or use native to select a "fix".
4. For a generated atom with two observed bonded neighbours, report a conservative
   necessary distance-shell compatibility test and a conservative cosine-interval
   test for the angle at the missing atom, under the existing bond/angle tolerances.
   This is scalar interval arithmetic on fixed distances and authority values,
   NOT sphere-intersection point generation, minimization, or complete feasibility.
   Passing these necessary tests does NOT prove a jointly feasible pose with all
   other angular/stereo/plane/context constraints. A failed chosen pose alone does
   NOT prove that fixed observations are inconsistent or must be moved.

Preserve actual source/authority incompatibilities, fixed-only failures, unmeasured
or degenerate terms and every old native-score worsening. Diagnose the chain-only
assumption, not a claimed PDBFixer upstream cause. No proposed general ring-closure
algorithm is scientifically qualified by this analysis.

## Guards, tests and outputs

Disable all coordinate constructors (including graph internal/template, extend/
two_angle_tip, PDBFixer/shared rebuild, projection optimizer, H builders, search and
rotators) during diagnostic execution and independent saved report readback.
Metadata-only Topology for the existing template global-L audit remains permitted.
No new Candidate/ProjectionResult/ConformerSearchResult or fake strategy records.

Tests use synthetic/tracked authority/request fixtures, not required ignored cohort
outputs. Cover generic cycle/noncycle edge classification, exact selected-frame
transport, both phases, fixed/generated roles, distance/cosine interval degeneracy
and sound conservative rejection, incomplete identities, duplicate/nonfinite data,
unknown constraints, constructor-free execution and input/report tampering. Use
the existing scalar measurements rather than a second chemical/variant authority.

Write a separately tagged diagnostic JSON with source/request/plan/candidate hashes,
per-term tables and explicit limitations; reproduce it without constructors. Keep
WP2P's 9/7/2 accounting and its 525 replay-file evidence unchanged. Run focused,
boundary/full-fast, Ruff, scoped/source-bound typing, maps/generated/version/fixture
checks. Verify historical source/artifact hashes; do NOT rerun scientific controls
merely for diagnostics. Preserve ambient 814 errors, GUI dependency state, WP2F
timeout, old refusals and native worsenings as separate open gates.

## Decision boundary

Deliver the measured cause(s), what the constructor actually enforced, and what
the necessary compatibility tests can/cannot establish. Any changed closure
constructor, frame selection, rejected-baseline projection admission or broader
cohort requires a separate prospective protocol and independent validation. Do
not retrofit a fix or pass into the frozen WP2P records during this slice.
