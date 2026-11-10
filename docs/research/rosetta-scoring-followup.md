# Predeclared atom-fixed scoring follow-up: crambin, lysozyme and barstar

Status: **predeclared cohort and frozen scoring-rule follow-up with an explicit
post-reconstruction serialization amendment**. Not a pristine unamended protocol,
randomized/blinded benchmark or broad qualification. Follows the
[retrospective scoring pilot](rosetta-fixed-scoring.md).

Continuation: [exact two-target controls](rosetta-joint-controls.md) now enumerate
all restricted joint combinations, with canonical sidecars and both target graph
views verified. Their mostly degenerate, numerically additive score tables do
**not** demonstrate coupled-packing benefit; keep the scientific scope separate.

## Enrollment before new proposals or scores

`scripts/prospective_branched_cohort.py` enrolled **1CRN (crambin), 1AKI (hen
egg-white lysozyme) and 1BTA (barstar)** before downloading/inspecting their inputs
or running reconstruction. These named systems were not previously used in this
reconstruction/scoring study; this is deliberately chosen diversity, not random
sampling or a formal family-independence analysis.

The rule remained: first two observed VAL and first two observed THR identities
in chain-A file order, four masks each (distal pair, whole sidechain, either distal
singleton). All **48 masks** were frozen without consulting source admission,
own success, template success, reference accuracy or Rosetta energy.

Projection selects the first coordinate model if multiple models exist, retains
every selected chain-A ATOM byte, including altloc/occupancy/H, and preserves
relevant SSBOND/LINK/CONECT annotations. Excluded link endpoints remain present in
metadata so admission can refuse rather than silently sever chemistry. This is
protein-chain context, not intact complexes/water/heterogen qualification.

Original RCSB downloads are tracked byte-for-byte under
`tests/fixtures/protein_context/scoring_followup/`; provenance and SHA-256 are in
the [fixture catalog](../../tests/fixtures/README.md) and
[manifest](../../tests/fixtures/MANIFEST.sha256). Neither raw structures nor input
hydrogens were repaired to improve eligibility.

The score function, all eligibility rules, H protocol, physical/virtual inventory,
terminal OXT allowance, topology scope, seed and relevant implementation hashes
are bound to the earlier selected policy **before these new outcomes**. No weights,
angles, clash cutoffs, candidate grids or selection tie-breaks were retuned.

## Explicit CONECT-order amendment

The initial projector incorrectly placed retained CONECT records before coordinates.
OpenMM failed at `_current_model.connects` before a model existed: **32 template
execution/protocol blockers**, not 32 failed scientific reconstructions. This was
our serialization error, not a newly established PDBFixer reconstruction defect.
Initial own outcomes and all logs remain retained.

A separate amended enrollment moved CONECT **after** coordinates. The audit proves:

- All 48 case identities/masks remain identical.
- Across all 54 original/projected/masked artifacts, **every non-CONECT byte** and
  **every CONECT record** remains identical; only those records' position changes.
- All 32 own results reproduce exact additions, audits and selected indices;
  all 16 source refusals also reproduce. Source/request digests necessarily change.
- The original scoring rule/code remains unchanged; new-cohort scoring had not begun.

This amendment happened **after initial reconstruction outcomes were visible**.
Preserve that limitation when discussing prospective evidence. The original
enrollment, broken-input run, matching old projector bytes, amendment and new run
are separate artifacts, not overwritten expected hashes.

## Complete outcomes

[Follow-up evidence](rosetta-scoring-followup-evidence.json),
[independent scoring readback](rosetta-scoring-followup-evidence-readback.json),
[reconstruction record](rosetta-scoring-followup-reconstruction.json) and
[native-reference audit](rosetta-scoring-followup-native-audit.json):

| System | Declared masks | Admitted | Own common passes | Template common passes |
|---|---:|---:|---:|---:|
| 1CRN | 16 | 16 | 15 | 9 |
| 1AKI | 16 | 16 | 16 | 5 |
| 1BTA | 16 | 0 | Not executed | Not executed |
| Total | **48** | **32** | **31** | **14** |

All 16 barstar refusals retain the source-scope reason: canonical heavy-atom
names/elements are required. Its **714 deposited H coordinates** remain in the source and
projection. They were not stripped to turn a refusal into a pass. No replacement
protein, chain/model substitution or outcome-driven re-enrollment was used.

After the serialization correction, template execution has **zero blockers**.
On the 32 admitted cases, both own/template pass on 13, own-only on 18,
template-only on one, neither on zero. These are unchanged prototype common gates,
not calibrated physical-quality or production safety scores.

All **32/32 native controls pass** that common gate. This does not establish broad
gate calibration or erase the previous 1TIM native 1–4 limitation. Own singleton
**1CRN VAL8 CG1** still refuses its CG1–C clash test; the template and native control
pass. Scoring does not rescue it or modify the old cutoff.

## Frozen scoring rule on the new cohort

All 64 own attempts, 32 template proposals and 32 native controls were scored:
**128/128 successful physical/topology snapshot audits**, with exact physical
heavy-coordinate preservation. No repacking, minimization, reflection or post-hoc
observation restoration. H is rebuilt separately for each fully specified rounded
heavy input, without H optimization. Extra OXT, virtual and non-target H coordinates
remain identical within each case; no comparison blockers were found.

Selection is still minimum total `ref2015` score within the **complete unchanged
own source-and-common-gate-feasible pool**, with stable index ties. Missing feasible
scores would block selection. Native controls are imported only after own selection
is committed, and reference errors are reporting only.

Of 32 admitted cases: 23 have one feasible own proposal, one has two, seven have
three, and one has none. Only **eight cases** offer a genuine ranking choice.
There are 31 selected cases under both methods; only **one** changes:

| Case | Old rounded exact-name RMSD | Score-selected RMSD |
|---|---:|---:|
| 1AKI VAL2, CG1/CG2 missing, fixed CB | 2.414 Å | **0.147 Å** |

Allowing only the existing **both-generated VAL methyl swap** reduces the old
error to **1.771 Å**, while the new error stays **0.147 Å**. Thus the improvement
is not explained solely by CG1/CG2 relabelling. This is a post-selection reporting
sensitivity check, not a new ranking criterion or permission to swap fixed atoms.

Matched all-own-feasible, PDB-rounded generated-coordinate medians:

| Class | Same-case denominator | Old selector | Score selector |
|---|---:|---:|---:|
| Paired / whole sidechain | 16 | 0.247 Å | 0.242 Å |
| Anchored singleton | 15 | 0.115 Å | 0.115 Å |

The other seven actual ranking choices remain unchanged. This is one additional
useful selection on one lysozyme residue, not a general accuracy improvement.
Do not merge these denominators with template-mutual-pass or earlier-pilot metrics.
For the original clearance selector's template-mutual-pass raw subsets, template paired accuracy remains
better (n=6 medians own/template **0.248/0.187 Å**); singleton n=7 is
**0.063/0.172 Å**. Scoring has not increased coverage or broadened chemistry.

Descriptive topology retains all **40 old rounded clash instances** across 128
inputs: ten known 1–4 and 30 with no short path in the partial graph. The latter
remain unresolved, not proven nonbonded. No contact class is newly exempted.

## Artifacts and next step

All retained paths are under `tmp/pdbfixer-research/`:

- `prospective-enrollment-{01,02}.json`, `prospective-inputs-01/` and
  `prospective-masks-{01,02}/`: original/amended enrollment and exact inputs.
- `prospective-reconstruction-{01,02}/`: initial protocol failures and corrected
  reconstruction, with all attempts and seeded template states.
- `prospective-conect-amendment-audit.json` and
  `prospective-own-replay-after-conect-fix.json`: serialization and own-result proofs.
- `prospective-score-plan-01/` / `prospective-score-run-01/`: final scoring plan,
  archived legacy builder plan, pre-native selections and full scored snapshots.

The final scoring plan explicitly carries prospective enrollment and amendment
provenance; the archived base builder keeps its old retrospective label. The
generic independent readback's benchmark-qualification flag remains false. These
are explicit provenance distinctions, not rewritten historical results.

New research scripts are `prospective_branched_cohort.py` and
`prospective_scoring_followup.py`; the reconstruction/scoring kernels and old
independent readers remain byte-unchanged. Tests cover full unfiltered enrollment,
model/chain selection, H retention, metadata endpoints/order, real OpenMM reads and
disulfides, and rejection of policy retuning.

Verification: **591 focused tests passed in 24.61 s**, targeted Ruff passed. This
is not the full repository suite or broad scientific qualification.
Package mypy, generated CLI/GUI checks, agent maps and all 49 fixture checksums
also pass. Saved scoring and native-audit readbacks replay without new proposals;
older evidence and unchanged scientific policies remain valid.

Next: exact **joint fixed-atom** controls on separately declared small cohorts,
before introducing an interaction-graph approximation. More independent systems,
candidate diversity, calibrated geometry/protonation and production integration
remain separate unfinished gates. No upstream code/table redistribution or new
PyRosetta package dependency is introduced.
