# P2 publication, exact reuse and anchor label sensitivity

Date: 2026-11-10. Continuation of the
[private P2 and calibration controls](gap-portfolios-and-anchor-calibration.md).
P3 is explicitly excluded. Public CLI/Homology/GUI portfolio admission and higher
production anchor ceilings remain disabled.

## Post-assembly bounded-anchor continuation

The original generated-only P2 workflow below is retained unchanged. A new
explicit private escalation, `post-portfolio-bounded-anchors-v1`, now connects
the already implemented bounded-anchor kernel and independent staged acceptance
to a **completed patchwork assembly**. This is not selection of moved-anchor
donors: every donor still shares the unchanged observed source frame.

The order is assembly → generated-only issue-local attempt → residual-driven
anchor authorization → bounded-anchor attempt → independent whole-complex audit.
As in P0/P1, the anchor attempt restarts from the assembled projection, **not**
from a failed generated-only result. The latter is preserved as eligibility
evidence. The source and frozen sampling request are never rewritten.

Only direct flanks named by residual actionable findings may move; preexisting
source-only findings cannot authorize motion. Generated mobility covers gaps
attributed to assembly or residual defects. Clean generated gaps, other observed
atoms, all metadata and noncoordinate records remain byte-identical. Unattributed
geometry findings retain the existing conservative generated-gap scope; they
cannot independently authorize an observed anchor. Geometry is scored after PDB
rounding, and both kernel walls and independent acceptance use the unchanged
per-anchor limits in the original source frame:

- backbone RMSD ≤ **0.5 Å**;
- Cα displacement ≤ **0.5 Å**;
- maximum heavy-atom displacement ≤ **1.25 Å**.

`portfolio_anchors.py` reuses `refine_generated_region`, experimental-B-aware
`build_anchor_restraints`, and `audit_staged_anchor_repair`. The portfolio audit
additionally rechecks exact clean-gap/nonmobile-record custody and global all-L.
Neither a local pass, a restraint strength nor a runner receipt grants movement
permission or replaces final geometry acceptance.

### Private workflow

For a new generated-only local attempt with optional anchor escalation:

```bash
python scripts/repair_gap_portfolio.py --portfolio <selection> \
  --output <new-repair-directory> --relax-anchors \
  --cohort <original-retained-cohort> --case 9CQD \
  --bundle <new-anchored-bundle>
```

For an already retained failed local attempt, do not rerun it:

```bash
python scripts/publish_gap_portfolio.py --portfolio <selection> \
  --repair <retained-generated-only-attempt> --cohort <original-retained-cohort> \
  --case 9CQD --work <new-evidence-workspace> --evidence-only
python scripts/repair_portfolio_anchors.py --workspace <new-evidence-workspace> \
  --output <new-anchor-attempt> --bundle <new-anchored-bundle> \
  --refinement-seed 7 --max-iterations 250 --restart-count 4 --platform CPU
```

`--evidence-only` creates complete donor/preparation custody, not an accepted
bundle. The repair independently audits that custody before calling the kernel;
only the expected geometry failure is permitted as input. Existing directories
are never adopted. Raw assembly, failed generated-only result, trial coordinates,
execution failures and final reports remain retained. An already clean local
attempt is not escalated. Unsupported or unjustified anchor motion fails closed.

### Versioning, publication and reuse

Anchored evidence uses **schema 2**, with explicit stage-1 bytes, separate anchor
operation/report and adapter/kernel/conditioning/acceptance/policy snapshots.
Schema-1 evidence and old bundles retain their original semantics and exact
resume identity. The new reader explicitly dispatches evidence versions; an old
schema-1-only reader cannot silently adopt anchor motion. The private bundle
index remains the disjoint `experimental-gap-portfolio-bundle` kind.

Publication re-audits the original donors, selected projection, generated-only
local attempt, authorized generated/observed masks, movement limits and final
geometry. Rehashed masks, increased limits or forged audit receipts fail. The
`.dat` still lists **only diffusion-generated atoms**; moved observed anchors are
recorded separately in validation/provenance and are never relabelled generated.
Resume remains read-only and binds every operation/source snapshot. Repeated
anchor escalation from an anchored result is refused; bounds cannot accumulate.

This is private integration, not public CLI/Homology/GUI admission or broader
scientific qualification. The historical 8/10 P2 outcome and original 27/50 seed
outcomes below are not replaced by these additional controls. P3 and larger
movement limits are not involved.

Verification of this continuation: **2,289 fast tests passed**, three existing
optional Torch/Biotite skips, 33 slow cases deselected; all 20 new anchor-portfolio
tests are included. All eight previously published schema-1 bundles passed
independent readback and exact resume under the extended reader with every file
byte unchanged. These are software/custody checks, not new scientific samples or
evidence of general scientific qualification. The separate CPU controls completed
with the outcomes below. Evidence: `tmp/portfolio-anchors-v1/fast-tests.json` and
`tmp/portfolio-anchors-v1/old-bundle-compatibility-01/summary.json`.

### Executed saved-assembly controls

Both cases used optimizer seed 7, CPU, 250 iterations and four restarts, with
no new neural sampling, thresholds or donor substitution. These are additional
routes on retained candidates, not replacements for historical seed/P2 outcomes.

| Case | Mobile generated gaps | Authorized anchors | Final private result |
|---|---|---|---|
| 9CQD | W26–28 (3 residues) | W25, W29 | Geometry/movement/custody passed; schema-2 bundle published, independently read and exactly resumed |
| 1HE8 | A322–352 (31), A896–901 (6), A968–980 (13) | A895, A981 | Rejected: A976 Ramachandran and A895 movement above limits; diagnostics retained, no accepted bundle |

The 1HE8 generated mask includes A322–352 because the anchor stage restarts from
the assembled projection, where that gap was still defective, not the later
generated-only coordinates. It does not silently discard that repair task.

| Observed anchor | Backbone RMSD (Å) | Cα displacement (Å) | Maximum heavy displacement (Å) | Within unchanged limits |
|---|---:|---:|---:|---|
| 9CQD W25 | 0.1995 | 0.1027 | 0.3215 | yes |
| 9CQD W29 | 0.4850 | 0.3373 | 1.0484 | yes |
| 1HE8 A895 | 0.5645 | 0.5091 | 1.2844 | **no** |
| 1HE8 A981 | 0.3372 | 0.3569 | 0.5572 | yes |

1HE8 retains a gross A976 Ramachandran outlier (φ 162.6°, ψ −23.2°); passing
clash/junction checks alone cannot grant acceptance. Its A895 overrun is retained
as a genuine failure, not a reason to expand bounds or conceal the last proposal.
9CQD demonstrates that the missing post-assembly integration, not an inherent
ban on moving observed anchors, prevented this generated-only P2 route from
passing. This single saved-assembly control is not prospective qualification.

Evidence: `tmp/portfolio-anchors-v1/retained-controls-01/`, including declarations,
complete original donor/preparation evidence, failed stage-1 coordinates, trial
coordinates/reports, final anchor audits and the accepted 9CQD bundle. The old
P2 and original cohort results remain unchanged; the combined route has an
accepted result for nine of the ten retained structures, with 1HE8 still open.

A further separate-process audit confirmed the 9CQD bundle's readback and exact
resume with unchanged bytes and **532 generated / zero observed atoms** in its
DatRecord. It independently refused 1HE8 publication before destination
reservation and verified unchanged historical summary digests and selected/local
coordinates. Evidence: `tmp/portfolio-anchors-v1/independent-readback-02/summary.json`.
The first extra audit attempt hit the outer 120-second orchestration timeout;
its directory remains retained separately. The successful second audit used no
new refinement or sampling and did not overwrite that interrupted attempt.

## A disjoint, private publication format

`portfolio_contract.py`, `portfolio_evidence.py` and `portfolio_publication.py`
define schema-1 `experimental-gap-portfolio-bundle`, not schema-4 single-seed
sampling and not atom-completion publication. Existing readers reject it.
Index/evidence contracts cannot claim scientific qualification or a public
profile. Each whole gap keeps its actual donor candidate and original seed.

The bundle retains:

- the exact frozen common request and normalized source bytes;
- every donor's stage-1 coordinates, raw proposal, original runner result,
  staged sampler trace and bound refinement evidence;
- the selected coordinate projection, whole-gap choices, bounded search
  protocol/report and selection source snapshot;
- any issue-local repair report and selector/kernel/conditioning snapshots;
- original/prepared input bytes, original repair report and frozen request
  snapshot exactly when the common request has preparation bindings;
- the final candidate, exact generated-only `candidate.dat` and an independently
  recomputed validation receipt.

Historical traces/evidence remain byte-identical, including their old internal
paths. The private manifest relocates its own flat references; readback matches
the historical links by digest instead of traversing old runtime paths.
Prepared snapshots need no rebinding because P2 preserves the common donor
request rather than silently changing its policy.

Publication audits live input custody, source preparation, lineage, mobility,
geometry and global L-chirality before reserving a destination. It durably stages
flat payloads and re-audits those exact bytes, then exclusively reserves the new
destination and commits `bundle.json` last using `commit_manifest_last`. A crash-
left directory is not a result and is never adopted/overwritten automatically.
Original run artifacts and failed staging diagnostics remain available.

Readers enforce bounded single-linked regular files, exact flat inventory,
digest bindings, strict versions/kinds and recomputed geometry/mobility. Neither
updated hashes nor an edited `passed` receipt grants scientific authority.
Every donor stage must itself be globally all-L and exactly a generated-coordinate
projection of its bound raw proposal, including when that donor's affected gap
is not ultimately selected. Donor generated metadata must also agree. Two
adversarial tests exposed the initially missing independent unused-donor checks;
the hardened auditor now rejects both cases even with a rehashed trace chain.
The `.dat` uses `DatRecord` and contains exactly `request.generated_atoms`, keyed
by case-sensitive chain/residue/insertion-code/atom identity. Template-built
atoms in the prepared conditioning source are not relabelled diffusion-generated.

### Exact read-only resume

`resume_portfolio_bundle` independently reads a complete bundle, then compares an
identity built from every artifact digest, parent label/seed, gap choice,
operation and software snapshot. Paths may be relocated; inputs/lineage/budgets
may not. It invokes neither selection, sampling nor minimization and does not
rewrite payloads or receipts. A directory without a valid committed index fails.

The research entry point is deliberately separate from `dvbfixer model`:

```bash
python scripts/publish_gap_portfolio.py --portfolio <selection> \
  --repair <optional-local-repair> --cohort <retained-cohort> --case 1KTZ \
  --work <new-evidence-workspace> --bundle <new-bundle-directory>
# Exact independent reuse; a fresh evidence workspace is still reserved:
python scripts/publish_gap_portfolio.py --portfolio <selection> \
  --repair <optional-local-repair> --cohort <retained-cohort> --case 1KTZ \
  --work <another-new-workspace> --bundle <complete-bundle> --resume
```

Omit `--repair` for an already clean selected portfolio. No neural sampler or
refinement is run by either publication mode.

### Executed publication controls

| Artifact | Generated `.dat` atoms | Preparation custody | Readback / exact reuse |
|---|---:|---|---|
| 1KTZ multi-parent + local repair | 710 | no preparation binding | passed, including a second process |
| 1BRS homogeneous selection | 31 | independently re-audited | passed |
| 1AHW homogeneous selection | 272 | no preparation binding | passed |
| bounded 1HE8 selection | — | — | geometry rejected before accepted bundle publication |

A rehashed mutation of 1BRS's original preparation bytes was rejected by the
source-binding audit even after updating evidence and resume-key hashes. Original
bundles were unchanged. Unit controls cover forged `.dat`/validation receipts,
expanded mobility, execution budgets, unknown kinds, symlinks/hardlinks, duplicate
JSON keys, payload bounds, index-commit failure and orphan destinations.

A separate actual CLI smoke copied the 1KTZ PDB/DatRecord outside the immutable
bundle and ran `minimize --no-solvent --platform CPU --max-iter 50 --atom-naming
standard`. It loaded/resolved all 710 generated atoms, produced all-L output,
and left the bundle unchanged. This is a downstream integration smoke, not a
replacement P2 candidate or a new fixed-coordinate acceptance result.

Evidence: `tmp/portfolio-publication-v1/`, including `controls-01/`,
`preparation-mutation-01/` and `minimize-handoff-01/`. A further fixed-budget
historical ten-case execution completed in `retained-panel-01/`:

| Historical case | Private P2 workflow |
|---|---|
| 1BRS, 1G9M, 1GC1, 1JPS, 1I2M, 1AHW | Existing generated-only donor selected; bundle/readback/resume passed |
| 1KTZ | Multi-parent selection + issue-local repair; bundle/readback/resume passed |
| 9CQD | Generated-only portfolio/local repair still rejected (amide/Ramachandran) |
| 1HE8 | Generated-only portfolio/local repair still rejected (clashes/junction/Ramachandran) |
| corrected 1N8Z | Separate composed-input-preparation software arm; prepared-source bundle/readback/resume passed |

There are seven private passes in the nine-case baseline software arm, plus one
in the separate 1N8Z arm. These are **one selected portfolio workflow per case**,
not 50 new samples. The original 27/50 seed outcomes remain unchanged. The 9CQD
negative result does not contradict its prior bounded-anchor P0/P1 success:
this P2 version admits only generated-only donors and repairs, not moved anchors.

After freezing/executing that panel, the independently hardened donor auditor
read all eight completed bundles again and performed exact read-only resume.
All eight passed with unchanged bundle bytes; the two rejected cases still had
no accepted bundle. `retained-panel-revalidation-01/` records original publication
versus stronger verifier source hashes separately, without rewriting the frozen
panel or inventing a repeated scientific execution.

The final fast suite passed **2,269 tests**, with three existing optional
Torch/Biotite skips and 33 slow tests deselected. All 31 focused publication and
symmetry tests passed; Ruff and agent-map validation passed. Slow integrations
and GUI tests were not rerun for this private Python route.

## Factorial anchor analysis without relaxing limits

[The new declaration](anchor-calibration-factorial-v1.json) specifies 36 separate
backbone/Cα/heavy triplets on the original five strata and both leave-anchor-out
frames. Raw inputs are verified against original download custody. Every exact-
name measurement is recomputed before analysis. Source bond/chirality QC is
recomputed too; the defective original 1LYZ stratum is retained and flagged.

`anchor_symmetry.py` adds a conservative label-sensitivity metric. Permutations
must be automorphisms of the canonical Amber heavy graph preserving element,
atom type, charge and bonded-H degree. N/CA/C/O/OXT/CB names are fixed. Atom-type
or protonation ambiguity is not guessed: e.g. canonical HIS without a unique
tautomer template is explicitly unavailable. Inputs are never relabelled or
optimized. The graph/template digest, complete mapping, exact maximum and
minimum possible maximum under permitted permutations are reported together.
This is not a new atom-naming policy, a full chemical-equivalence oracle, or a
production movement exception.
In particular, a label-sensitivity control between deposited experimental
models does not authorize a runner to swap fixed input identities or relabel
its optimized coordinates to pass the existing movement audit.

For a fair comparison, unavailable templates are excluded only from the **paired
metric comparison**, with denominators reported; exact-name results remain for
all original measured residues.

### Current-grid coverage, local leave-anchor-out frame

| Stratum | Paired residue count | Exact names | Added symmetry-sensitive metric |
|---|---:|---:|---:|
| 1UBQ–1UBI X-ray pair | 71 | 60 | 69 |
| 1D3Z models 1–2, one NMR ensemble | 71 | 54 | 54 |
| 1UBQ–1D3Z cross-method | 71 | 33 | 38 |
| original 1LYZ–193L, source-QC flagged | 124 | 71 | 75 |
| 1AKI–193L follow-up | 124 | 98 | 99 |

These comparisons retain the **same** 0.5/0.5/1.25 Å limits. The X-ray ubiquitin
pair shows that some apparent heavy-atom movement is atom-label ambiguity, not
evidence for a larger physical displacement radius. The within-NMR result does
not improve under those same automorphisms; substantial conformational variation
is not universally explained by label swaps.

Separately increasing only backbone RMSD from 0.5 to 0.75 Å, leaving Cα and heavy
ceilings fixed, adds no jointly passing residues in the clean X-ray paired
controls above. This is sensitivity of these constrained triplets, not proof
that backbone motion is negligible or that 0.5 Å is universally calibrated.

Raw exact contact changes remain descriptive controls: an improved symmetry
metric alone does not establish meaningful contact preservation, density support
or energetic validity. The five strata are not pooled, NMR models are not
independent experiments, and two proteins are not broad calibration. No
production ceiling, request revision, historical byte sequence or B-factor
interpretation changes here.

Evidence: `tmp/anchor-calibration/factorial-symmetry-v1/`. Larger independent
experimental/interface panels and prospective policy qualification remain open.
See the [canonical checklist](../plans/implementation-checklist.md) and
[implementation plan](../plans/diffusion-geometry-repair.md).
