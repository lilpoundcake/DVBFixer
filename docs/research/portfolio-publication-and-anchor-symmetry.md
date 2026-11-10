# P2 publication, exact reuse and anchor label sensitivity

Date: 2026-11-10. Continuation of the
[private P2 and calibration controls](gap-portfolios-and-anchor-calibration.md).
P3 is explicitly excluded. Public CLI/Homology/GUI portfolio admission and higher
production anchor ceilings remain disabled.

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
