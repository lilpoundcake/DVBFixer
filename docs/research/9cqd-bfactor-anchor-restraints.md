# 9CQD: B factors and restrained anchor relaxation

Date: 2026-10-27. Research branch: `research/9cqd-anchor-relaxation`.
Status: all twelve minimizations completed; not a production policy.

## Question and literature

Can positional restraints on the immediate observed anchors repair gap geometry
with less source drift than the earlier unrestrained anchor minimization? Does
B-dependent weighting improve that trade-off over equal weights?

There is no established universal B-factor threshold that authorizes moving
loop anchors. B factors include local and collective disorder, refinement and
data-quality contributions; they are not coordinate-error estimates or allowed
displacement radii. The conversion `Uiso = B/(8*pi^2)` describes isotropic ADPs,
not a justified coordinate restraint width. Carugo's independent-structure study
estimates B reproducibility around 6–9 Å² in its lysozyme dataset; this is not
an error bar calibrated for individual 9CQD atoms or its 3.10 Å resolution.

Sources accessed 2026-10-27:

- [Carugo, 2022, B-factor accuracy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8725162/).
- [PHENIX refinement documentation](https://phenix-online.org/documentation/reference/refinement.html):
  ADP decomposition, coordinate refinement and geometry/data targets.
- [RCSB entry metadata](https://data.rcsb.org/rest/v1/core/entry/9CQD):
  X-ray diffraction, 3.10 Å, deposition-wide mean B 74.5 Å², Wilson B 67.15 Å²,
  collection at 100 K. These entry-wide values are not the selected-complex mean.

Relative backbone B within a comparable chain/domain is a useful ranking signal.
P75–P90 and above-P90 are proposed exploratory priority bands, not literature
standards, automatic mobility selections or acceptance thresholds. Domain/TLS
motion, occupancy, alternative conformations, map support and local chemistry
must also be considered. High B alone is insufficient; low B does not prohibit
repair of an independently established model error. B columns containing pLDDT,
placeholder values or nonexperimental predictions cannot use this interpretation.

## Source observations

Input: `tests/fixtures/protein_context/9cqd_split.pdb`, SHA-256
`6bcc9638c3d299ea12151c3894b7e6c3482ebd2fdfe1d8e54dcc1f6f6316a281`.
Only observed ATOM N/CA/C/O are counted; all 457 residues have all four atoms.
Units are Å². Atom-weighted statistics:

| Chain | Residues | Minimum | Mean | Median | Maximum | Population SD |
|---|---:|---:|---:|---:|---:|---:|
| All | 457 | 35.09 | 73.64 | 70.60 | 140.70 | 18.32 |
| M | 27 | 37.68 | 65.35 | 63.95 | 101.73 | 12.09 |
| W | 211 | 35.09 | 76.75 | 74.02 | 128.92 | 17.91 |
| X | 219 | 35.15 | 71.66 | 67.76 | 140.70 | 18.76 |

Across all atoms P5/P25/P75/P95 are 49.02/60.76/83.44/110.25.
Fractions in <40, 40–60, 60–80, 80–100, 100–120 and ≥120 are respectively
0.66%, 22.76%, 47.05%, 20.02%, 6.89% and 2.63% (rounded).

Residue statistics below first average N/CA/C/O within each residue. Percentile
is the fraction of observed residues in that chain with mean B <= the anchor's.

| Anchor | Mean B | Chain percentile | Chain P75 | Chain P90 |
|---|---:|---:|---:|---:|
| W25 SER | 88.73 | 77.73 | 86.34 | 100.92 |
| W29 VAL | 89.99 | 79.62 | 86.34 | 100.92 |
| X139 THR | 97.33 | 91.78 | 79.05 | 94.58 |
| X142 GLY | 94.00 | 89.04 | 79.05 | 94.58 |

Anchor occupancy is 1.00. W anchors are elevated but not extreme. X anchors
rank higher within their chain. Neighboring W30–33 average about 98–114;
X135–138 about 89–102. After X142 the means decrease to X143 71.73,
X144 58.86, X145 51.15. This argues against indiscriminately releasing ±3
observed residues. Domain-specific normalization remains future work.

## Earlier evidence and scope

The [earlier factorial study](cc89-lrd-ablation-9cqd.md) found geometry failures
with exact anchors and geometry passes in all six generated+anchor OpenMM arms,
but all six failed strict source-coordinate preservation. Raw → anchor-mobile
OpenMM without LRD had W-aligned combined anchor backbone RMSD 0.512 Å (off)
and 0.455 Å (on). These are measured outcomes, not experimental error tolerances.
Observed-mobile LRD increased drift and did not remove the need for minimization.

## Frozen new experiment

Use both saved seed-7 raw candidates from the corrected-frame study, without
LRD or new neural sampling. Targets M/W/X 27/214/221; gaps W26–28, X140–141;
anchors W25/W29/X139/X142. Keep the original request's generated mask unchanged
for independent validation. Mobility is a separate research-minimizer argument.

Twelve new runs: two reinjection arms times six policies:

1. Generated atoms only, anchors exactly fixed.
2. Generated plus anchor heavy atoms, no positional anchor restraints.
3. Same mobile mask, uniform k = 1000 kJ mol⁻¹ nm⁻².
4. Same mobile mask, B-weighted k with mean 1000 kJ mol⁻¹ nm⁻².
5. Uniform k = 10000 kJ mol⁻¹ nm⁻².
6. B-weighted k with mean 10000 kJ mol⁻¹ nm⁻².

For each observed anchor heavy atom the research-only harmonic force is
`E = 0.5*k_i*|x_i-x_source_i|²`, with OpenMM coordinates in nm.
All existing anchor heavy atoms, not only backbone, receive restraints.
Each residue uses its original N/CA/C/O mean B. Raw weights are
`median(four anchor B)/B_residue`, clipped to [0.5, 2.0], then divided by their
mean over the 24 restrained atoms. Thus both policies have the same average k;
generated atoms have no positional restraints. No direct ADP-to-k calibration
is claimed. The narrow range of these four B values may yield only a small
uniform-versus-weighted difference. Two strengths are a predeclared sensitivity
check, not a reference-guided search.

The mean-one coefficients calculated before inspecting outcomes are W25 1.040823,
W29 1.026220, X139 0.948852, X142 0.982389. Thus the B-weighted policy is only
about 5% softer or 4% stiffer than uniform at individual residues, not an on/off
release decision. At uniform k=1000, a 1 Å single-atom displacement costs
5 kJ/mol; at k=10000 it costs 50 kJ/mol. These are chosen regularization strengths,
not experimental B-derived uncertainty estimates.

Reuse boundary refiner v8, AMBER force field, CPU one thread, no solvent,
500 maximum iterations per attempt, up to eight restarts, perturbation 0.25 Å,
seed 7 and unchanged geometry-first/energy-second restart selection. Nonfixed
policies use the same generated+anchor scope for amide/angle guards and restart
scoring. The fixed control retains the original generated-only scope. Runtime
injection adds a restraint force only in a research process; production source
files, launcher behavior, validators and thresholds are not modified.

Save inputs, code snapshots and hashes before executing. Never overwrite earlier
experiments. Report failures as failures. Check nonmobile atom records byte for
byte. Report strict validation separately from geometry with source-drift gates
excluded; inspect all mobile residues plus one sequence neighbor at outer edges,
including clashes against the remainder of the complex. Report per-anchor
backbone RMSD and maximum heavy-atom movement, in the original frame and after
one W-chain observed-backbone Kabsch fit applied to the whole complex. Generated
atoms have no experimental source counterpart: these are preservation metrics,
not missing-loop accuracy. Require finite all-L coordinates.

The W fit includes observed anchors, excludes generated atoms, and is shared
with chain X (no independent X alignment). RMSDs are calculated before PDB
re-rounding of the aligned view. Literal source-frame displacement remains the
metric relevant to the unchanged production preservation contract.

## Results

All 12 runs completed without execution errors (about 377 seconds summed run
time). All 10 mobile-anchor candidates pass the tested geometry, including the
expanded neighborhood and outer junctions. Both fixed controls fail amide
planarity at W29. All candidates preserve nonmobile atom records exactly.
Zero candidates pass strict original-coordinate acceptance: both fixed controls
fail geometry, and all mobile controls fail the two unchanged source-drift gates.
No experimental-density validation is included in this pilot.

Combined RMSD is over 16 anchor N/CA/C/O atoms after the same W-chain fit;
maximum movement is over all 24 anchor heavy atoms in the original frame.
Units are Å; k is kJ mol⁻¹ nm⁻².

| Reinjection | Policy | Mean k | Combined anchor RMSD | Maximum anchor movement | Geometry pass |
|---|---|---:|---:|---:|---|
| Off | Fixed | — | 0.0000 | 0.0000 | No |
| Off | Free | 0 | 0.5015 | 1.6448 | Yes |
| Off | Uniform | 1000 | 0.4519 | 1.5743 | Yes |
| Off | B-weighted | 1000 | 0.4547 | 1.5776 | Yes |
| Off | Uniform | 10000 | 0.2379 | 0.8240 | Yes |
| Off | B-weighted | 10000 | 0.2326 | 0.8306 | Yes |
| On | Fixed | — | 0.0000 | 0.0000 | No |
| On | Free | 0 | 0.4804 | 1.9920 | Yes |
| On | Uniform | 1000 | 0.4346 | 1.5950 | Yes |
| On | B-weighted | 1000 | 0.4352 | 1.6127 | Yes |
| On | Uniform | 10000 | 0.2324 | 0.8838 | Yes |
| On | B-weighted | 10000 | 0.2324 | 0.9085 | Yes |

For the strongest tested positional restraints, per-anchor backbone RMSD:

| Reinjection | Policy, k=10000 | W25 | W29 | X139 | X142 |
|---|---|---:|---:|---:|---:|
| Off | Uniform | 0.1944 | 0.4080 | 0.1228 | 0.0838 |
| Off | B-weighted | 0.1907 | 0.3950 | 0.1207 | 0.0973 |
| On | Uniform | 0.1821 | 0.3982 | 0.1173 | 0.1025 |
| On | B-weighted | 0.1715 | 0.4052 | 0.1285 | 0.0771 |

### Interpretation

Increasing positional stiffness from 1000 to 10000 reduces drift substantially
without losing the tested geometry passes. This supports studying restrained,
not fully free, immediate-anchor motion. It does not identify a maximum safe
stiffness, a B threshold, an optimal schedule, or experimentally accurate loops.

B weighting has no consistent advantage here: at k=1000 combined RMSD is slightly
higher in both arms; at k=10000 it is lower by about 0.0053 Å off and nearly
unchanged on. Maximum heavy-atom movement is slightly higher than uniform in
all four matched comparisons. With one source, one seed, narrow B variation and
no replicate uncertainty estimate, these differences do not justify promoting
the weighted policy. W29 remains the dominant backbone displacement even though
the X anchors have larger relative B values. B is not a reliable predictor of
the displacement actually needed to close a particular junction.

The new free controls qualitatively reproduce the earlier geometry/drift trade-off,
but not identical RMSDs (earlier 0.512/0.455 versus new 0.502/0.480 Å). The new
fixed-on control has one amide finding, versus two in the historical run.
The archived candidates remain unchanged; the new matrix uses its own freshly
rerun controls. Seeded protocols are not asserted to give bit-identical minima
across scientific/runtime revisions or executions. The recorded runtime is
OpenMM 8.5.2, PDBFixer 1.12.0, NumPy 2.4.6.

### Evidence and reproducibility

Research code: [runner](../../scripts/cc89_anchor_restraints.py) and
[unit tests](../../tests/test_cc89_anchor_restraints.py). Run from repository root:

```bash
micromamba run -n dvbfixer python -m scripts.cc89_anchor_restraints \
  --baseline test/errors_diffusion/9cqd-lrd-factorial-seed7-fixed-frame \
  --output test/errors_diffusion/9cqd-anchor-restraints-seed7-v1
```

The existing output directory is never overwritten; a replay needs a new output
path. This script creates no solvent box and never publishes accepted artifacts.
The ignored evidence directory retains all twelve PDBs, per-run JSON, full
comparison, source/request copies, launch-time code snapshots and manifest.
Full candidate SHA-256 bindings and compact measurements are also retained in
the tracked [summary table](9cqd-anchor-restraints-seed7-v1.tsv).

| Evidence file | SHA-256 |
|---|---|
| manifest.json | bd49e1bb40cf78380190e40933770b4057704b1a7fba92518b35ee46a43a5150 |
| comparison.json | e977fa5f63517483c9d8296ba14806e250438f81f8fead149c713e8cf86ff943 |
| completed.json | 86f4df17ca1cd1b9d791293c384c1d9efb61766daf3a73372f9993d3bc520387 |

The launch-time script snapshot precedes non-scientific import/style cleanup
and explicit closure binding in the working script. Runs used a synchronous
closure within each iteration; the cleanup does not alter the declared force,
mask, strength or minimization settings. The launch snapshot and its hash,
rather than an assumption that current source bytes are identical, are evidence.

## Visual-review export

All 35 review PDBs are saved under the ignored directory
`test/errors_diffusion/9cqd-anchor-restraints-visual-review-v1/`: one source,
14 original-frame outputs (two raw plus twelve minimized), 14 W-aligned copies,
and six motif files. Each gap has a source-identity conditioning subset and
both byte-identical archived engine motif inputs. Engine motifs use proxy chain
IDs/local numbering; their coordinate/atom-name/residue-name/occupancy/B content
was verified against the source-identity subset, independent of relabeling.

Use source plus `kabsch-W/` for overlays; use `original-frame/` to inspect literal
source-coordinate preservation. The export includes `structures.tsv`, full
metrics, rigid transforms, context spans, a README and `SHA256SUMS`. All 43
exported-file checksums were independently verified. Review PDBs are not added
to tracked fixtures. The [exporter](../../scripts/export_cc89_anchor_review.py)
refuses to overwrite an existing output directory and runs no new modeling.

## Implementation pilot: geometry-gated stages and chain-normalized B

Executed 2026-10-27. This is a separate research pilot, not a qualified default
change. See the [handoff](../plans/geometry-gated-anchor-relaxation-handoff.md).
The new core provides pure B provenance/fallback and per-anchor movement checks,
a staged controller and an explicit observed-anchor restraint API. Production
runner/CLI/GUI and publication acceptance have not yet been migrated.

The [pilot script](../../scripts/run_cc89_staged_anchor_pilot.py) independently
starts both minimizations from raw. It compares uniform k=10000 with
`k=10000*clip(chain_backbone_median/anchor_backbone_mean,0.5,2)`, releasing only
residual-geometry-named immediate anchors. Eight complete chain backbone samples
is an engineering pilot population criterion, not a calibrated scientific cutoff.

All four stage-2 outputs pass the tested generated/anchor/outer geometry checks;
all still fail the old strict source-coordinate gates. No density agreement or
production anchor-drift qualification is claimed. Units below are Å.

| Reinjection | Policy | Selected anchors | W25 backbone RMSD, original frame | W29 backbone RMSD, original frame | Max anchor heavy movement |
|---|---|---|---:|---:|---:|
| Off | Uniform | W29 | 0.0000 | 0.4271 | 0.8256 |
| Off | Chain-B | W29 | 0.0000 | 0.4819 | 0.9023 |
| On | Uniform | W25, W29 | 0.1932 | 0.4179 | 1.0369 |
| On | Chain-B | W25, W29 | 0.1821 | 0.4375 | 0.9379 |

X139/X142 and every source atom record outside the explicit mobile masks remain
exact. Every attempt's input digest equals the original raw digest. Chain-B
coefficients are 8429.6985 for W25 and 8311.4321 for W29 (kJ mol⁻¹ nm⁻²).

Combined W-aligned RMSD over all 16 anchor-backbone atoms is 0.2120/0.2390 Å
for off uniform/Chain-B and 0.2295/0.2357 Å for on uniform/Chain-B. These combined
values include untouched anchors: use per-released-anchor metrics for acceptance.
Comparison against an earlier all-four-mobile mask is not a controlled test of
B weighting alone. There is again no consistent B advantage: off increases drift;
on slightly increases W29 backbone drift but reduces maximum heavy displacement.

These results do not establish safe movement bounds. A 1 Å maximum inferred from
earlier 9CQD results would reject the new uniform-on output; do not silently tune
such a limit to admit a preferred candidate. Qualification remains blocked on a
broader cohort and a versioned observed-mobility contract.

Evidence: `test/errors_diffusion/9cqd-staged-anchor-pilot-v1/` (ignored), retaining
raw/both stage outputs, full diagnostics, masks, coefficients, source/request and
launch-time code snapshots. Publication is explicitly disabled.

| File | SHA-256 |
|---|---|
| manifest.json | 4b6956ee67450dc92cdf0b80aabdef1940531414d32a604b9d4f8194e3a8fe53 |
| comparison.json | 6c5152b14d47043675672ca3e4236e1c3c620f3304da211a6edce43f4aa57f6e |
| completed.json | db9e1c8a383ab68a3fa43132b45597d656aade99133d3d7521e4f4fde6a4ac3d |

## Independent acceptance-audit follow-up

The [research auditor](../../src/dvbfixer/model/diffusion/anchor_acceptance.py)
independently recomputes eligibility from stage-1 residual geometry rather than
trusting a runner-provided mobile mask. It preserves the historical strict
validation result and applies explicit caller-supplied per-anchor bounds in a
separate research summary. All reports retain `publication_qualified=false`.
This is not a replacement for runner trace validation or a shipped public
mobility contract. Typed stage-1 numerical-failure evidence remains future work.

Checks include source SHA binding, both stage input digests matching raw, clean
raw/stage-1 rejection of unnecessary escalation, unchanged observed stage-1 records,
unselected final records, selected-anchor identity/B/occupancy fields, observed
ANISOU/SIGATM/SIGUIJ records, full geometry, and new outer junction distances even
when topology inference omits the broken bond. Pre-existing source findings alone
do not justify releasing an anchor. Backbone RMSD and maximum heavy displacement
are tested per anchor in the original coordinate frame, never averaged over the
complex or hidden by alignment. Excessive movement is not treated as a warning.

Two explicit sensitivity settings were evaluated on the same four saved outputs:

| Limits: per-anchor backbone RMSD / heavy maximum, Å | Off uniform | Off Chain-B | On uniform | On Chain-B |
|---|---|---|---|---|
| 0.5 / 1.0 | Pass | Pass | Fail: movement | Pass |
| 0.6 / 1.25 | Pass | Pass | Pass | Pass |

Neither setting is adopted as a default or experimentally justified by these
results. The second setting is not a correction of the first to obtain a pass.
It illustrates the importance of predeclared bounds; the cohort remains too small
to establish general permissible motion.

The [audit CLI](../../scripts/audit_cc89_staged_anchors.py) requires both limits
explicitly, verifies saved coordinate SHA bindings and creates a new evidence
directory without rerunning minimization. Results are retained under
`test/errors_diffusion/9cqd-staged-anchor-audit-bb05-max10-v2/` and
`test/errors_diffusion/9cqd-staged-anchor-audit-bb06-max125-v2/`. Audit report hashes:

- 0.5/1.0: `848688421d862a0f4343d05c99ce9ceceddc19298ad228a79ef35e4d94285296`.
- 0.6/1.25: `5c93b2e65b0f4546df25d6ec1d7356703839fbe74a266858ddef68fadebcfd17`.

Regression evidence: 120 combined narrow tests pass, including twelve independent
audit tests. Two additional [1UBQ native controls](../../tests/test_diffusion_anchor_controls.py)
remove A6–A8 to form a source request and restore the exact experimental atoms
as a clean raw control; no neural reconstruction or minimization is performed.
Both the original experimental header and absent-method variant retain identical
raw coordinates with zero minimizer calls. The tracked 1UBQ source SHA-256 is
`d4a6812d8951cf6594e6a0763f089e35f5a80b62acb3c117b2c5565228a7b161`.
These controls demonstrate unnecessary-motion avoidance, not broader drift-bound
calibration or improved missing-loop accuracy.

The subsequent full fast-lane snapshot passed 1520 tests with two optional skips
and 31 slow tests deselected; the two 1UBQ controls added after collection also
pass separately. Ruff and documentation checks pass. These local results are
not represented as a completed remote research-branch CI run.
