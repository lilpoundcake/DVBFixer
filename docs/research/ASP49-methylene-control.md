# Local CH₂ initialization prevents ASP49 inward HA but creates a different failure

Status: **prospective one-protein initial-condition mechanism control**, not
production repair or whole-system acceptance. Follows the
[saved-only ASP49 force attribution](graph-methyl-ASP49-attribution.md). All
production parameters, safeguards, earlier adverse results and source refusals
remain unchanged.

## Frozen matched intervention

[Complete independent evidence](ASP49-methylene-control-evidence.json) retains
**18/18 completed cells**, six existing 1TIM contexts × unchanged / local LEU48
CH₂ / all-carbon CH₂. The six contexts are all three earlier methyl source arms
and both existing seeds, not six independent proteins. There are no new seeds,
retries, timeouts, execution blockers, replacements or budget nonexecutions.

The local arm changes **only initial LEU48 HB2/HB3 coordinates**. Both are newly
generated, mobile H on CB with heavy neighbors CA/CG. It never rewrites ASP49 HA,
an observed/kept H, any heavy atom or a force. The all-carbon arm applies the same
construction to all **310** eligible CH₂ parents in each context. No methyl group,
nitrogen H group, force parameter or ownership mask is silently broadened.

Construction is fixed before sampling. The two lexicographically ordered heavy
bond unit vectors define a bisector `b` and cross-product normal `n`. For observed
heavy angle θ, the H directions are `a*b ± sqrt(1−a²)*n`,
`a = −1/(3*cos(θ/2))`; the first lexicographic H takes +normal. Each **original
C–H radius** is preserved. Both H–C–heavy cosines are −1/3 while the heavy–C–heavy
angle remains unchanged. **H–C–H is implied**, not independently forced to an
ideal angle when heavy geometry is nonideal. Undefined/collinear frames,
`cos(θ/2) <= 1/3` or degenerate H separation are refused, without fallback or
heavy-coordinate repair. No energy/native-based slot or phase selection occurs.

Every arm reuses its original H System and **Reference / tolerance 1.0 kJ/mol/nm /
maximum 50 iterations**, with a fresh minimizing Context. All six unchanged
endpoints exactly reproduce the prior source-arm endpoints in nm. Independent
auditing checks ownership, both geometric cosines, slot sign, radii, unchanged
nonowned starts and all fixed-atom trajectory coordinates. All mobile-H forces
equal negative native callback gradients at **900 callbacks**, maximum component
error **0.0 kJ/mol/nm**. All **936 initial/callback/final nodes** remain reported.

Plan SHA256:
`36131c49174fb83ef637b908ba4d8bd2399a11dbf3ac28ffbca85d0673b4f74f`.
Private plan/run/audits are `ASP49-methylene-control-*` below
`tmp/pdbfixer-research/`. Budgets are 1,800 seconds per cell and 7,200 seconds total.
No candidate PDB or production output is written; raw runtime paths do not
transfer through Git.

## Focal effect and full denominator

| H-stage descriptor | Unchanged | Local LEU48 CH₂ | All-carbon CH₂ |
|---|---:|---:|---:|
| Executed matched contexts | 6 | 6 | 6 |
| Any ASP49 inward callback | 3 | **0** | **0** |
| Final ASP49 inward HA | 2 | **0** | **0** |
| Final global all-acute HA atoms | 2 | **1** | **0** |

The local two-H intervention prevents the focal ASP49 endpoint and transient
crossings in every matched context without rewriting ASP49 HA or changing native
forces. This gives controlled support for a **local initial-condition mechanism**,
beyond the preceding force/path association. It is not a single-pair force
ablation, a unique-root-cause proof or a general safe initializer: changing two
initial H also changes their interactions and the subsequent global optimizer path.

## New adverse outcome: ALA151, not ASP49

In **seed 2 / graph-tetrahedral-zero source / local LEU48 CH₂**, a new final
all-acute HA appears at **1TIM/A/ALA151**. Its N/C/CB-to-HA angles are
**63.558/87.155/71.232°**, CA–HA **1.073012 Å**. Its original starting HA angles
were 86.193/109.615/120.729°, and its initial coordinate is untouched by the local
intervention. The matched unchanged source arm did not have this final failure.

This is spatially distant: original LEU48 CA/CB are **27.429/26.700 Å** from
ALA151 CA. The two newly placed local H are **26.450/25.783 Å** from initial
ALA151 HA, beyond the native **10 Å** nonbonded cutoff, with no native pair
exception between them. That excludes a direct initial nonbonded force from those
two H on this distant HA; it does **not** distinguish indirect environment
coupling from the global quasi-Newton/line-search history. Preserve the regression
rather than calling a local fix globally safe or claiming nonlocal physical action.

The all-carbon CH₂ arm avoids this particular final all-acute count, but it also
rewrites **620 generated H**, not two, and is not an isolated local-causality test.

## Large H strain remains in every construction

Maximum final H-bond deviations range across the six contexts:
**0.133–0.228 Å** unchanged, **0.157–0.195 Å** local, **0.149–0.167 Å** all-carbon.
Maximum H-angle deviations reach **72.281/70.937/71.256°**, respectively. For seed
1 / unchanged methyl source, that maximum grows **52.365→70.160/71.256°** under
local/all-carbon construction despite the rescued focal HA. Those adverse
descriptors remain visible; zero focal crossings or zero final all-acute H is not
complete H-geometry acceptance or convergence at the fixed iteration cap.

No heavy atom moves in this H-stage comparison, so the later ASP49 or ALA151 Cα
geometry cannot be inferred. The next separately frozen whole-native-phase
continuation must include **all 18 cells**, reproduce all six corresponding prior
whole endpoints exactly in the unchanged arm, retain the original H-label repair/
object-index bridge, and report both sites plus global H/heavy outcomes and
original-observation motion. No extra seed, iteration, phase or favorable-case
selection is justified by this mixed H-stage result.

Validation: **22 focused geometry/native-worker/independent-audit tests** pass,
including a later explicit nonideal-heavy/implied-H–H test. Full fast suite through
the initial control/audit implementation: **2,999 passed**, three optional
Torch/Biotite skips and 33 slow deselected; that later test also passes separately.
Ruff passes and independent publications regenerate byte-identically. Original
THR2 Cβ causality, broader qualification and empirical-table production provenance/
license/approval remain unresolved; production/public repair stays disabled.

Code: `scripts/research_methylene_geometry.py`,
`scripts/research_ASP49_methylene_control.py`, `scripts/audit_ASP49_methylene_control.py`.
