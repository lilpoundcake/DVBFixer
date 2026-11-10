# Rule-selected 48-mask comparison and native-reference gate control

Status: **predeclared, own-success-unfiltered case study**, not broad scientific
qualification. Follows the [retrospective matched comparison](protein-completion-template-comparison.md).
Production, source policies, thresholds, seeded retries and chirality guards are
unchanged. Only research harnesses/tests were added in this continuation.

## Selection frozen before reconstruction

`scripts/unfiltered_branched_comparison.py` separates `plan` and `run`. Planning
selects the **first two VAL and first two THR identities in chain-A file order**
from each of 1EMV, 1TIM and 1DQJ, with four masks per residue:

1. Both distal branches missing, observed CB.
2. Both distal branches and CB missing.
3. First distal atom missing, observed CB/partner.
4. Second distal atom missing, observed CB/partner.

All **48 masked inputs**, complete declared ATOM projections, originals, source
hashes and case identities were written under
`tmp/pdbfixer-research/unfiltered-plan-01/` before any proposal execution. The
runner verifies these hashes and cannot silently substitute changed code/inputs.
Planning never calls admission or generation; tests enforce this. Runner tests
use mock errors rather than previewing real reconstruction outcomes.

Selection does not filter by occupancy, altloc, source stereo or own reconstruction
success. The chain-A projections retain **every original ATOM record**, including
unsupported observations, with explicit bare TER/END. No occupancy rewrite, altloc
choice, source repair or hidden cropping. These are **not intact complexes**:
other chains and heterogens/waters are not in the controlled context. 1EMV is
homologous to previously tested 1FR2; these files are not independent protein
families. This is deterministic rule selection, not randomized/blind validation.

Both arms share admission and exactly the same admitted source. Own rejection
does not prevent template execution. Environment errors stop preflight; later
execution/protocol errors would be retained separately, never counted as geometry
refusals. The actual run had **zero template execution/protocol blockers**.

## Executed outcomes

[Full comparison evidence](protein-completion-unfiltered-evidence.json), retained
under `tmp/pdbfixer-research/unfiltered-comparison-01/`:

| Projection and mask class | Planned | Shared admission | Own common pass | PDBFixer common pass |
|---|---:|---:|---:|---:|
| 1EMV paired | 8 | 8 | 8 | 6 |
| 1EMV singleton | 8 | 8 | 8 | 6 |
| 1TIM paired | 8 | 8 | 8 | 2 |
| 1TIM singleton | 8 | 8 | 7 | 4 |
| 1DQJ paired | 8 | 0 | Not executed | Not executed |
| 1DQJ singleton | 8 | 0 | Not executed | Not executed |
| Total | **48** | **32** | **31** | **18** |

All 16 1DQJ masks refused shared admission because source altloc/partial occupancy
was retained. They remain part of the declared denominator, not dropped masks or
16 scientific proposal failures. Of 32 admitted masks: both arms pass on 17,
own alone on 14, template alone on one, neither on zero.

On the **same mutually passing subset**, exact-name generated-coordinate median
RMSD remains a separate question:

| Class | Shared subset | Own median | PDBFixer median |
|---|---:|---:|---:|
| Paired | 8 | 0.184 Å | **0.130 Å** |
| Singleton | 9 | **0.055 Å** | 0.137 Å |

Reference coordinates were consulted only after both proposals/selection. No
truth-based reranking, additional refinement or reflection was used. The gate
still uses prototype illustrative angles and a coarse steric threshold: it is
**not an independent calibrated geometry/packing standard**. This run does not
establish general superiority or full legacy prepare/minimize performance.

## Native-reference sanity control exposes a gate limitation

`scripts/audit_unfiltered_reference.py` subsequently reinserts only withheld native
coordinates for every admitted mask and runs the **unchanged** common gate, with
no new reconstruction or threshold tuning. The
[native-reference audit](protein-completion-unfiltered-reference-audit.json)
shows **29/32 native controls pass**. All three failures are masks exposing the
same **1TIM VAL39 CG1–backbone C** contact, a topological **1–4 pair**:

- Native distance: **2.762756 Å**, correctly assigned CA=S.
- Our singleton proposal distance: **2.767208 Å**, native error only **0.018689 Å**,
  but rejected by the common C–C 2.8 Å cutoff (VDW sum minus 0.6 Å).
- The template singleton passes the gate, with native error **0.198586 Å**.

Thus this is not evidence that a 0.019 Å reconstruction is worse than a 0.199 Å
one, nor three independent source stereo defects. It reveals a **local steric-gate
limitation**, shared by both arms. The gate excludes bonded/angle-connected pairs
but treats the remaining 1–4 contact by the same coarse VDW cutoff. Do not silently
relax that threshold, exclude every nearby pair or reinterpret native failures
as reconstructed D/allo chemistry. No policy was changed after seeing this result.

The native audit also reports all own attempts against truth **without choosing
among them**. A rejected one-atom proposal remains rejected in the experiment.
This control is a reference-relative diagnostic, not calibrated external validation;
independent chemistry-aware packing/geometry qualification remains open.

## VAL symmetry and labels

Reporting now includes both exact-name RMSD and a generated-only methyl-swap RMSD.
Swapping VAL CG1/CG2 is allowed **only when both methyl atoms are generated**;
observed CB and partners are never exchanged or moved. Singleton masks cannot
borrow a fixed partner to improve RMSD. THR OG1/CG2 are chemically distinct and
never permuted. Named χ1 error is retained separately.

No executed proposal improved under the allowed swap, so symmetry-minimized
medians here equal exact-name medians. Unit tests use deliberately swapped methyls
to verify that reporting normalization works and cannot change fixed ownership.
Neither scoring option changes candidates, named-orientation gates or selection.

## Verification and reproduction

**556 focused tests passed in 24.42 s**; seventeen-file package mypy and targeted
Ruff passed. Tests enforce outcome-independent planning, frozen-input checks,
baseline execution despite own failure, fixed-partner symmetry constraints and
the zero-RMSD native-reference gate failure without mislabelling it wrong stereo.
Historical controls, hashes and thresholds remain unchanged.

```bash
# Scientific environment required; use distinct exclusive destinations.
PYTHONPATH=src:. python scripts/unfiltered_branched_comparison.py plan --output tmp/pdbfixer-research/new-unfiltered-plan
PYTHONPATH=src:. python scripts/unfiltered_branched_comparison.py run --plan tmp/pdbfixer-research/new-unfiltered-plan --output tmp/pdbfixer-research/new-unfiltered-run --summary tmp/pdbfixer-research/new-unfiltered-results.json
PYTHONPATH=src:. python scripts/audit_unfiltered_reference.py --results tmp/pdbfixer-research/new-unfiltered-results.json --output tmp/pdbfixer-research/new-native-audit.json
```

Next: independently specify/test chemistry-aware 1–4 steric handling using native
and deliberately perturbed controls, without erasing genuine generated clashes.
Version any changed admission/audit policy and retain this unchanged baseline.
Broader families/masks, calibrated external metrics, ILE/multi-chi, licensed
rotamers and production consumer integration remain separate unfinished gates.
