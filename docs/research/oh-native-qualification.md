# Real-input O–H postflight qualification: matched replay, incomplete cohort

## Outcome

The [prospective cohort](../plans/oh-native-qualification.md), frozen at
`5d99286`, executed **all 12 declared attempts exactly once**. Eleven preparations
completed; one retained upstream naming failure prevents a complete-cohort pass.
No native attempt was repeated, replaced, extended or dropped.

All **12 saved production postflight states** replay exactly with the candidate:
same nm coordinates, return/error and repair count. Heavy coordinates, already
valid H, atom identities and explicit bonds are preserved; there are **zero new
postflight refusals**. The patch corrects **12 O–H upper-bound violations** left
by the old repair. All twelve occur on **forcefield-free fallback**, not native
AMBER/CHARMM H placement. This is useful fallback evidence, not native FF repair
efficacy or force-field accuracy qualification.

Primary scalar records:

- [Full matched audit](oh-native-qualification-evidence.json)
- [Coverage and upstream-blocker summary](oh-native-qualification-summary.json)
- [Final engineering validation](oh-native-validation.json)
- [Frozen execution plan](../plans/oh-native-qualification-plan.json)
- [Offline audit transport amendment](../plans/oh-native-audit-amendment.md)

The overall machine-readable `qualification_passed` remains **false**. An empty
postflight `blockers` list means no matched-repair defect, not 12/12 preparation
success. `main` remains untouched; no release or merge is authorized.

## Actual routes and every outcome

Numbers in the repair column include legacy nonoxygen repairs, not just O–H.
All native H minimizer boundaries reported automatic OpenCL / Apple M3 Pro /
single precision; no reporter or changed minimizer parameter was introduced.

| Case | Source | Requested FF | Successful H-placement route | Old / new repairs | Preparation |
|---|---|---|---|---:|---|
| 0 | 1UBQ + water | amber | AMBER19 | 4 / 4 | completed |
| 1 | 1UBQ + water | charmm | CHARMM36 | 7 / 7 | completed |
| 2 | 1DQJ + water | amber | AMBER19 | 2 / 2 | completed |
| 3 | 1DQJ + water | charmm | CHARMM36 | 2 / 2 | completed |
| 4 | 8CDE, pH 3 | amber | forcefield-free after C-terminal GLU template failure | 0 / 7 | completed |
| 5 | 8CDE, pH 3 | charmm | CHARMM36 | 12 / 12 | completed |
| 6 | 3RY6 | amber+glycam | AMBER14/GLYCAM | 33 / 33 | completed |
| 7 | 3RY6 | charmm | AMBER14/GLYCAM bridge, then CHARMM naming | 26 / 26 | completed |
| 8 | 1VCU | amber | forcefield-free after DAN template failure | 4 / 6 | completed |
| 9 | 1VCU | charmm | forcefield-free after DAN template failure | 8 / 11 | completed |
| 10 | GLYCAM hinge/CH3 | amber+glycam | AMBER14/GLYCAM | 0 / 0 | completed |
| 11 | GLYCAM hinge/CH3 | charmm | forcefield-free after NLN template failure | 0 / 0 | **naming failure** |

The eight successful native-FF H routes cover AMBER19, CHARMM36 and
AMBER14/GLYCAM. Four successful H calls are forcefield-free, one in an ultimately
failed preparation. Successful H placement is not successful preparation.
Native FF routes did not supply a newly stretched O–H example for this patch;
their evidence here is preservation and absence of new refusals.

## What changed geometrically

The original postflight leaves all twelve added-contract violations unresolved:

- Case 4: seven stretched O–H, including three acidic ASP/GLU hydroxyls and
  SER/TYR. One distance is just above the strict bound (1.2000098 Å); the largest
  is 1.2881384 Å. Do not move the bound to hide near-boundary findings.
- Case 8: THR/SER at 1.2167674/1.2218767 Å.
- Case 9: THR/SER/water at 1.2126643/1.2027448/1.2029470 Å.

The candidate repairs these on the **identical saved H coordinates**. Other
already-valid H and all heavy atoms stay fixed. The independent post-repair
checks find no upper-bound/nonparent-heavy violations, short-O–H (<0.8 Å)
warnings or ≤1.0 Å cross-residue-heavy contacts in these saved postflight states.
Those descriptive checks are not a general all-atom-clash or torsion validation.

## Retained upstream blocker, case 11

The CHARMM request on the already-GLYCAM-named hinge input does not use the
PDB-standard-sugar bridge. Its native FF attempt fails on NLN and the existing
pipeline runs forcefield-free H placement. The snapshots locate duplicate atom
identities **before postflight**:

- `H/UYB232/H1` and `K/UYB242/H1` are each present once and unbonded immediately
  before the fallback H call.
- Immediately afterward each identity occurs twice: a new H1 bonded to C1 plus
  the original unbonded H1.
- Candidate postflight performs **zero repairs** and preserves those atoms;
  baseline also performs zero repairs on the same state. Naming subsequently
  raises `NamingConversionError: Duplicate source atom identity`.

This supports the bounded conclusion that the observed failure predates O–H
postflight and is not introduced by this patch. It does not prove the upstream
reason the original sugar H bonds were absent or qualify a fix for that route.
No production change or compensating repeat was made.

A partial PDB exists but no valid `.dat` is published. OpenMM readback warns
about the duplicated H1 records and merges them; a superficially O–H-clean PDB
parse is **not evidence of a valid final structure**. Keep the nonzero exit,
duplicate-identity inventory, partial output and missing `.dat` visible. Never
promote the partial structure as a successful preparation.

## Auditor engineering amendment

The first auditor stopped on the legitimate A/B/C protein blocks followed by
A/B/C water blocks in case 2. Its small-graph helper merged equal chain IDs.
The original frozen source and failure record remain retained. Amendment
`6080fe0` introduced a separate offline adapter before re-auditing: ordered ID
runs, unique contiguous residue ownership, explicit unchanged bonds and exact
nm floats. It neither changes native artifacts nor re-executes sampling.

All twelve retained states satisfy that adapter's admission and exact replay.
Original native chain-occurrence indices were not serialized; exact ownership
of adjacent same-label native chain objects is therefore not claimed. Seven
observer controls and three ordered-transport controls passed. An additional
private toy audit detected an altered actual heavy coordinate and refused to
call missing-output evidence qualified.

## Coverage limitations and next gate

- Input bytes and [fixture provenance](../../tests/fixtures/README.md) remain
  unchanged; 1DQJ's 84 altloc records were not silently removed from the source.
- Protonated acid O–H is actually observed at pH 3, but the AMBER request used
  forcefield-free H placement. Native CHARMM pH-3 preservation is separate.
- GLYCAM and water O–H preservation are measured. Late RDKit/OpenBabel ligand H
  is not direct matched-postflight ligand coverage. Final PDB-inferred bond
  graphs may also omit nonstandard ligand/sugar bonds.
- The CHARMM-labelled 3RY6 case is **not native CHARMM glycan qualification**.
- One seed label per route on six sources cannot estimate failure rates or
  prove native convergence, protonation/torsion accuracy, lower O–H bounds,
  THR2 HB handedness, hardware reproducibility or production readiness.
- Closing full-cohort acceptance requires a separately scoped fix/qualification
  for the retained GLYCAM/CHARMM fallback duplicate-H route, not an unchanged
  retry. Broad ligand/native-CHARMM glycan coverage also remains absent.

Private coordinates, RNG arrays, XML, checkpoints, failed structures, complete
logs and artifacts are retained below
`tmp/hydrogen-postflight-contract-01/native-run-01/`. Only scalar evidence,
identities and hashes are tracked. The historical lost THR2/SER trajectories
remain lost; this cohort does not replace them.

## Engineering validation

The final full fast suite passed **3,184 tests**, with three optional
Torch/Biotite skips and 33 slow deselections. All ten observer/transport controls,
Ruff, strict-target mypy, generated CLI/GUI specs and agent-map checks passed.
Independent offline audit readback is byte-identical. The original SER test and
production source remain unchanged from `3ffa1f4`; this stage adds research
observation, audit and evidence only. These engineering passes do not override
the retained case-11 failure or the false complete-cohort qualification status.
