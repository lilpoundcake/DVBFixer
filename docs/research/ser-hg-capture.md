# SER-HG: unavailable original artifacts and bounded fresh evidence

## Scope, source and retention

The user authorized access to the failed pytest directory. The read tool still
rejected external-directory access; the user's own copy command then reported
that the source did not exist. Thus the original failed `repaired.pdb` and `.dat`
are **unavailable**, not inspected, recovered or replaced. The retained original
[failure record](thr2-hswap-validation.json), project-local full-suite log and XML
remain authoritative: **3,082 passed / one failed**, OG–HG **1.245813 Å** against
the unchanged **<1.2 Å** assertion. New successful runs cannot turn that run green.

After that loss the user instructed autonomous continuation. This stage uses only
the tracked, synthetic two-SER input literal in `tests/test_prepare_broken_geom.py`.
It is not another protein cohort or a reconstruction of the lost trajectory.
Production code, scientific thresholds, test assertions and safeguards stay
unchanged. Documentation scope: this protocol/results note, a continuation link
from [H-swap](thr2-hswap.md) and the canonical
[handoff](../plans/pdbfixer-reconstruction-handoff.md); validate links and evidence.

## Prospective protocol

The [frozen plan](ser-hg-capture-plan.json) is a byte-identical public metadata
copy of `tmp/pdbfixer-research/ser-hg-capture-plan-01.json`; raw payloads stay private.

Freeze source/runtime/input hashes before scientific execution. Enroll exactly
**six fresh preparation invocations**: one unseeded and H Python seeds **1–5**,
each once. The same `prepare main([input, '-o', output])` call as the tracked test
uses otherwise unchanged defaults and automatic platform routing. Seed and
restore Python RNG only around native `Modeller.addHydrogens` for seeded labels;
the unseeded label preserves natural RNG evolution.

Observe native H minimizer boundaries without adding a reporter, energy query,
force query or changed optimizer arguments. Retain System/Integrator/State XML,
checkpoint, actual platform, topology/bonds and H before/after coordinates. Wrap
the actual production H-repair call to retain its input/output, pure detector
findings and return count, forwarding arguments/results unchanged. Preserve full
RNG state privately. ForceField wrappers retain the native signature/defaults.

All artifacts are written from the start into exclusive directories below
`tmp/pdbfixer-research/ser-hg-capture-run-01/`, not disposable pytest directories.
Retain every outcome, error, timeout and budget nonexecution. Bounds: **eight H
calls and eight native H calls per case**, **600 s per case / 3,600 s total**.
No extra seed, repeat-until-failure, iteration change, repair or acceptance filter.
Every final distance assertion is reported separately from execution success.

## Deterministic detector-contract control

Before examining fresh outcomes, freeze a static local-graph sweep at O–H
distances **0.97, 1.199, 1.2, 1.2458133889150464, 1.551, 1.552, 1.553 and
1.7 Å**. These toy graphs have one H bonded to O, one heavy neighbor and no close
nonparents. Report the existing detector's findings, actual repair count and
post-repair distance, without changing thresholds or running a minimizer.

Current production's distance guard is **0.97×1.6 = 1.552 Å**, whereas the regression
requires **<1.2 Å**. The synthetic 1.245813 Å O–H graph reproducibly fails that
regression predicate yet produces no finding or repair. This establishes a
**contract gap**, not the complete geometry or cause of the lost failure. Other
detector predicates can still repair an atom below that length if it coincides
with a nonparent; do not infer those predicates were absent in the lost output.

## Interpretation boundary

Fresh native captures can show whether an independently generated endpoint lies
in the gap and identify its observed H/repair lineage. They cannot prove the lost
run's stochastic path or qualify a general change in O–H tolerances. No native
sampling is repeated merely to obtain a pass or failure. Production integration,
upstream issue publication and threshold changes are outside this stage.

## Results: no reproduction of the lost sub-threshold endpoint

Protocol checkpoint **`89aaf27`** precedes all scientific execution. All **six
preparations / six native H calls** completed on automatically selected
**OpenCL / Apple M3 Pro / single precision**. No additional seed or preparation
was attempted after seeing the outcomes.

| Label | Native returned OG–HG (Å) | Native HG–OXT (Å) | Production repair count (all H) |
|---|---:|---:|---:|
| unseeded | 1.680236 | 0.002686 | 1 |
| seed 1 | 1.674850 | 0.004394 | 1 |
| seed 2 | 1.747696 | 0.091563 | 6 |
| seed 3 | 1.677872 | 0.000914 | 1 |
| seed 4 | 1.678057 | 0.000371 | 5 |
| seed 5 | 1.669148 | 0.029154 | 8 |

Every native H return has HG near OXT and an over-long bond to its actual OG
parent. The existing production detector catches HG as `parent_distance` in
every case. Repair restores OG–HG to **0.970000 Å**, and every final PDB has
**0.970200 Å** after output rounding. All six final outputs satisfy both original
test assertions. Different repair counts include other H atoms, not six HG
repairs in one residue.

The [complete boundary/contract audit](ser-hg-capture-evidence.json) retains these
negative native intermediates as well as the successful final repairs. It also
records the static contract gap at 1.245813 Å. **None of these six final outcomes
reproduces the logged original endpoint**, and they neither erase nor explain
the entire lost failed trajectory. The original full fast suite remains failed.

## Saved-force evidence: current collisions do not require an index swap

The separate [saved-force audit](ser-hg-force-evidence.json) examines all six
native Systems and both boundaries in both captured OpenCL and Reference,
**24 readback states without any new minimization**. All six serialized native
Systems are identical.

For the explicitly named HG/OG/OXT atoms in these captures:

- HG is bonded to **OG**, with a native harmonic equilibrium length **0.96 Å**;
  it is not bonded to OXT. The atom map, topology and force-term indices agree.
- HG is mobile (**1.008 Da**); OG and OXT are pinned (**0 Da**) during H placement.
- HG and OXT are **five topology bonds apart**, with **no nonbonded exception**.
- HG charge is **+0.4474 e** and its LJ epsilon is **zero**. OXT charge is
  **−0.8132 e**. The mixed HG–OXT LJ epsilon is therefore zero.

This pair has attractive Coulomb `1/r` behavior without a repulsive pair LJ core.
The finite cutoff reaction-field correction does not remove that short-range
singularity. With the heavy atoms pinned, approaching OXT can dominate the finite
O–H stretch penalty. This is a documented local mechanism under the captured
H-placement conditions, **not a blanket claim that AMBER is defective** or that
every observed H problem shares this cause.

At the saved returned positions, Reference HG force projected toward OXT ranges
from **6.32×10⁵ to 3.68×10¹⁰ kJ/mol/nm**, dominated by nonbonded force. None of
the six returns satisfies the native gradient criterion in either readback mode.
These are not converged, chemically correct H minima. No callback path or native
error enum was recorded in this protocol, so the exact termination reason and
the presence/absence of terminal leakage are **not established**.

The older source/test commentary attributing such a collision to HG inheriting
OXT's index is not supported by the **present six captured Systems**. It must not
be treated as proven causality for these captures. Nor can the missing original
System be reconstructed from that older commentary.

Force-group sum discrepancies are retained rather than called exact: maximum
energy discrepancy **0.000771 kJ/mol**, maximum HG-projection discrepancy
**37.55 kJ/mol/nm** at extremely large single-precision force magnitudes. Native
runtime-declared source provenance is not binary attestation.

## Verification and next authorization boundary

The focused research suite passed **30 tests**. Ruff, generated CLI/GUI specs,
agent-documentation and changed-document link checks passed. The six captures and
their input hashes are immutable; independent boundary and force audits were
repeated with byte-identical reports and without new minimization. Original logs
and the unavailable-artifact status remain retained. No full-suite retry was used
to manufacture a green run. The [validation record](ser-hg-validation.json)
keeps focused checks, fresh outcomes and the older full-suite failure separate.

The research stage is complete at its declared budget. A sensible next objective
is to reconcile the production H-detector/post-repair contract with the unchanged
regression predicate, **not to change the native force field or retry seeds**.
That is a behavior-changing implementation scope and needs explicit user
authorization on a separate branch. It is not automatically approved by the
research instruction.

Before any default change, freeze qualification on identical retained starts,
including hydroxyl/phenol/protonated-acid/solvent and ligand cases; preserve
existing identity, insertion-code, chirality and nonprotein routing rules. Test
over-long-but-noncoincident H separately from gross collapse, avoid unconditional
movement of already valid H, retain adverse outcomes, and do not weaken the
existing test or apply a blanket absolute cutoff to every parent element. The
six synthetic preparations alone do not qualify a general replacement threshold.

**Subsequent approved implementation:** the user authorized a separate fix and
qualification branch. The [O–H implementation record](../plans/hydrogen-postflight-contract.md)
retains the first candidate's failed validation and the amended candidate's
distance/contact qualification. Its full fast suite passed 3,174 tests. This
follow-up does not alter the research results or recover the lost original run;
the implementation is not merged or released by this note.
