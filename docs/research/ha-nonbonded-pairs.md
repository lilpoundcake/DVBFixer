# The initial inward HA force is dominated by experimental THR24 methyl contacts

Status: **checked static pair-force attribution**, not a force-ablation experiment,
production repair or explanation of the original 1CRN THR2 inversion. Follows the
[native H-minimizer path](ha-minimizer-path.md). No new minimization, H generation,
coordinate perturbation, force modification or native-reference scoring is used.

## Exact complete-cohort decomposition

[Complete evidence](ha-nonbonded-pairs-evidence.json) retains all **24 cells /
1,248 path nodes / 1,190,592 HA-other-particle pair observations**. Every pair is
counted, including excluded exceptions and beyond-cutoff zeros. These remain four
masks × two H seeds × three source-geometry arms of **one protein**, not an
independent population. Each node publishes H/heavy and same/other-residue class
totals and exactly three most-negative outward-axis pairs; complete private terms
are retained, so this inspection view does not remove failures or small terms.

Every original System uses unswitched **CutoffNonPeriodic**, cutoff 1 nm,
reaction-field dielectric 78.3, no particle/exception offsets. The decomposition
uses native Lorentz–Berthelot LJ mixing and reaction-field Coulomb for regular
pairs. Exceptions use their own charge product/LJ parameters and plain Coulomb,
without the regular cutoff. Other force modes/switching/offsets are refused.
The Coulomb MD-unit constant and formulas are calibrated against native Reference
two-particle forces, including exceptions beyond the regular cutoff.

The complete three-dimensional pair-force sum matches the previously measured
**unmodified native NonbondedForce** at every node, with maximum component error
**8.731×10⁻¹⁰ kJ/mol/nm**. A separate streaming auditor independently checks all
1,190,592 private terms, complete identities, class/vector/projection sums and
the declared three-pair inspection rule. Pair-force attribution is therefore not
based on naming the nearest atom or selectively deleting interactions.

The plan was frozen after observing the large nonbonded-group force, but before
identifying individual contributors: `tmp/pdbfixer-research/ha-nonbonded-pairs-plan-01.json`,
SHA256 `7bb42ad94178f7d4817bfeaa5d1c6d1c025871aceb2d2618d2d0b492a9294735`.
Private complete pairs are under `ha-nonbonded-pairs-run-01/` beside that plan;
no force-field parameter tables or complete raw coordinate/gradient/pair archives
are published. These local artifacts do not transfer through Git.

## Initial seed-1 force contributors

Known THR24 / **full-generated-heavy rebuild** / original H seed 1, before the
native H-only minimization. The focal atom is generated **GLN16 HA**.
Outward-axis force is along the same fixed heavy-derived direction as the prior
path audit; a negative value points toward the inward side.

| Partner | Distance | Total outward projection | LJ contribution | Coulomb contribution |
|---|---:|---:|---:|---:|
| THR24 HG23 | **1.356899 Å** | **−43245.578** | **−43207.680** | −37.897 |
| THR24 CG2 | 1.923824 Å | −6494.745 | −6571.439 | +76.694 |
| THR24 HG21 | 1.635608 Å | −3687.309 | −3661.691 | −25.618 |

Force units are kJ/mol/nm. All three are ordinary **inter-residue** pairs, not
exceptions. The total native nonbonded projection is −56105.474: H partners
contribute −49620.742, heavy partners −6484.732; other residues contribute
−56102.196 and the same residue −3.279. The dominant starting impulse is thus
short-range **Lennard–Jones repulsion involving the experimental THR24 methyl
group**, especially HG23, not a same-residue HA naming error or an isolated
heavy-atom angle term.

This cell's THR24 heavy sidechain was experimentally reconstructed before H
placement. Do not call that modified candidate deposited-source geometry or
attribute it to all ordinary production inputs. Locally regular heavy-atom bond
geometry alone did not guarantee a safe neighboring H environment. This is a
negative qualification result for that experimental reconstruction route, not a
new default method or a license to rewrite the neighboring observed residue.

## Force contributors change after the large optimizer displacement

At callback 0, HA has already crossed inward with CA–HA radius 2.506358 Å. Its
largest inward nonbonded contributor is now **its own GLN16 O** at **1.408459 Å**:
outward-axis projection **−51601.067**, comprising LJ −51737.285 and Coulomb
+136.219 kJ/mol/nm. This is an original **exception pair** (the local HA–CA–C–O
path), handled with the original exception parameters rather than regular mixing.
The same-residue class now dominates (−51626.056), whereas other residues sum to
+442.231. The initial THR methyl contact is not falsely described as the sole
force source at every later optimizer state.

## Seed 2 is an important counterexample to a universal outcome

For full-generated THR24 seed 2, the initial three largest inward pairs are also
the target methyl hydrogens:

| Partner | Distance | Outward projection, kJ/mol/nm |
|---|---:|---:|
| HG23 | 1.310981 Å | −68732.623 |
| HG21 | 1.353550 Å | −45085.763 |
| HG22 | 1.422913 Å | −23323.419 |

Its initial nonbonded projection is even more negative (−142933.015), but that
H-stage replay eventually returns outward and never has an all-acute callback.
After its first large displacement, the largest inward contributor is LYS18 NZ
at 1.388442 Å (−54846.394). Thus neither the sign/magnitude of initial force nor
a close methyl-H contact alone predicts the final inward basin. Complete path
and pair evidence, rather than a favorable single seed, remain necessary.

## Interpretation, validation and next boundary

The validated contributors explain **which pairs supply the recorded force**.
They do not prove that deleting a selected pair would prevent the failure, isolate
a unique L-BFGS line-search mechanism, establish convergence or qualify a repair.
The earlier one-HA intervention remains the controlled whole-system result:
one GLN16 failure prevented, total raw/rounded geometry passes only 10/24→11/24.
Original THR2 causality and the remaining target THR/VAL failures stay separate.

Nine formula/native-Reference tests and five independent private-record audit
tests pass; combined angle/H-stage/pair controls pass 73 tests and targeted Ruff
passes. Both complete private pair runs and independent publications regenerate
byte-identically. Full fast suite: **2,863 passed**, three optional Torch/Biotite
skips and 33 slow tests deselected; agent-map and generated CLI/GUI checks pass.
This is not broader scientific/platform qualification. Production and
empirical-table license/approval gates are unchanged.

Next meaningful work needs a separately declared environment-aware initial-H or
candidate-geometry mechanism control, with all matched cases/failures retained.
Do not silently ablate nonbonded pairs, tune exclusions, choose seeds by outcome,
or promote a global HA rewrite into production.

Code: `scripts/audit_ha_nonbonded_pairs.py`,
`scripts/audit_ha_pair_readback.py`.
