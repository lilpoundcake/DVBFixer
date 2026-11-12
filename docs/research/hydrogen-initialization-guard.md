# Private selective H-initialization guard

## Static outcome

On `feat/hydrogen-initialization-guard`, the contract was frozen at `1302ba6`
and implementation/cohort at **`1da6783`**, before any scientific proposal.
See the [contract](../plans/hydrogen-initialization-guard.md),
[static plan](hydrogen-initialization-guard-static-plan.json),
[complete outcomes](hydrogen-initialization-guard-static-evidence.json) and
[independent audit](hydrogen-initialization-guard-static-audit.json).

All **15 native initialization states** were admitted and evaluated once:
three original THR2 starts and all twelve fixed real-input cases. The new-H mask
was derived from positive H-only System masses and independently matched to
absence from pre-H observations; heavy and observed-H particles must be pinned.
Occurrence metadata is preserved; ambiguous legacy keys refuse rather than merge.

- Only original **unseeded THR2/HB** moves. Its parent–H radius is preserved;
  all other initialized atoms remain exact.
- The other **14 states are exact no-ops**, including all twelve real-input
  states. No admission/proposal refusals or execution failures occurred.
- Across these correlated states, 34,403 new H were accounted for: **4,880**
  supported single-H carbon centers and **29,523** explicitly unsupported H.
  These counts are not independent biological samples or general-H coverage.
- All 19 historical H invocations remain in the enrollment, including four
  failed FF invocations before forcefield-free fallback. Actual successful
  routes are four forcefield-free, three CHARMM XML, three AMBER/GLYCAM XML,
  two AMBER19 XML and three original THR2 routes with **unrecorded XML provenance**.
  Requested CHARMM on a glycoprotein can actually execute AMBER/GLYCAM or fallback;
  this is not native CHARMM glycan qualification. Late ligand-H placement is
  outside this matched Modeller-initialization coverage.

The independent NumPy audit rebinds proposal inputs to original State coordinates,
atom/graph identity and masks, then checks enrollment, selectivity, formula,
radius, sign and exact preservation of all unmoved atoms. It shares the source
admission adapter, not the proposal function. Its second offline readback is
byte-identical. No Context, energy/force query, minimization, native H placement
or new RNG stream was used in this static stage.

## Engineering and interpretation

The prototype lives only in `scripts/hydrogen_initialization_guard.py` and is
never imported/registered by production. Twenty-two engineering tests passed,
covering observation preservation, exact no-ops, all 24 proper rotation controls,
permutation invariance, chain case/insertion/occurrence separation, atomic
refusal, invalid data, gross overlap, observed-H mobility and ambiguous metadata.
The shared variant and canonical atom inventories remain the chemistry authority.

Final engineering validation is **3,237 fast tests passed**, three optional
Torch/Biotite skips and 33 slow deselections, including both offline-accounting
regressions. All 24 focused tests, Ruff (production and new research files),
strict-target mypy, generated CLI/GUI specs and agent maps passed. The
[validation/quality record](hydrogen-initialization-guard-validation.json) binds
the evidence, failed parent summary, offline recovery and full-suite digests.
It deliberately keeps local physical readback, orchestration failure, missing
process exit and production readiness as separate statuses.

The rule detects **spatial heavy/H incompatibility**, not CIP, D/L or rotamer
accuracy. It cannot correct wrong heavy stereochemistry. The 0.05 nm overlap
guard prevents gross collisions, not all steric problems. Noisy unsupported H
and compatible H retain their original initialization exactly.

The new selective, radius-preserving rule is not the older unconditional 0.1 nm
HB-only intervention. Static success does not transfer that intervention's
dynamic result. The independent real contexts require no changes, so they test
noninterference, **not protection against natural conflicts elsewhere**. Broader
activation/context coverage remains open; no production integration is approved.

Private arrays remain in
`tmp/pdbfixer-research/hydrogen-initializer-guard-static-01/`; tracked artifacts
contain scalar measurements, identities and hashes only. Original trajectories,
missing artifacts, prior failures and empirical-table licensing gates are intact.

## Next gate

Before claiming dynamic protection even in THR2, separately freeze a native replay
of the **actually changed** proposal. Reuse exact previous unchanged controls and
no-op starts by hash instead of re-running identical cases for favorable outcomes.
Do not infer additional independent scientific trials from such reuse. Native
FFs, restraints, masks, old bridge and stopping rules must stay unchanged.
Production needs a separately approved integration boundary: a global minimizer
monkeypatch is not a supported initializer API.

## Actual-proposal native result

The separate [native protocol](../plans/hydrogen-initializer-guard-native.md) and
[native plan](hydrogen-initialization-guard-native-plan.json) were frozen at
**`a430c86`** before execution. Only the actually changed unseeded proposal ran:
two new paths / six native calls. All six phases returned and saved their States.
The parent runner **failed** while creating its second summary: the helper opens
files exclusively, and `write(evidence.json)` was incorrectly inside the loop.
The frozen research runner is retained with that defect; do not reuse it as a
general orchestrator or silently edit its source binding.

The explicit [offline recovery amendment](../plans/hydrogen-guard-offline-summary-recovery.md)
was committed at `8d31530`. Recovery created a new audit view and index, verified
all saved return flags/counts/States, and changed no original files. See the
[recovery accounting](hydrogen-initialization-guard-summary-recovery.json) and
[independent native audit](hydrogen-initialization-guard-native-audit.json).
There was **no native retry**. The second worker's process exit code is unknown;
its status is `phases-returned-exit-unrecorded`, **not completed/exit zero**.
Keep that distinction from its independently verified returned native States.

| New trajectory | H-only returned Cβ | Whole phase 1 Cβ | Whole phase 2 Cβ |
|---|---|---|---|
| captured OpenCL | R | R | R |
| Reference (process exit unrecorded) | R | R | R |

Both reused unchanged unseeded controls give final S; the four exact no-op seed
1/2 proposal trajectories reuse their earlier unchanged R controls by hash. These
reused paths are not new runs, seeds or independent-context evidence.

All **356 new callbacks** match independent same-platform physical mobile
gradients and objectives exactly (zero maximum error). Every phase System is
byte-identical to its original; initial unselected coordinates, H-only pinned
heavy atoms and own-return phase continuity remain exact. The original bridge
makes **zero repairs** on both H returns. No negative Cα node or returned/last-
callback position mismatch occurs. Maximum new OpenCL transport component error
is **3.83935e-8 nm**, inside the predeclared rounding bound.

Neither H-only return meets the native gradient criterion at the unchanged
50-iteration limit. All four whole-system returns meet their native criterion;
this is not global convergence or comprehensive geometry qualification.
Independent native audit readback repeats byte-identically. The local result is
favorable for this actual selective candidate, while the failed parent runner,
unrecorded process exit and nonconverged H-only returns remain explicit.

## Production integration decision remains separate

The installed OpenMM reports `8.5.2.dev-36a30cb`, revision
`36a30cbca54e727b216b606f3c011b67201eb8b4`. Its public signature is
`addHydrogens(forcefield=None, pH=7.0, variants=None, platform=None, residueTemplates={})`.
The [exact upstream revision, lines 753–1107](https://github.com/openmm/openmm/blob/36a30cbca54e727b216b606f3c011b67201eb8b4/wrappers/python/openmm/app/modeller.py#L753-L1107)
constructs coordinates, creates its System/Context and minimizes internally before
return. No public coordinate-initialization callback is exposed there. The API
and upstream source inspection do not independently attest the installed body
or GPU binary. No external installed files were modified.

Keeping this **pre-H-minimization** intervention would need an explicit, maintained
upstream/local hook (or a larger version-bound implementation fork), plus further
qualification. A post-`addHydrogens` guard could avoid that hook, but it is a
**different intervention timing** needing a new contract and native qualification;
the present result does not establish it. Never silently replace this boundary
with a global minimizer monkeypatch. The maintenance/API direction needs user
agreement before public integration work; independent natural activation and
broader scientific qualification still remain open either way.
