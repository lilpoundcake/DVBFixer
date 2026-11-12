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
