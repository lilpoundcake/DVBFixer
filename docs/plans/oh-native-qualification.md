# Prospective real-input O–H postflight qualification

This extends the [implementation qualification](hydrogen-postflight-contract.md)
on `fix/hydrogen-postflight-contract`. The user approved real-input work and
autonomous continuation, **not a merge into main**. Production physics, existing
tests, parameterization and the approved coarse O–H contract stay unchanged.

## Cohort and budget (frozen before native execution)

Each source runs once for each requested FF, with H seed label **1**:

| Tracked fixture below `tests/fixtures/` | FF requests | pH |
|---|---|---|
| `protein_context/1UBQ.pdb` | amber, charmm | 7 |
| `regressions/1DQJ.pdb` | amber, charmm | 7 |
| `c_glh/8cde_t_u.pdb` | amber, charmm | 3 |
| `3ry6.pdb` | amber+glycam, charmm | 7 |
| `1VCU.pdb` | amber, charmm | 7 |
| `hinge_CH3_glycosylated.pdb` | amber+glycam, charmm | 7 |

Twelve preparations, six source structures. Run order is table order, first FF
then second. Per-case limit **1,200 s**, aggregate **7,200 s**, at most **eight
H calls and eight native H minimizer calls per case**. Keep every result,
upstream error, refusal, timeout and budget-not-executed case; no retries,
replacement cases, additional seeds or outcome-dependent enrollment. Killing a
timed-out process group must also terminate its external-tool descendants.

Command: `prepare INPUT -o OUTPUT --ff REQUEST --ph PH --keep-water
--atom-naming standard`. Other legacy defaults, including PROPKA/Reduce, remain
unchanged. Each worker receives a private copy; original fixture bytes and
manifest checksums must remain unchanged. See
[fixture provenance](../../tests/fixtures/README.md).

## Observation and matched comparison

`scripts/qualify_oh_native.py` freezes source/runtime/fixture/XML digests before
execution. It captures native H boundaries without a reporter, energy/force
queries or changed minimizer parameters; records actual FF XML provenance and
forcefield-free fallbacks; and retains native XML/checkpoints, RNG states and
raw/pre/post-repair coordinates privately inside the project. Per-call seed
initialization/restoration is explicit research orchestration, not production
determinism. Cross-platform native runs are not bit-identical by assumption.

`scripts/audit_oh_native.py` feeds archived pre-fix and current postflight the
**same exact saved nm floats**. Topology reconstruction must not substitute an
Å-round-tripped coordinate array. Require exact replay of the actual candidate
positions, return/error and repair count. Independently measure O–H <1.2 Å and
same-residue nonparent-heavy distance >1.0 Å; check heavy/valid-H/identity/bond
preservation. Existing detector findings define selected H for the broader legacy
repair contract; independent O–H measurements do not depend on those findings.

Load `.dat` only via `DatRecord`. Inspect final PDB-readback O–H as a separate
serialization check, while reporting that inferred output topology need not
recover every ligand/sugar bond. Short O–H (<0.8 Å) and cross-residue heavy
contacts (≤1.0 Å) are descriptive warnings, not newly enforced patch thresholds.

Any new refusal, unmatched replay, unexplained valid-H/heavy/identity motion or
remaining contract violation blocks qualification. Upstream execution failures
and absent scientific coverage also remain visible; an incomplete cohort is
not a broad qualification pass. Preserve all artifacts before any amendment.

## Coverage and limitations

- Record actual XMLs, native platforms and measured hydroxyl classes. Requested
  FF labels alone are not evidence of execution under that FF.
- CHARMM requests on PDB-standard glycans use the existing Amber/GLYCAM bridge
  and output renaming: **not native CHARMM glycan qualification**.
- Unknown 1VCU ligands can trigger forcefield-free Modeller H placement. Later
  RDKit/OpenBabel ligand H is outside the matched postflight; do not count it as
  direct ligand-OH repair coverage.
- Acidic pH intends to exercise protonated acid sidechains, but actual variants
  and bonded O–H must be measured rather than presumed. The original 8CDE input
  contains no explicit GLH.
- 1DQJ has 84 altloc records. Preserve them in the input; ordinary pipeline
  alternate-state selection is not a pristine-deposition guarantee.
- Six sources and one draw label per FF route do not establish error rates,
  hydroxyl accuracy, native convergence or THR2 HB stereochemistry. The separate
  research-only empirical-table licensing/production gate remains closed.

## Engineering checks and retention

Seven toy observer controls cover exact nm transport, matched regression behavior,
normal execution, template fallback, typed-error propagation, both call-budget
aborts, no-reporter boundary capture and RNG/hook restoration (some checks share a test). These are not native
cohort runs. No real input is sampled before the protocol/source freeze commit.

Retained evidence lives below the already ignored
`tmp/hydrogen-postflight-contract-01/`. Tracked reports contain only hashes,
identities, bounded scalar measurements, statuses and coverage summaries, not
coordinates, RNG arrays, native states, checkpoints or empirical tables.
