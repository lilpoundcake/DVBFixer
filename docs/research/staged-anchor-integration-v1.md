# Staged anchor refinement: integration verification

Date: 2026-10-27. Branch: `research/9cqd-anchor-relaxation`.
Status: experimental software integration, not broader scientific calibration.

## Delivered contract

`geometry-gated-anchors-v1` explicitly binds new requests to geometry-first,
generated-only stage 1 and default-enabled raw-restarted stage 2. Both user
opt-outs are request-bound. The limits are per observed anchor, without fitting:
backbone N/CA/C/O RMSD ≤0.5 Å, Cα ≤0.5 Å, maximum heavy displacement ≤1.25 Å.
B factors adjust restraints, not acceptance limits. Missing/untrusted method
provenance uses uniform fallback; observed masks never become generated masks.

The independent validator reads verified raw/stage-1/selected artifacts and
recomputes eligibility, expanded geometry, unselected-record preservation and
movement. The old strict result remains available in metrics. The research
auditor still does not independently grant publication authority; a versioned
public adapter validates request-bound evidence before granting its narrow
observed-motion exception. Only typed numerical stage-1 failure can replace a
stage-1 PDB; arbitrary engine/configuration errors cannot authorize release.

Publication retains the diagnostic report (effective coefficients and B fallback
reasons), relocates and re-hashes evidence and traces. Per-seed subrequests
must match the exact parent-derived request; publication binds the verified
evidence back to the full parent request. Resume copies all evidence into a
private workspace and recomputes validation. Digest tampering, unsafe paths,
missing/misbound evidence and source/artifact path collisions fail closed.

Model, ZBS, Homology and the GUI expose the two switches. GUI outputs and Homology
show explicit Python anchor notices in visible warning banners, not only logs;
clean runs show no empty banner. Node verifies all manifest-declared artifacts
as transport custody without duplicating scientific acceptance. Existing frozen
requests keep their historical semantics and byte-for-byte serialization.

## Native 9CQD replay

Both archived raw candidates were replayed with uniform and experimental-B
policies through the real CPU refiner, public validation, publication and resume:
four runs total. No solvent and no new neural sampling.
The replay script explicitly marks sampler evidence incomplete and uses an
archive-replay profile. Sampler-only context plans are preserved in the copied
archived request, not asserted as newly exercised sampling evidence. Source,
generated identities, masks and seed 7 remain unchanged. This is a refinement
and custody test, not new neural-profile/hardware qualification.

| Archived arm / restraint policy | Released anchor | Backbone RMSD Å | Cα Å | Max heavy Å |
|---|---|---:|---:|---:|
| Reinjection Off / experimental B | W:29 | 0.481456 | 0.397784 | 0.915494 |
| Reinjection Off / uniform | W:29 | 0.421562 | 0.300450 | 0.816287 |
| Reinjection On / uniform | W:25 | 0.195467 | 0.090846 | 0.532246 |
| Reinjection On / uniform | W:29 | 0.428733 | 0.244641 | 0.937126 |
| Reinjection On / experimental B | W:25 | 0.212060 | 0.099303 | 0.604649 |
| Reinjection On / experimental B | W:29 | 0.460908 | 0.339252 | 0.972898 |

All choose generated-and-anchors, pass all public gates and resume validation,
and preserve all unselected observed records. The B-weighted W29 coefficient is
k=8311.4321; uniform uses k=10000 kJ mol⁻¹ nm⁻². Slight differences from earlier research pilot coordinates
are not treated as regressions: seeded minimization does not promise identical
coordinates across numeric backends/builds. No reference-loop RMSD is used.
The earlier v4 On/B generated-only assessment authorized W29 only, unlike the
v5 run's W25/W29 selection. These independently executed minimizations are not
a fixed-stage-1 causal B-effect comparison. No threshold or selection was tuned
to make the runs pass.

Optional local evidence (ignored, not CI inputs):

| Directory below `test/errors_diffusion/` | `report.json` SHA-256 |
|---|---|
| `9cqd-production-staged-off-b-v5/` | 264ab7fc21ac6ebfa0f2e6d5606183bcc366eca2810a725f7c827768f2b9f1a8 |
| `9cqd-production-staged-on-uniform-v5/` | fcdca68ef3417977cc9df68e1c80b0b0040900c3a1fc97d1e9ba57bb0becd8a7 |
| `9cqd-production-staged-off-uniform-v5/` | 1fff35c2f481b10e9322f355298c02da73dbf90b614c7927f72214a76f946782 |
| `9cqd-production-staged-on-b-v5/` | 2ed1286da32ea27ac1df3e4308df93d6105aa05f832fd3358cbaccae05c3195e |

Earlier v1/v2 replay setup failures are retained; they did not perform any
minimization. Earlier v3/v4 native successes are also retained. Every v5 run
captures 25 launch-time Python source files and their digests before refinement;
bundles retain the digest-bound effective-restraint diagnostic report.
The replay command is [scripts/replay_staged_refinement.py](../../scripts/replay_staged_refinement.py).

## Regression coverage and limitations

[Integration tests](../../tests/test_diffusion_staged_integration.py) exercise
raw no-op, generated-only repair, observed-anchor repair, typed numerical
recovery, opt-out rejection, public validation, bundle relocation and resume,
tampering, CLI defaults, and per-seed custody. Pure codec tests cover bounded
input, duplicate/unknown keys, versions and contradictory evidence.

The earlier [eight-case corrected replay](anchor-repair-corrected-cohort-v1.md)
provides five clean raw controls and three generated-only repairs, with zero
observed motion. It does not calibrate the movement limits or establish B
weighting superiority. A broader independent real-gap motion cohort and target
GPU acceptance remain scientific/hardware follow-up, not claims of this change.

Verification: 1587 fast tests pass; two optional-dependency tests skip and 31 slow
tests are deselected. The no-OpenMM/PDBFixer lightweight lane passes 430 tests
with 61 expected skips. GUI: 403 tests pass across 57 files; TypeScript typecheck
passes. Ruff, strict mypy targets, user/agent documentation and generated
CLI/GUI consistency checks pass. The complete external-tool/slow suite and a
new neural GPU sampling run were not executed for this refinement integration.
