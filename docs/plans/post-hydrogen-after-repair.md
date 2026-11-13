# Prospective offline post-repair admission control

The public-API pair frozen at `94b8098` is **negative** for the current ordering:
baseline completed two whole phases, guarded stopped atomically after a returned
native H call with `PostHGuardError`. No guarded whole phase or FF fallback ran.
Raw-return readback locates one out-of-admission supported center per arm:
baseline GLU228/HA (~2.137 Å), guarded PRO379/HA (~1.651 Å). Both raw SER0/HA
states remain spatially incompatible despite genuinely FF-provided H placement.
Baseline also reports astronomical starting energy (~2.345e13 kJ/mol); preserve
this warning, rather than call the full workflow clean from its exit code.

Before any changed candidate is considered, perform **offline-only** admission
on BOTH complete retained raw API returns. Never select only the easier arm,
change distance/sign thresholds, weaken old tests, repeat H sampling or use the
baseline trajectory to relabel the guarded run successful. This is a new timing
control: ordinary production H repair first, unchanged selective guard second.
The failed before-repair candidate and all its evidence remain immutable.

Freeze the two raw journals, original plan/index, current production geometry
source and helper code before applying the repair. Reconstruct each exact atom/
graph/bond-metadata table in a new private directory. Source scope is single-model
1VCU with no ambiguous component occurrences. Derive new-H ownership from each
journal's own before/raw identity difference; source observed H remains zero.
Run `repair_misplaced_hydrogens` exactly once on a reconstructed copy, then the
unchanged `propose` exactly once if ordinary repair returns. No Context, force
query, minimization, native H call or new RNG draw. Preserve repair errors and
guard refusals. Check both Python and NumPy RNG neutrality, exact heavy positions,
all per-atom ordinary-repair changes and the guard's separate new-H-only mask.

Report every out-of-admission supported center before/after ordinary repair,
SER0 HA heavy/four-ligand volumes and bond length, ordinary-repair counts and
guard proposals/refusals. Label these diagnostics as a prospective follow-up to
an observed failure, not an original predeclared acceptance endpoint.

For baseline only, compare shadow ordinary-repair coordinates with its captured
first whole-System start using exact canonical variant identities/graph and the
already fixed float32 transport bound. Preserve any mismatch and any heavy motion
as blockers. A successful match is an offline handoff check, not a new native run.
Raw API completion and caller-adapter commit must be separate: the failed
guarded wrapper had a real native return but committed nothing to its caller.

Budget: exactly two shadow ordinary-repair attempts and at most two strict guard
evaluations, 120 seconds total, zero native optimization/context/H-sampling calls.
Write all private raw/repaired/proposed coordinates once and scalar/public
outcomes once. No source parameter or threshold changes.

Do not start another native continuation unless **both** retained states admit
under this unchanged rule and ordinary repair preserves heavy atoms. Even then,
require a separately frozen matching whole-System protocol and retain the failed
public-API experiment as failed. If either refuses, stop the native gate and
report the remaining defect, not a pass over an admitted subset. This control
does not qualify a live after-repair integration or any public/default capability.

Documentation scope: this explicit follow-up contract, frozen JSON, complete
scalar outcomes, then update the result note/handoff. Pending links: none.
