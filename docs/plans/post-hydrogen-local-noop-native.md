# Retained-state shared-physics downstream control

## Scope and prospective enrollment

Offline local-policy plan `666c745` admitted BOTH retained post-repair states
with independent audits. This separately frozen experiment tests downstream
behavior only. It does not qualify a live public-API hook, resample H, repeat
ordinary repair, or retroactively pass the old strict gate. Public capability,
production, main and numerical-gradient certification remain disabled.

Use the EXACT captured baseline 1VCU whole-System (including ligands, 213 water
constraints, native positional restraints and chirality wall), Integrator,
state parameters and captured platform properties. This is explicitly a
**shared-physics counterfactual**, NOT an independently constructed guarded FF
or a historical replay of its nonexistent whole-System. Before execution bind
the complete atom identities/elements/canonical residue names, graph and bond
metadata to the baseline topology; compare both raw API selected-variant lists.
Validate all restraint anchors reference heavy atoms, use identical phase-1
System bytes, then the exact captured same-anchor divide-ten phase-2 System.
Do not recenter anchors, change constraints, drop DAN/EPE or select a new platform.
Both repaired states' heavy coordinates must agree with the captured System
start within the previously fixed float32 transport bound. Native transport is
reported separately from exact raw-to-proposal ownership preservation.

Exactly THREE new paths, each two native whole phases, without extra H-only
relaxation, reflecting or post-hoc projection:

1. baseline retained state + its audited local proposal;
2. guarded retained state unchanged after its previously saved ordinary repair;
3. guarded retained state + its audited local proposal.

The original baseline unchanged trajectory is reused as an existing control,
NOT rerun until green. Guarded unchanged was never run previously; it is needed
to distinguish a policy effect from the two raw native H placements' differences.
Each new path uses its own phase-1 return as phase-2 input. Both H-only source
states stay enrolled, even if one fails. Six minimizer calls maximum, native
tolerance10/maxIterations1000 and 4096 nonterminating callbacks/phase. Set no
additional minimizer options or forces. At most 1800 s/path, 5400 s total;
subprocess timeouts kill the complete process group, retain partials and account
all paths. No seed replacements, platform substitutions, budget extensions,
additional phases or native retries. Requests are inherited from the original
`minimize --no-solvent --keep-water` run; add no solvent/ions here.

## Audits and acceptance boundaries

Retain exact requested/applied positions, Systems/Integrator/States/checkpoints,
all callbacks and failures privately in a new ignored project directory. On a
global D-Cα outcome mark failure; do not perform unplanned emergency repair or
publish a structure. These are private downstream experiments, not full CLI
outputs. Verify own-return continuity, native rounding, source-bound starting
coordinates and phase System byte identity independently.

Use the unchanged original constraint-aware `audit_whole` and its frozen
component-wise energy/gradient/penalty envelopes for ALL callbacks. No new
overflow retuning or alternative CPU attribution is enrolled. If physical
attestation fails, retain that failure independently of any positive geometric
endpoint. Constrained stationarity is not certified. Independently measure
SER0/HA and HAs283/286 compatibility and normalized heavy/four-ligand volumes,
bond lengths and PRO379/HA at every initial/callback/return node. Report ALL
sign conflicts and degenerate nodes (same 1e-4 normalized threshold), not just
selected successful final states. Unsupported centers remain outside the guard's
claim. A final compatible target alone is not proof of correct native gradients,
global tetrahedral regularity, CIP accuracy or public readiness.

Compare guarded raw control to guarded projection under identical captured
physics; preserve each outcome. Full local geometric success requires both
proposal paths to return both phases, pass transport/System/ownership audits,
have no negative Cα nodes, and preserve compatible nondegenerate signs for all
three moved H at every recorded node. Report the unmodified control regardless
of whether it improves. Native numerical qualification separately requires all
original physical envelopes; do not conflate these statuses.

Documentation: this contract and frozen plan before execution; all-path scalar
audit and native validation after; update the local-policy note/handoff last.
No future links; pending links none. Coordinates/RNG/Systems/checkpoints remain
private, only hashes/scalars/identities in public records.
