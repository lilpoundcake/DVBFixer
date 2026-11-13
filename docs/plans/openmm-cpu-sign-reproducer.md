# Isolated synthetic OpenMM CPU-sign reproducer

## Authorization, isolation and documentation

User instructed continued work after the explicit proposal for an isolated
OpenMM reproducer, following `330b40a`. Work is private on
`research/openmm-cpu-sign-reproducer`; no installed dependency patch, new package
installation, production change, main merge, public issue or release. Earlier
negative structural results and all frozen contracts remain unchanged.

Isolation means independent child processes, exclusive project-owned result/
temporary/cache directories and synthetic Systems. Use the existing scientific
interpreter read-only, without modifying its environment. Metadata reports
8.5.2.dev-36a30cb / revision 36a30cbca54e727b216b606f3c011b67201eb8b4.
No full installed-binary hash attestation is available: external native-library
directory discovery was permission-denied and is NOT bypassed. The preliminary
public introspection attempt also found Platform.getLoadedPlugins unavailable;
neither is a scientific trial. Preserve these limitations, do not relabel runtime
metadata or source inspection as binary provenance.

Documentation order: this contract, frozen plan and synthetic engineering tests,
all-case evidence/audit, result note/validation and canonical handoff. User's
continued-work instruction covers these bounded records. No future links before
creation; pending links none. Validate digests/local links and absence of private
coordinate/RNG/State payloads before publication. Engineering, observed gradient
inconsistency and any claim of internal CPU fallback causality are distinct.

## Fixed synthetic system and factorial controls

Three positive-mass particles (12 Da), no molecule/FF templates, solvent,
nonbonded interactions, virtual sites or biological inputs. Particles 0/1 have
a harmonic bond k=10000 kJ/mol/nm², r_eq=0.12 nm. Their initial separation is
0.100005 nm, centered on the x axis. Particle 2 starts at (1,1,0) nm with
external energy `0.5*k_stress*(x²+3*y²+z²)`. Its force is analytically known,
independent of the bonded pair. Constraint, when enabled: particles0/1 at0.1 nm.
The initial relative constraint residual is 5e-5, inside native working tolerance
1e-4 but outside requested final tolerance1e-5. This may expose a finite penalty
gradient during the optimizer's early steps; verify, never assume it survives
native initial projection. Masses are all positive, avoiding the separate
source mass-index concern.

Exactly 12 cases: constraint on/off × stress k=1000/1e14 × OpenCL-single/CPU/
Reference, one native minimizer call each. OpenCL properties are explicitly
DeviceIndex0/OpenCLPlatformIndex0/Precision single/UseCpuPme false/
DisablePmeStream false, matching the earlier recorded platform; no substitution
if unavailable. CPU uses defaults, Reference uses defaults. This is a new
factorial synthetic study, NOT a re-execution of 1VCU or a native FF change.

Integrator Verlet step0.001 ps, constraint tolerance1e-5; do not integrate.
Minimizer tolerance10 kJ/mol/nm, maxIterations8, at most256 callbacks per case.
Reporter is nonterminating and only copies public arrays/statistics. No global
minimizer hook, private Context introspection or modified native source.
Save raw trajectories, Systems/Integrator/States/checkpoints privately. Preserve
initial projection, errors, no-callback outcomes and returned motion. No random
draws, H generation, protein preparation or extra relaxation.

## Prospective analysis and criteria

Independent analytic bond/stress energies and gradients at every reporter point;
positive constraint-penalty energy/gradient use a separately tested Cartesian
derivative. At every callback compare three fixed hypotheses: physical+penalty,
physical-penalty, float32(physical-penalty) transported back to float64. The last
two describe inspected source behavior, NOT correct objective derivatives.
Apply original absolute/relative bounds: gradient1e-3 +2e-5*max(1,abs(g1),abs(g2));
physical energy1e-3 +2e-6*max(1,abs(E)); penalty energy1e-3 +2e-5*max(1,abs(Ep)).
For OpenCL only, also report the original separate constraint-arithmetic
allowance: 32*eps32*abs(K)*(r+r0), summed at constraint endpoints. Report both
strict and allowance envelopes for every hypothesis; do not choose or widen a
bound after results. No allowance for unconstrained particles or CPU/Reference.

Termwise central finite differences (h=1e-6 nm) validate the analytic physical+
penalty gradient at every callback. Difference each energy term separately
before summing to avoid cancellation of small pair energies by huge independent
stress energy. Fixed finite-difference bound1e-5 +1e-6*max(1,abs(analytic),abs(FD));
this derivative-control bound is separate from unchanged native-gradient bounds.
Retain all failures and scalar per-component witnesses; do not use the wrong
source sign as a substitute physical derivative.

Primary reproducer signal: at least one high-stress constrained OpenCL callback
with reported gradient above2^31 fixed-point capacity and nonzero penalty
gradient rejects the physical derivative even with the original constraint-
arithmetic allowance, while a predeclared source-sign hypothesis agrees. Require
matching physical energy and finite differences. Control interpretation requires
all reported CPU/Reference gradients and low-stress OpenCL gradients to agree
with physical+penalty within their original allowance rules; high-stress
unconstrained OpenCL must also agree. Missing callbacks/inactive penalty are
inconclusive, not passes. Report control failures even if primary signal appears.

This can establish installed runtime public-reporter/objective inconsistency
on a synthetic System. It cannot independently witness a specific internal
fallback event, prove source/binary identity, show global minimizer correctness,
or explain structural HA286 inversion. No corrected dependency, trajectory retry
or public issue follows automatically.

## Budget and failure preservation

Twelve child workers, each at most30 s and total360 s; kill timeout process
groups, preserve partial artifacts and account every case, including spawn/
platform/serialization failures and nonexecutions. No retries, alternative
constants, altered starting points, additional seeds or adaptive budget changes.
Freeze code/tests and retained upstream exact-revision source before execution.
No empirical structures, RNG arrays or fitted tables enter tracked artifacts.
Public output consists only of scalar/hashes/case settings; toy coordinate and
State arrays remain in the ignored project directory.
