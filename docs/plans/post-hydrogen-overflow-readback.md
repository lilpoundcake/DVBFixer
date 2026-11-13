# Prospective CPU/Reference attribution of failed GPU readback

The frozen 1VCU audit retains **26 failed phase-1 callback comparisons** and zero
failed phase-2 comparisons. Early reported gradients reach ~4e12 kJ/mol/nm,
while direct GPU force readback is limited by the fixed-point buffer. Do not
change the original audit, thresholds, native trajectory, guarded refusal or
whole-complex qualification status.

The exact runtime revision's CommonMinimizeKernel.cpp explicitly retries GPU
force overflow on CPU (Reference if CPU unavailable). Test that explanation
read-only, not by optimizing again. Freeze ALL 26 originally failed callback
ordinals, exact System/Integrator/State/path bytes, original audit and this source
reference before queries. Enroll both CPU and Reference for every callback, not
whichever platform looks favorable. Use default CPU properties as in the cited
fallback, the original System, state parameters, integration force groups and
the independently calculated existing constraint penalty. No FF change or force
added to any Context.

Reuse the original prospective energy/gradient/penalty envelopes exactly. These
CPU/Reference queries get no single-precision constraint-arithmetic allowance.
Record all errors, unavailable backends, timeouts, nonexecutions and numerical
residuals, even when CPU agrees. Include whether the reported gradient exceeds
the 2^31 kJ/mol/nm signed fixed-point force capacity (64-bit /2^32 scaling), as
source-consistent supporting evidence, not independent proof of an internal event.

Bound this to two readback Contexts, 26 states per backend, 300 seconds/backend
and 600 seconds aggregate. No native minimization, H generation, random draws,
CPU reruns of native trajectories, extra seeds or adaptive tolerance changes.
Keep all native coordinate/gradient arrays private. Public output is scalar/hash
only. A favorable secondary attribution does NOT make the original direct-GPU
audit pass, erase the astronomical starting-energy warning, certify constrained
stationarity, or qualify the blocked guarded path.

Documentation scope: this amendment, frozen JSON, complete scalar attribution,
then result/handoff update. Pending links: none.
