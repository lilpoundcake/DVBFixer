# Isolated synthetic OpenMM constraint-gradient research

## Scope and first prospective series

User instructed private continuation after `330b40a` and the explicit isolated-
reproducer proposal. Work is on `research/openmm-cpu-sign-reproducer`. No installed
OpenMM patch, new environment/package installation, production/main change,
external issue or public capability. Independent child processes use exclusive
project-owned output/temp/cache directories. Native-library directory inspection
was permission-denied and not bypassed; runtime metadata is not binary-hash
attestation. A preliminary unsupported Platform.getLoadedPlugins introspection
was an engineering failure, not a scientific case.

The [first contract](../plans/openmm-cpu-sign-reproducer.md) and
[12-case plan](openmm-cpu-sign-reproducer-plan.json) were frozen at **`287eacb`**
before native execution. Three positive-mass particles, one harmonic bond,
independent quadratic stress, optional pair constraint; no biological structures
or force-field templates. OpenCL-single/CPU/Reference × constraints on/off ×
low/high stress all ran once, with public nonterminating reporters. Analytic
energy/gradients and termwise finite differences are independent of OpenMM.

## Complete first result: desired trigger not reproduced

[All 12 outcomes](openmm-cpu-sign-reproducer-evidence.json) have worker exit0 and
returned native calls, **162 callbacks** total. Finite differences pass all162.
The predeclared primary source-sign signal is **false**; all-controls gate is
also **false**. These are retained negative results, not discarded tuning runs.

- All eight CPU/Reference controls pass physical gradient/energy envelopes.
- OpenCL free-low passes. OpenCL constrained-low passes its gradient envelope
  with the predeclared constraint-arithmetic allowance, but one of43 callbacks
  fails physical energy: error~0.00148483 kJ/mol. Keep this failure, not a wider
  absolute bound. The allowance is not strict gradient equality.
- Both OpenCL high-stress cases produce one callback, with no reported gradient
  above the fixed-point capacity criterion. Their stress-particle derivatives
  are wrong even WITHOUT constraints: an analytic x gradient~6.70954e13 is
  reported as~-6.16563e8. None of the three declared gradient hypotheses matches.
  This is consistent with force-buffer aliasing at excessive forces, not an
  observed CPU-sign reproduction. No internal buffer/fallback event is directly
  witnessed.
- Positive penalty energy/gradient finite differences remain correct. Neither
  source-sign matching in old 1VCU readbacks nor new synthetic native completion
  rehabilitates physical qualification or HA286 geometry.

The huge-stress trigger is therefore not a minimal sign reproducer on this
runtime. Do not rerun these unchanged cases until green, alter their recorded
classification, or modify the old code/tolerances after results. A separately
declared representable-force threshold control is needed to distinguish the
source's overflow detector from out-of-range force aliasing. No dependency patch
or corrected native gradient has been tested.
