# Probe / xTB follow-up: decision protocol

Starting point: Linux audit commit `6d84ca2`; the installed adapters fail seven
contracts and isolated xTB optimization worsens complex geometry. This study
finishes the bounded dependency/integration investigation, not a general force
field qualification.

## Questions and stopping rules

1. Can a version-checked, correctly invoked Probe research adapter resolve
   identities and contacts on the existing controls and four real baselines,
   while rejecting malformed/native-failure cases? Preserve both radius arms,
   case-sensitive chains, insertion codes, blank chains and formatting controls.
   A pass supports repairing the optional adapter; counts do not establish
   comparative sensitivity/specificity or angular H-bond quality.
2. Does including the entire fixed complex prevent the degradation observed
   with isolated heterogens? Use the fixed xTB library's GFN-FF gradients with
   external L-BFGS-B over mobile coordinates only. Match the optimizer between
   isolated and full-context arms. Full context preserves all input atoms and
   cuts no covalent bonds; mobile heterogens and anchor-sidechain atoms are
   selected by the production subsystem extractor. Its frozen anchor mask and
   all other surrounding atoms remain fixed.
3. Does a predeclared 5 kcal/mol/Å² harmonic restraint on mobile heavy atoms
   preserve geometry better? Compare zero-restraint isolated, zero-restraint
   full-context and restrained full-context on hinge, PLM and DAN/EPE baselines.
   No parameter search or outcome-selected input substitution.

First validate API gradient units by finite differences and match an ethanol
single-point energy against the same native binary. Use one numerical thread,
500 optimizer iterations, bounded worker execution and retained native logs.
Report timed-out/unconverged jobs as such. Native ANC/L-ANC allocates dense
full-system Hessians even with many fixed atoms; an external mobile-coordinate
driver avoids conflating that allocation with the environmental question.

Runtime amendment before reviewing full-context geometries: the Netlib hinge
context arms exceed 600 seconds after only a handful of evaluations. Replay
the unchanged nine-arm matrix with a separately built OpenBLAS backend from the
same xTB commit, first at one thread, then a labelled eight-thread sensitivity
arm. Both optimizers/model selections and geometry gates remain unchanged.
Record actual in-process library mappings, not just standalone `ldd`: Python
may have loaded BLAS dependencies before the xTB library. This is a feasibility
sensitivity, not permission to reinterpret timeouts as chemical passes.

The loader control subsequently confirms that the first Netlib-linked API arm
resolves BLAS through the scientific environment's preloaded MKL-backed
`libcblas.so.3`. Label it the **host-BLAS API reference**, not an isolated Netlib
API qualification. The explicit OpenBLAS build resolves BLAS symbols directly
to its OpenBLAS library. Earlier standalone Netlib binary/unit results are
separate and retain their recorded provenance.

One supplemental convergence check replays the hinge context/k=5 arm with
1800 seconds, retaining the same 500 iterations, gradient threshold, inputs,
restraint and geometry gates. Selection is based on its near-threshold residual
gradient in the bounded pilot, not an observed final geometry improvement.
It cannot substitute for failed/unfinished members of the original matrices.
An eight-thread native full-hinge single point independently checks whether the
API overhead, rather than the full model, explains the runtime.

After that long hinge run converges with contacts retained but 35 versus 33
baseline clashes, add the logically necessary restraint-only ablation: replay
isolated/k=5 for all three inputs, with the same source, eight-thread backend,
500 iterations, 600-second limit and unchanged screens. This separates the
combined context+restraint result from a claim that protein context alone caused
the improvement. No new force constant is selected. Stop after this ablation.

Screen candidates for convergence (including projected maximum gradient
≤1e-4 Hartree/Å), finite coordinates, stable atom mapping, zero frozen drift,
no new relative bond changes >0.1 Å, angle changes >10°, tetrahedral sign flips
or new mobile-involving chemical ERROR findings. Require non-increasing Python
clash count and retention of at least 80% of each residue's original 4 Å protein
contacts when any exist. Record pair-level gained/lost clashes as well as counts.
These are predeclared comparative screens, not affinity or absolute chemistry
certification. Frozen protein atoms cannot relax reciprocally; GFN-FF still
infers its own chemistry from coordinates. Preserve the production route's
default charge handling (API zero total charge) rather than guessing a formal
charge from protein sequence; this is an additional limitation of the context
model, especially for protonated proteins. GFN-FF's own size-dependent accuracy
defaults also apply. No explicit-charge GFN2 context qualification is implied.

The study stops after controls and this fixed matrix, preserving failures.
Conclusions must distinguish adapter feasibility, constraint mechanism, missing
environment, optimizer dependence and evidence limits. Production/default/
dependency changes require an implementation decision after the report.
