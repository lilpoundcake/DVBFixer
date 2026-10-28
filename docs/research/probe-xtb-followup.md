# Probe / xTB: bounded follow-up and decisions

Subsequent implementation: [ADR 0012](../adr/0012-remove-xtb-refinement-and-repair-probe.md)
authorizes removal of production xTB refinement and repairs the optional Probe
adapter. Statements below about unchanged adapters/known failures describe the
historical study at commit `3b4ac0b`, not the current implementation.

Continuation of the [Linux audit](probe-xtb-linux-audit.md), starting at commit
`6d84ca2`, under the [fixed decision protocol](../plans/probe-xtb-followup.md).

**Completed: 60 Probe executions, eight failure controls, 31 xTB optimization
jobs and the API/native controls. No optimization arm passes all declared
comparative screens.** This closes the bounded investigation of the current
integration and proposed context/positional-restraint variants.

## Decisions supported by the combined studies

| Question | Decision |
|---|---|
| Remove Probe because Python already detects clashes? | **No.** Its additional intramolecular coverage is demonstrated. Retain an optional diagnostic backend. |
| Trust the current production Probe adapter? | **No.** The research replacement is feasible, but production invocation, parsing, formatting and engine-status handling still need implementation. |
| Make Probe mandatory or rank detectors by finding counts? | **Not supported.** Real-input precision/recall has not been adjudicated; radii and exclusions differ. |
| Is upgrading xTB enough to make the current refinement safe? | **No.** #1429 fixes tested native constraints, not missing environmental interactions or adapter validation. |
| Can another optimizer rescue isolated GFN-FF polishing? | **Not in this study.** A converged matched L-BFGS-B ligand calculation still loses pocket contacts and increases clashes. |
| Is GFN2 an automatic replacement? | **No.** The explicit charge/spin/solvent arms in the Linux audit also fail interface-retention screens. |
| Is there a qualified full-context production profile now? | **No.** This is an experimental model with runtime, charge/perception and reciprocal-environment limitations; the full three-input cohort is not qualified. |
| Does the long full-context result establish an advantage over simple restraints? | **No for the tested hinge.** A converged restraint-only ablation gives identical clash pairs, retained contacts and very similar coordinates, much more cheaply. |

These conclusions concern structure preparation and complex polishing. They do
not assert that xTB is generally unsuitable for isolated molecules, nor that a
different future context model cannot work. No production/default/dependency
changes were made by these experiments.

## Probe: integration is feasible, with additional identity rules

Research implementation: [`audit_probe_adapter.py`](../../scripts/audit_probe_adapter.py).
It requires the tested Probe 2.18 version, uses `-self ALL`, excludes waters to
match the Python default, parses the correct minimum-gap field, validates every
record and resolves both atom identities against the exact input. One physical
overlap pair is aggregated while all its contact classes are preserved.

**60/60 valid-input executions complete**, covering 15 inputs, original/normalized
atom-field formatting and electron/nuclear radius arms. Inputs include ASN, the
four fresh real baselines, the established covalent/intramolecular/H-bond controls,
and additional case-sensitive chain, insertion-code, blank-chain, negative-number
and lowercase-atom controls.

**Eight failure controls reject safely and explicitly fall back to Python:**
missing executable, timeout, nonzero exit, malformed output, nonfinite gap,
unsupported version, mixed valid/malformed output and ambiguous atom-case aliases.
Every fallback retains requested/actual engine, failure status and nine findings
on the established clashing-methane control. Successful empty results also carry
explicit engine/status metadata.

The first strict prototype rejected GLYCAM `0fA` output rendered as `0FA`.
Upstream `readPDBrecs.c:58–80` explicitly uppercases atom/residue names, while
preserving chains, insertion codes and altloc. The final prototype uses an
injective emitted-name-to-original-identity map. It preserves the original
`0fA` and lowercase atom spellings in identity resolution, preserves `H` versus
`h`, and rejects input atom names which collapse onto one emitted spelling.
This is a parser identity rule, not a conversion of GLYCAM chemistry.

### Formatting alone can dominate the result

The misaligned hinge control changes only columns 13–16. Identity, element
columns and coordinates are unchanged and checked before execution.

| Radius arm | Misaligned overlap pairs | After normalization | Matches baseline pairs and metrics |
|---|---:|---:|---|
| Electron | **3465** | **55** | Yes |
| Nuclear | **3429** | **28** | Yes |

Atom-field normalization is therefore a scientific input-interpretation
requirement for this backend, not cosmetic cleanup. Choose/report radii according
to hydrogen-coordinate provenance; a lower finding count does not select the
correct radius convention.

Probe still detects octane C1–C5 at 2.62675 Å / four graph bonds, missed by the
Python same-residue exclusion. The linear and 90° H-bond controls both retain
`hb`; some pairs also carry `so`. This does not qualify angular H-bond assessment
or establish that every reported overlap is an independently validated clash.

## xTB: separate optimizer, environment and numerical runtime

Research implementation: [`audit_xtb_context.py`](../../scripts/audit_xtb_context.py).
It calls the fixed upstream xTB C API for GFN-FF energy/gradients and uses SciPy
L-BFGS-B over only the selected mobile coordinates. Fixed atoms are never
optimizer variables, rather than being restored after an unconstrained move.
A regression test uses a backend which actively pulls protein atoms away and
checks their exact immobility during every evaluation.

The three inputs are the existing hinge, 7X35/PLM and 1VCU/DAN+EPE baselines.
Each matrix contains isolated/k=0, full-context/k=0 and full-context/k=5 arms.
The full model retains every input atom and cuts no peptide/covalent bonds.
The production extractor selects mobile heterogens plus its mobile protein
anchor-sidechain atoms; all remaining atoms stay fixed.

The restrained objective is
`E_GFNFF + 0.5 × 5 × Σ|r_heavy − r_reference|²` in kcal/mol, with consistent
conversion to Hartree and Angstrom derivatives. Parameters were fixed before
reviewing final geometries; no restraint-strength search was performed.

### Instrumentation controls passed

- Ethanol API/native energy difference is about **1.6e-13 Hartree**.
- Six heavy-atom/H gradient finite-difference checks have maximum absolute error
  **1.47e-10 Hartree/Å**. A separate test checks harmonic-restraint units/masks.
- A native full-hinge single point matches the API initial energy to about
  **4.4e-11 Hartree**: approximately −703.6696820356 Hartree.
- Both separately built native Netlib and OpenBLAS configurations pass **30/30**
  upstream unit suites; disabled optional tblite/CPCM-X backends remain unqualified.

### Python library loading changes the numerical provenance

The initial API arm loads the Netlib-linked xTB library into the scientific
Python process, which has already loaded MKL-backed `libcblas.so.3` through
NumPy. `dladdr` on the xTB handle's BLAS/LAPACK symbols confirms that resolution.
It must be labelled **host-BLAS API reference**, not a clean Netlib API run.
This does not change the earlier standalone Netlib binary/unit provenance.

The separately built OpenBLAS arm resolves those symbols directly to OpenBLAS
0.3.34. Both NumPy/SciPy and vendor libraries can appear in process mappings;
standalone `ldd` alone does not establish the Python process's effective linkage.
The evidence records symbol resolution, loaded-library paths and digests.

The three bounded matrices are host-BLAS/one thread, explicit OpenBLAS/one thread
and explicit OpenBLAS/eight threads. Thread/backend sensitivity is labelled;
bit-identical final coordinates across numerical configurations are not assumed.

### Isolated optimization fails even with a different optimizer

| Input | Isolated external L-BFGS-B outcome, one-thread runs |
|---|---|
| Hinge | 500-iteration limit; maximum projected gradient ≈0.00402 Hartree/Å; diagnostic clashes **33 → 58** |
| PLM | **Converges** in 414 iterations, gradient ≈8.85e-5 Hartree/Å; clashes **57 → 66**, pocket-retention gate fails |
| DAN/EPE | 500-iteration limit, gradient ≈6.50e-4 Hartree/Å; diagnostic clashes **89 → 114**, angle-change and pocket-retention gates fail |

The eight-thread isolated PLM arm converges in 406 iterations and also ends at
66 clashes. No isolated arm passes all comparative screens. Unconverged geometry
is diagnostic only. Native ANC and external L-BFGS-B reach different basins;
neither unlike optimizer energies nor unlike total-system energies are ranked.

### Full-context outcome

All **18 full-context arms** in the three main matrices exceed the 600-second
worker limit. They are completed bounded failures, not accepted structures or
zero-finding results.

The supplemental hinge/context/k=5 run converges in **306 optimizer iterations**
and **1711.6 seconds total elapsed time** (about 28.5 minutes). Its projected
maximum gradient is **9.58e-5 Hartree/Å**, raw frozen drift is exactly zero, all
28 original protein contacts are retained, and relative bond/angle/stereovolume
screens plus the new-mobile-chemical-error screen pass. However, clashes increase
**33 → 35**, so the declared clash gate fails. Convergence and contact retention
are not sufficient to call this an improvement over the OpenMM baseline.

### The restraint-only ablation changes the interpretation

To avoid attributing the combined result to protein context, the same k=5
restraint is applied to the isolated subsystem on all three inputs. All three
ablation optimizations converge with the unchanged gradient threshold.

| Isolated/k=5 arm | Iterations | Observed elapsed time | Clashes | Original contacts retained |
|---|---:|---:|---|---|
| Hinge | 190 | 13.7 s | **33 → 35** | **28/28** |
| PLM | 185 | 6.1 s | **57 → 61** | **26/28** |
| DAN/EPE | 240 | 5.3 s | **89 → 90** | DAN **80/90**, EPE **46/48** and **29/34** |

The restraint-only hinge and the long full-context hinge have **identical
physical clash pairs, severities and retained original contacts**. Their 248
mobile heavy atoms differ by RMS **0.0315 Å**, maximum **0.175 Å**.
Thus this experiment establishes no extra geometric benefit from the expensive
full-context model over the cheap positional restraint. It does not establish
that environment never matters; the restraint may dominate this local basin.

**Restraints reduce the damage but do not produce a passing polish:** every
ablation still increases clashes, and DAN/EPE additionally fails the >10°
relative angle-change screen. Do not adopt k=5 as a validated production fix,
or tune it post hoc until a finding count happens to pass.

The independent native single-point control takes **243.6 seconds** on the
4117-atom full hinge, including setup; its final energy/gradient evaluation is
about 3.8 seconds. Native logs show topology-distance/Floyd–Warshall construction
and six π-occupation fallback warnings. Those warnings require perception review;
they are not, by themselves, proof that the final fallback chemistry is wrong.
Observed timings overlap other jobs and are not calibrated scalability benchmarks.

Limits of the model remain even if an individual geometric screen passes:
the surroundings cannot respond reciprocally, GFN-FF infers its own topology,
the API retains default zero total charge rather than a validated protein
formal charge, and GFN-FF's size-dependent accuracy defaults apply. The model is
not a qualified charge-aware QM/MM calculation or an explicit-charge GFN2 model.

## Engineering disposition and evidence

Production adapters remain unchanged, including seven strict known-failure
contracts. Implementing the tested Probe invocation/parser/identity/fallback rules
is a well-bounded next engineering task. For xTB, output/convergence/identity/raw
drift checks are necessary but cannot manufacture missing environmental physics.
Do not promote an xTB package upgrade or isolated refinement into a recommended
complex-preparation step on these results. Both package and fixed source identify
as 6.7.1, so version strings alone cannot establish #1429 provenance.

Evidence root: `.artifacts/probe-xtb/`. Final matrix/prototype summaries are
generated by [`summarize_probe_xtb_followup.py`](../../scripts/summarize_probe_xtb_followup.py).
The machine-readable [results](probe-xtb-followup-results.json) distinguish
execution completion, native/optimizer status and comparative screening.
Missing geometry from timeouts is `null`, never a zero-error result. Pair-change
summaries distinguish physical atom-pair changes from severity-only changes.

The locally retained, checksummed bundle is
`.artifacts/probe-xtb/probe-xtb-followup-evidence.tar.gz`, with an adjacent
`.sha256` file. It includes frozen inputs, raw logs, reports, script/source
snapshots and build provenance. Large native runtime caches remain locally and
are listed as archive omissions; generated coordinates/reports are research
evidence, not accepted output structures. The bundle is ignored by git and must
be transferred separately from the branch.

The final fast regression run has **1500 passed, 7 skipped, 31 deselected,
7 strict xfailed**, with the previously recorded six SWIG/closed-stream warnings.
The final focused run has **41 passed, 7 strict xfailed**, including two focused
regressions added after the full run (fixed-coordinate projection and pair/severity
distinction). Ruff, `git diff --check` and agent-document checks pass.
