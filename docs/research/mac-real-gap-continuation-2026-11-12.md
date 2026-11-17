# Natural-gap reconstruction continuation: execution log

Date:2026-11-12. Implements the
[200-entry continuation handoff](../plans/mac-reconstruction-real-gaps.md).
Primary objective: complete authentic missing regions with usable geometry and
observation custody, not recovery of a crystallographic reference pose. No new
gap/tail-length cap and no new public capability.

Current state: all declared restoration calculations are settled; earlier
launch/queue statements below are chronological history. See
[settled consolidation](#settled-seed101-results-and-all-input-consolidation) and
the [current WP tracker](../plans/mac-reconstruction-real-gaps.md#authoritative-current-state--retained-data-consolidation).
No additional two-case restoration runs are queued.

## 2026-11-14 autonomous terminal geometry continuation

Current outcome supersedes earlier launch/status statements, not their evidence.
The approved native oxidized TGF-β3 target remains 21 pairs; original 18 explicit
annotations, frozen request and failed native graph-transport result are unchanged.
No new structure, neural sample, resource-limit increase or observed-coordinate
movement was used in this continuation. Evidence lives in
`tmp/mac-real-gap-continuation-2026-11-14/`.

| Retained attempt | Outcome / limitation |
|---|---|
| `1KTZ-terminal-torsion-study-attempt01` | 24 graph-bounded proper rotations; selected ordinal 3 / +45°, amide restored to trans; residual packing rejected |
| `1KTZ-inward-wall-refinement-attempt01` | All four returned native stages within 3 Å and fixed motion zero; geometry still rejected |
| `1KTZ-upstream-two-torsion-study-attempt01` | 144 same-case φ7/φ4 trials selected unchanged input; conditional repeat in `1KTZ-upstream-wall-refinement-attempt01/not-executed.json` was not run |
| `1KTZ-selected-MD-api-control-attempt01` | Actual two-step selected MD interface/fixation control passed; not a scientific restoration sample |
| `1KTZ-selected-terminal-annealing-attempt01` | Four MD and six CG calls; all native observations fixed/global L, final geometry rejected |
| `1KTZ-observed-symmetry-audit-attempt01` | Original A/C, B/D and all 326 waters match one assembly operation within 0.0008 Å; observed-only evidence |
| `1KTZ-assembly-terminal-copy-attempt01` | Deposited BIOMT2 transfers 101 missing A-terminal proposal atoms to C; two Leu11 contacts and two terminal S–S twist warnings remain |
| `1KTZ-assembly-SS-wall-refinement-attempt01` | Explicit twofold torsion priors (k=100 kJ/mol) only for A7–16/C7–16, inward 2.8 Å walls and independent 3 Å gates; native fixed motion zero, final max 2.801 Å, narrow body/water/S–S checks pass |
| `1KTZ-restored-heavy-export-attempt01` | Actual full CLI exposes THR C60 reversed β stereochemistry, absent from narrow checks; preserved failed export |
| `1KTZ-missing-THR60-beta-review-attempt01` | Reviewed reflection of missing CG2 fixes β-R but introduces two weak contacts; retained negative result |
| `1KTZ-THR60-rotamer-review-attempt01` | 24 proper χ1 trials; ordinal 6 / +90° selected, zero new warnings, combined 2.70035 Å motion relative to the pre-reflection parent |
| `1KTZ-restored-heavy-export-attempt02` | Full PDB/DatRecord export and independent custody readback pass; exactly two unresolved OXT errors remain |

All graph-based rotations refuse observations/cyclic cuts. Finite walls are not
hard constraints or continuous-path evidence; native returned-stage audits precede
coordinate materialization. BIOMT transfers a proposal, not measured missing
coordinates or an independently acquired donor structure. Grids are same-case
conformational trials, not independent structures or scientific seeds.

The first rotamer reader incorrectly compared post-rounded intermediate PDB
distances against unrounded prior coordinates. Its failing source/report remain
retained; the corrected independent reader separately proves exact reflection,
three-decimal intermediate output and subsequent exact proper rotation. There
was no resampling or relaxation of scientific gates.

Trailing original waters sharing protein chain IDs also exposed a PDBFixer
missing-terminal false negative. `diagnose.structural.check_missing_terminals`
now supplements only water-terminated otherwise-free protein ends, preserving
caps, cyclic bonds, unresolved/future SEQRES gaps and case-sensitive chains.
Nine controls include actual PDBFixer behavior. Old masked reports remain intact.

The [latest independent export receipt](../../tmp/mac-real-gap-continuation-2026-11-14/1KTZ-restored-heavy-export-attempt02/readback-receipt.json)
confirms 4,022 heavy atoms: 2,986 unchanged protein observations, 326 unchanged
water oxygens and 710 restored atoms in 92 residues. All original context lines
are an ordered byte-identical subsequence. Only new atom serials are reassigned;
three labelled target SSBOND additions supplement the 18 original records.
Canonical DatRecord contains only the 710 restored atoms; 454 Cα centres remain L.
All eight final CLI warnings are source findings; **B136/D136 missing OXT are two
real errors**, not baseline warnings. Two temporary FF OXT atoms were never
published and remain outside the frozen mask.

The candidate is incomplete/not MD-ready until free uncapped receptor C-termini
versus explicit amidation/capping is confirmed. This does not reopen the confirmed
native oxidized TGF-β3 target. No H/protonation/charge or expanded completion mask
has been silently assigned. Static backbone screen has zero alarms, not knot or
trajectory qualification. Broader cohort work remains pending.

Software: final fast suite **2,848 passed, three optional skips, 33 slow deselected,
two warnings, 299.81 s**. An earlier focused launch lacked `dvbfixer` on PATH and
failed two subprocess controls; the corrected PATH launch and final suite passed.
Completed logs are retained in `software-verification/`; scientific receipts do
not bind live logs. No commit/push or destructive cleanup was performed.

## WP0 — publication check repaired

The post-publication license-inventory failure came from three newly tracked
backend patch files, not changed dependency versions. Regeneration adds their
declaration hashes and pending-source-review records (810→813 records), without
claiming legal clearance. `--check` and all12 inventory tests pass.

The complete fast lane subsequently passed **2474 passed, three optional
Torch/Biotite skips,33 slow deselected**,300.88 seconds. Evidence:
`tmp/mac-real-gap-continuation-2026-11-12/wp0-fast-suite-attempt01.log`.
New census/source-inventory tests added afterwards are checked separately; this
full-suite count must not be presented as including tests created after collection.

## WP1 — real metadata obtained, not merely a proposed list

`scripts/collect_real_gap_candidates.py` queried original public RCSB metadata,
retaining exact requests/responses, timestamps, hashes and manifest-last receipts.
The selection uses a fixed hash rank across five **total-entry-missing-count proxy**
strata plus previously inspected controls. These query strata are not falsely
reported as longest-gap strata. No private sequence was uploaded to an MSA service.

All **1000 distinct entries** were inventoried;995 have reported protein natural
missing spans. Longest annotated spans:366 entries at1–8,381 at9–24,131 at25–49,
52 at50–99 and65 at100+ residues. Five entries with no matching protein-gap feature
remain recorded rather than silently replaced. These are deposited-instance
annotations, not verified biological assemblies or scientific admission.

The old `rcsb_cluster_membership` field returned null throughout. A separately
retained modern `rcsb_polymer_entity_group_membership` query supplies30%-identity
groups for all protein entities in977 entries. Exact sequence plus shared30%
clusters produce **715 conservative connected entry components**; the largest has81
entries. This is not proof of independent folds/training absence. Incomplete
membership remains explicit; no final80-entry holdout has been assigned.

The metadata inventory contains difficult cases up to an annotated1066-residue
missing span. Whole unobserved chains, long tails, expression/construct ambiguities
and source/feature disagreements must be distinguished during coordinate audit;
they are not automatically supported reconstruction tasks.

### Original source and assembly inventories completed

A frozen **412-entry reserve pool**, spanning346 conservative components, was
downloaded as unmodified original gzipped mmCIFs. Reserve counts are75/110/
110/52/65 in the same increasing longest-span bins. All52 available50–99 and all65
100+ candidates are retained, not replaced with easy short loops. Small1BRS does
not count toward the primary cohort. This reserve pool is not the final200 entries.

All412 downloads and label-axis inventories completed, with no inventory refusal.
Deposited canonical sequences and absent-coordinate spans agree with metadata for
every chain; all sources have one coordinate model. Construct boundaries, chemistry
and source geometry are not yet admitted.170 entries have partial/zero occupancy
before residue-coherent altloc selection. The50–99 bin has only39 entries without
such observations, so additional screening may be needed rather than rewriting
occupancy or lowering the40-entry collection target.

7N8X has an812-residue internal span (831 observed/1643 target residues);7ZZ1
includes a1066-residue N-terminal span (618 observed/2286 target residues). These
need construct/sequence review: deposited absence does not establish that the full
terminal sequence was present in the experimental construct. Huge assemblies such
as9LY8 stay in the ledger even if PDB limits or host resources prevent admission.

`scripts/real_gap_source_inventory.py` inventories exact deposited label positions
against coordinate records: no greedy sequence matching, PDB conversion, chain-ID
truncation, artificial deletion, occupancy rewriting or coordinate movement.
Zero/partial occupancy, multiple models, author chain case/insertion codes,
unobserved whole chains, terminal construct review and assembly operators remain
explicit. Final assembly rendering, missing-heavy/chemistry checks and scientific
admission are still required.

Resource declaration: one network request at a time,60s/request, up to three
retained attempts,16MiB compressed/128MiB decompressed per-source read ceilings,
2GiB evidence cap and2GiB disk headroom. These are acquisition resource guards,
not maximum gap lengths. No original source is edited or overwritten.

Evidence under `tmp/mac-real-gap-continuation-2026-11-12/`:

- `candidate-census-attempt01/`:1000-entry metadata census and HTTP receipts.
- `candidate-groups-attempt01/`: modern grouping responses and connected components.
- `original-sources-attempt01/`: frozen reserve selection, original downloads and
  per-source inventory; completion is not inferred from directory existence.

`scripts/preflight_real_gap_assemblies.py` completed seven guarded CPU batches over all
412 originals using the shared CIF boundary, metadata-driven assembly1 renderer
and residue-coherent altloc selection. It retains original chemical context and
reports missing observed heavy atoms/explicit protein–nonprotein links, without
repair, occupancy rewriting or generation. The executor's `seed-0` suffix is only
a reservation key, not an ML seed. Source refusals are retained per entry; missing
assembly1/PDB capacity/chemistry/occupancy problems are not bypassed by cropping.

Of412 entries,380 produced complete assembly/altloc inventories and32 retain
refusals (18 normalization,11 altloc selection, three assembly). These are source
pipeline/format outcomes, not failed gap generations. Examples include modern
five-character nonprotein residue codes that cannot fit the current PDB boundary,
assembly serial capacity and ambiguous residue-coherent altloc selection.

Among380 assembly inventories,140 retain selected partial/zero-occupancy protein
atoms,78 need noncanonical polymer/ATOM chemistry review and190 have explicit
protein/nonprotein links. These sets overlap. Only118 have none of these initial
inventory obstacles; this is **not** permission to call the118 the intended200-entry
benchmark or discard the remaining source reviews.

The1HE8 smoke confirms915 observed canonical protein residues, no missing observed
heavy atoms,39 HETATM records and two explicit SER17/THR35-to-Mg201 links. LINK/CONECT
existence alone does not establish covalent bonding. The first smoke's overbroad
wording is superseded by separate attempt02 coordination-review details; both
records and unchanged original source remain retained.

`scripts/classify_real_gap_context.py` retained original mmCIF `struct_conn`
annotations for all412 sources:166 entries annotate protein–metal coordination
and49 have other cross-boundary link types (overlapping sets). Cross-boundary
records comprise3087 `metalc`,428 `covale` and two `hydrog` annotations. These are
deposited labels, not new chemistry/force-field inference. Metalc is not blanket
rejected as covalent; absent annotations do not prove absence of a bond.

### Source preparation:29 completed, one retained timeout, continuation declared

The118 currently obstacle-free canonical projections are separately declared for
guarded preparation using the shared sequence-placement/request builder and
`prepare_input` hook. Original full context remains retained; protein projection
and omitted-heterogen counts are explicit, not whole-entry conditioning/FF coverage.
Observed-heavy completion/bounded source repair preserve their existing evidence,
movement gates and ownership, separate from future generated residues. Existing
public-profile scope results remain recorded; the private input-readiness marker
is not a registered capability or backend admission.

Each case has an1800s/RSS8GiB serial CPU guard and the study has a10GiB output cap.
No checkpoint loads; request seed7 and executor reservation seed0 are not consumed
neural samples. Final200 selection and full-context/backend admission remain
pending. Evidence: `source-preparation-attempt01/`.

The first29 cases completed source preparation with valid geometry/global L,
including1AHW (1294 targets),1I2M (39-residue gap),9CQD (40-residue tail) and9XWG
(3093 targets). This is private protein-input readiness, not generated-gap success.
Assembly1 can be much smaller than the deposited unit: cases such as1U5G
(122 target residues),2BEB (186) and3EBK (176) are small controls, not evidence of
large-assembly readiness or automatic members of the primary200.

Case6MYX hit the unchanged1800s limit (peak RSS1313587200 bytes) while seeded
PDBFixer was completing372 missing heavy atoms in **observed residues**. The full
assembly has8400 targets/8076 observed residues in12 chains; its48 actual gaps are
only3/4/7/13 residues. Thus this is a source-preparation/context scaling failure,
not evidence for a maximum reconstructable gap length. Partial evidence remains;
no accepted result or continued launches were inferred from the workdir.

`source-preparation-continuation-attempt01/` explicitly binds/reuses all29 completed
results and retains6MYX as blocked without retry or timeout increase. It declares
the remaining88 unchanged first executions. This continuation is waiting for the
1HE8 trials to release the shared executor; no completed source is recomputed.

Independent request placement readback (`source-placement-readback-attempt01/`)
matches all29 prepared entries to exact deposited mmCIF label positions via
sequences and author residue/insertion-code evidence. Assembly chain renaming is
allowed, arbitrary repeat placement is not. All29 match; three small assemblies
are explicitly marked controls. One resource block and88 unexecuted jobs remain
separate, not imputed passes. New source/placement, continuation-safety and passive
reporter regression checks pass23 tests separately from the earlier2474-test full
lane; Ruff passes.

## WP2 — actual retained 1HE8 outcomes independently checked

`scripts/audit_1he8_joint_relaxation.py` rereads original masks/frames, final records,
geometry, chirality and movement, without sampling or optimization:

- Full patchwork: all1712 generated atoms plus148 observed atoms in19 anchor
  residues were mobile together. The remaining Ramachandran finding is **A581**;
  three anchor residues exceed the historical movement bounds. Aggregate and
  residue findings are not counted as two independent physical defects.
- Richer-context control: all1712 generated atoms plus18 atoms in A895/A902 were
  mobile together. The remaining amide finding is A901–A902, with two movement
  violations under its separate tighter envelope.

Both preserve every nonmobile record and global L-chirality; neither passes its
unchanged combined criteria. No crystal-RMSD gate caused rejection.

Historical records have phase energy/iteration ceilings, but **not actual iteration
counts, termination reasons, final gradients or a saved H-containing force state**.
Thus convergence and an exact historical force diagnosis cannot be inferred.

A separately declared guarded inspection reconstructs temporary hydrogenated
systems around saved final heavy atoms, audits actual particle masses and
decomposes force components. It aborts before the main heavy-atom minimization;
temporary hydrogen placement still includes OpenMM's internal H optimization.
This is not a historical H-state replay or a newly repaired candidate. Both guarded
inspections completed: `1HE8-force-state-attempt01/`. Actual temporary particle
masses confirm1712 mobile generated atoms and148/18 mobile observed anchor atoms;
7247/7377 other observed heavy atoms respectively are massless/fixed. This rules
out a missing generated-atom mobility mask in the reconstructed systems, not every
possible historical topology/numerical issue.

Per-force energies/vectors, atom identity/index/mass roles, temporary topology,
original-group system XML and reconstructed state/force arrays are retained. The
two-anchor system has active auxiliary amide forces at A901–902 and nonzero anchor
boundary forces; it is not a system where those terms were absent. Forces from
different terms oppose one another, but this alone does not establish the cause
of the old failed endpoint. Rebuilt hydrogens can substantially change force
values; one full-patchwork reconstructed A901 N has a large nonbonded contribution,
which must not be misreported as the historical final gradient. A passive actual
iteration/gradient recorder and prospective controlled mobility trials remain
necessary before causal claims or changes to restraints/thresholds.

### Traced mobility comparison completed and independently audited

`1HE8-traced-mobility-attempt01/plan.json` now freezes **eight CPU trials**:
two retained pre-refinement protein-only proposals × fixed/direct/sequence±2/
spatial5Å-plus-one-neighbor policies. Every trial moves all1712 generated atoms
together; observed mobility uses the original source frame, uniform k10000,
the historical0.5Å backbone/0.5Å Cα/1.25Å heavy envelope, seed7, no restarts and
1000 full-phase iterations after the same preconditioning. No crystal scoring or
new neural proposal is used. These are new controls, not rewrites of the tighter
historical two-anchor attempt; they do not establish whole-entry cofactor coverage
or manually approve new RFD3 donors.

The passive OpenMM minimization reporter records actual callback iterations,
objective-gradient statistics and snapshots without returning an early-stop
request or restarting the optimizer. A separate analytic Reference-platform
control gives an identical minimum with/without the reporter; this does not claim
bit-identical GPU trajectories. Phase-return force states are saved, while an
engine termination reason is not invented from an API that does not expose it.
All eight trials completed and independent readback verified the final geometry,
nonmobile records, global L and the actual XML wall parameters/original reference
coordinates (`1HE8-traced-mobility-independent-readback-attempt02/`). Combined
acceptance is **0/8**; geometry passes **3/8**, all richer-context mobile arms.
Richer-context direct/sequence2/spatial backbone maxima are0.542530/0.519114/
0.516911Å versus the unchanged0.500Å limit. Solver-state comparisons confirm true
exceedances, not merely PDB rounding. All original outcomes remain rejected.

The boundary wall is a finite quadratic penalty, not an exact constraint. At the
richer/direct A902 solver backbone displacement0.542541Å its energy is only about
9.049kJ/mol. Four new engineering controls are declared in
`1HE8-buffered-walls-attempt01/`: original-patchwork/direct negative control and
richer direct/sequence2/spatial. Wall activation moves inward by0.05Å to
BB0.45/Cα0.45/heavy1.20Å; **outer acceptance stays0.5/0.5/1.25Å**. Same proposals,
seed7, k10000, CPU, no restarts and1000 full-phase iterations. This is not a new
physical uncertainty estimate or a retrospective relaxation of acceptance.
These controls await the shared executor, now running the88 source first
executions.6MYX is not retried.

### A581 retained-frame diagnosis

`scripts/analyze_1he8_ramachandran_trace.py` follows request-axis neighbors and
verifies the phi/psi bonds against each serialized temporary FF system. In all
four original-patchwork arms, A581 initially has defined angles157.47/−32.41°,
which pass the unchanged gross-outlier reference gate. The **first retained**
outlier is phase3 callback50: the full-strength soft-core phase (`soft_scale=1`),
before the final full-force-field phase. Final phi/psi is about163/−26°.
Thus the recorded defect was not an initially undefined dihedral merely becoming
measurable. Sparse frames do not locate its exact first transition or establish
why the local objective favors this basin. Richer-context A581 has no outlier in
any retained frame; this is not a continuous no-passage proof. No crystal reference
was used. Evidence: `1HE8-A581-frame-readback-attempt01/`.

### First source continuation stopped; buffered controls launched

The88-task source continuation executed26 completed cases before another retained
1800s stop:25 private protein-source preparations passed and5DW7 was refused
because bounded observed-source repair exceeded its movement envelope. The refusal
is not a failed neural sample and is not counted as ready. Independent deposited
axis readback matched all25 successes (`source-placement-readback-attempt02/`).
Together with the original29 there are54 ready projections, one source-repair
refusal, two time-blocked cases and61 unexecuted tasks of the original118.

7SFS stopped after1800.240s, peak RSS1417691136 bytes, during completion of24
missing observed heavy atoms. Its assembly contains10416 target residues in24
chains and61 gaps (lengths4/8/10/20/37), not one exceptionally long gap. It joins
6MYX as a retained source-preparation resource block; neither is retried.
`source-preparation-continuation-attempt02/` binds the preceding evidence and
declares only61 unchanged first executions, with no increased budgets.

The54 source-ready projections currently have longest-gap counts18/26/10/0/0 in
the1–8/9–24/25–49/50–99/100+ bins. Thus this partial prefix of the task order is
not the promised hard cohort or evidence of long-gap admission. Remaining cases
and original context/construct reviews must be completed before selection.

The shared executor now runs the four declared buffered1HE8 controls. Their
launch does not imply success; original0/8 combined outcomes and unchanged outer
acceptance remain preserved. The61-source continuation waits for executor release.

### Buffered controls completed: two joint geometry-and-movement passes

All four controls completed and were independently read back, including actual
serialized wall parameters/references, exact nonmobile records and global L.
Evidence: `1HE8-buffered-walls-independent-readback-attempt01/`. The newly declared
series has **2/4 combined passes,3/4 geometry passes**; historical0/8 is unchanged.

| Proposal / released observations | Geometry | Max BB RMS (Å) | Max Cα (Å) | Max heavy (Å) | Combined |
|---|---|---:|---:|---:|---|
| Original patchwork / direct | A581 outlier |0.503310|0.479884|1.211298|Reject|
| Richer context / direct |Pass|0.506037|0.452306|1.010349|Reject|
| Richer context / sequence±2 |Pass|0.486972|0.445481|0.903109|Pass|
| Richer context / spatial5Å+padding1 |Pass|0.473990|0.340359|1.019724|Pass|

These maxima are computed from final PDBs against the original observed frame.
The two passes also have zero movement violations at solver precision; direct
still violates at solver precision (A902 BB0.505964Å), not just PDB rounding.
All1712 generated atoms were mobile together; sequence/spatial released57/228
observed residues, respectively, under the declared restraints. Actual full-phase
callback counts were437/293/229/242 in table order; return alone is not an exposed
engine convergence reason or a continuous trajectory safety proof.

This establishes two usable endpoints under this private protein-only engineering
contract, not a general anchor-uncertainty policy, new manual donor approval,
whole-entry Mg/ligand coverage, topology/trajectory qualification or public product
readiness. Inward finite-wall activation is still not a guarantee of feasibility,
as the failed direct controls demonstrate. No additional guard-band tuning,
restarts or resampling is silently added. Next: complete the large real-input
admission ledger and design independent validation before broad method claims.

The61 unchanged source first executions started after these controls released the
shared executor (`source-preparation-continuation-attempt02/`). Completed outcomes
and both resource-blocked cases remain preserved and are not repeated.

### Second source continuation:69 ready, two refusals, three time blocks

The next16 completed outcomes contain15 ready protein projections and one retained
source-preparation refusal (5W3F: selected repair residues contain measured
hydrogens). The existing preparation contract refuses modifying such observations;
no silent H stripping or occupancy/identity rewrite was introduced. Independent
deposited-axis readback matches all15 ready results; one additional small assembly
is marked a control (`source-placement-readback-attempt03/`). Total verified ready
placements are69, including four small assembly controls, not69 large-case claims.

6PBF stopped at1800.082s, peak RSS655097856 bytes, while completing808 absent heavy
atoms in observed residues. It joins6MYX/7SFS as a retained time block without retry
or budget increase. `source-preparation-continuation-attempt03/` binds all prior
evidence and runs only the44 still-unexecuted tasks under unchanged limits.

A new cumulative read-only ledger verifies that settled outcomes plus the next
plan partition all118 original source jobs exactly once. It verifies ready source
and request hashes/containment, preserves refusals and guards, and does not inspect
running output directories or count their existence as success. Snapshot evidence:
`source-cumulative-ledger-attempt01/` —69 ready, two refusals, three time-blocked,
44 pending at the snapshot. Ready longest-gap bins:18/26/23/2/0. The two50–99
projections are input-readiness evidence, not sampling/long-gap qualification.
This remains insufficient for the200-entry primary cohort and its hard quotas.

### Third source continuation:90 ready,12 refusals, four time blocks

Of31 completed cases in the next series,21 passed private source preparation and
independent deposited-axis readback;10 were refused. Across all settled series the
118-task ledger now partitions90 ready,12 refused, four time-blocked and12 pending
at snapshot (`source-cumulative-ledger-attempt02/`). Ready longest-gap bins are
18/26/23/13/10. The fourth continuation runs only the final12 first executions.

The new refusals remain separately classified:

- Measured H in selected repair residues:1AK5,2JDW,7RFH.
- Source chain breaks beyond automatic repair scope:6VM2,6Y9V,8TUS.
- D-Cα reported after heavy-atom preparation attempts:4ACC (A/B GLN385),6KXW
  (B GLN30). No sampling occurred. This is an observed failure report, not proof
  that PDBFixer created D geometry or a resolved causal diagnosis; original and
  intermediate observations must be compared before assigning cause.
- Explicit links involving an incomplete residue:7A6C.
- Request residue-identity allocation:8ZU3. The current allocator preserves
  observed IDs, tries numeric space between anchors, then at most26 uppercase
  insertion codes at the left anchor. This is a numbering/representation
  limitation, not an ML gap-length or generated-geometry failure. It is not
  bypassed by silently renumbering observations or changing the frozen worker.

6RN2 reached1800.298s, peak RSS868188160 bytes, while completing44 missing observed
heavy atoms. It is retained without retry alongside6MYX/7SFS/6PBF.

Importantly,7N8X's authentic812-residue internal gap survived source preparation
and independent deposited-axis audit:1643 targets,6617 fixed atoms and6407 future
generated atoms, with one gap and global L. Its prepared observed derivative is
explicitly bound in the request and source repair evidence remains retained.
**No coordinates for the812-residue gap have been generated**. A generic private
scope warning is not per-backend admission or qualification at this length;
construct/context review and explicit sampler task/resource preflight still apply.

A read-only inspection of the retained pre-preparation protein projections further
limits interpretation of the two D-Cα refusals:4ACC A/B GLN385 each contain only
the N atom, while6KXW B GLN30 has CA/C/O/CB/CG/CD/OE1/NE2 but no N. Thus neither
affected source site has the complete N–CA–C–CB quartet required to define the
canonical signed-volume test. `find_d_residues` finds no D centre in these source
projections, but that cannot be described as proof these incomplete centres were
initially L. Their later rejected D geometry is real; the causal mechanism and
possible safe completion method remain unestablished. No source or failed repair
was changed and no D endpoint was admitted.

### Initial118-task source series fully accounted; recruitment expanded

The final12 first executions completed:11 prepared projections independently
matched deposited axes and20ZC was refused because bounded source repair exceeded
observed-residue movement limits. The complete ledger
`source-cumulative-ledger-final-attempt01/` accounts for all118 tasks exactly once:
**101 private source-ready,13 refused, four time-blocked, zero unexecuted**.
All101 ready placements match deposited axes across five readback series; four
remain explicitly small-assembly controls. This is not101 qualified ML inputs.

Ready longest-gap bins are18/26/23/13/21. Even before construct/context and group
admission, this cannot satisfy the200-entry quota30/50/50/40/30. No final80-entry
evaluation set or backend sampling matrix is frozen. The cohort target and hard
strata are not reduced to fit the easy canonical subset.

`scripts/expand_real_gap_source_reserves.py` declares and starts acquisition of
**181 additional originals** from unused entries of the existing1000-entry grouped
census:60/100/21 in the1–8/9–24/25–49 bins. Previously downloaded412 IDs are excluded;
fixed hash ranking favors less-represented conservative components. Selection,
code and prior ledger are hash-bound before requests. This is serial public source
screening only, with unchanged16MiB compressed/128MiB decompressed per-source read
ceilings,60s/request, at most three retained HTTP attempts and2GiB evidence cap.
No checkpoint or new neural task is declared. Evidence: `source-expansion-attempt01/`.

All50–99/100+ entries in the current census were already acquired in the original
reserve pool. The expansion therefore **does not replenish those strata** and
cannot be represented as doing so: their original chemistry/occupancy/context
reviews and targeted additional metadata recruitment remain necessary, especially
for the40-entry50–99 quota and group-safe holdout. Failed source cases and original
inputs remain preserved, without timeout retries or silent observation changes.

### Additional181 sources inventoried;41 preparation tasks launched

All181 additional original downloads and label-axis inventories completed without
inventory refusal. This brings original-source coverage to593 unique entries from
the existing1000-entry metadata census. Separate connection annotation readback
found71 entries with deposited protein–metal coordination and19 with other
cross-boundary links (overlapping;879 `metalc`,82 `covale` records). This remains
annotation evidence, not chemistry/force-field admission.

Assembly/altloc preflight completed all181:163 inventory successes,18 retained
refusals (15 normalization, two altloc selection, one assembly). Among successful
inventories,91 have explicit protein/nonprotein links needing classification,
28 have noncanonical chemistry and80 retain partial/zero-occupancy observations;
sets overlap.41 have none of these initial obstacles and are separately declared
for the unchanged shared source preparation (1800s/case, RSS8GiB, serial guarded
CPU). None is yet counted as prepared simply because it launched. Evidence:
`assembly-preflight-expansion-attempt01/`, `context-annotations-expansion-attempt01/`
and `source-preparation-expansion-attempt01/`. No original failure is retried.

### Targeted additional metadata recruitment declared and running

The existing long-stratum shortfall justifies the plan's expandable screening
pool. `additional-hard-metadata-attempt01/` declares at most500 **new public
metadata entries**, excluding all1000 previous IDs. No coordinates, neural tasks
or scientific seeds are added by this acquisition. Three **total-entry missing
count proxy** strata recruit250 at50–99,150 at100–199 and100 at200+; these are not
claims about longest natural-gap lengths. The upper stratum has no length ceiling.
Fixed hash ranking, exact HTTP receipts and a pre-request hash-bound declaration
retain failed requests and refuse duplicates/undersupplied strata rather than
substituting silently. Resources: serial requests,60s/request, three attempts,
32MiB response ceiling,512MiB evidence cap and2GiB disk headroom.

This is screening expansion to a potential1500-entry metadata pool, not a changed
200-entry target or authorization for a large sampler matrix. Source/construct and
backend admission, hard quotas and **union-wide** sequence-component grouping are
still required before any final development/calibration/evaluation assignment;
grouping only the new subset would not prevent leakage to prior entries.

### Additional metadata and union grouping completed;140 hard originals acquiring

All500 new metadata records completed without inventory refusal or previous-entry
duplication. Their actual longest-gap bins are19/200/141/63/77, demonstrating why
the total-missing search proxies cannot be reported as actual gap lengths.
Additional grouping supplied complete30%-identity memberships for486/500 entries.
`scripts/merge_real_gap_candidate_groups.py` then recomputed the connected graph
over the **union**, rather than concatenating either subset's local component IDs.
Regressions cover a new record bridging two old components and duplicate refusal.

`union-candidate-groups-attempt01/` retains1500 entries,1016 components and1463
entries with complete30% memberships. Largest component is176 (old largest81;
new-subset largest51). This is a conservative sequence-key graph, not a biological
family taxonomy or proof of model-training absence. The reduction from715+385
local components to1016 union components confirms cross-subset links; no split
was frozen before this recomputation. Incomplete membership remains explicit.

Union longest-gap bins are385/581/272/115/142 (five entries have no matching
protein-gap feature).50–99 spans cover115 entries in86 components;100+ spans
cover142 entries in76 components. These are metadata counts, not ready structures
or available independent scientific tasks.

`hard-original-sources-attempt01/` now freezes and acquires **140 unused50+ sources**:
63 at50–99 and77 at100+, excluding every previously acquired593 entry. Original
mmCIF bytes, source-axis inventories and failures use the unchanged acquisition
guards; the outer declaration binds union metadata, previous ledgers and worker
code before requests. Budget is140 originals (hard refusal if a declaration would
exceed200), serial60s requests, three attempts and2GiB evidence ceiling; no upper
gap length cap. This is source recruitment only, not neural sampling or backend
admission. The independent41-case source preparation expansion continues under
the serial CPU executor. No running output directory was polled for readiness.

### All140 hard originals acquired; assembly study declared and waiting

All140 original downloads and label-axis inventories completed with zero refusals
or source/metadata gap disagreements. Each has one coordinate model. All733
acquired original entries across the three acquisitions are unique by declaration;
source-ready counts remain separate. Before residue-coherent altloc selection,
46/63 of the new50–99 sources and66/77 of the100+ sources have no partial/zero
occupancy observations. These do not yet pass assembly, chemistry or geometry.

Original connection annotations cover51 entries with protein–metal coordination
and23 with other cross-boundary link types (overlap;832 `metalc`,512 `covale`
records). Link annotations do not establish force-field support or automatic
covalency. Evidence: `context-annotations-hard-attempt01/`.

`assembly-preflight-hard-attempt01/` declares three unchanged guarded assembly/
altloc inventory batches for all140 sources. They wait for the41-case source
preparation expansion to release the shared serial CPU executor; no second CPU
study is launched against its lock. No crop/renumber/occupancy changes or neural
sampling are performed by this acquisition.

### Source expansion stopped at8V3Z; hard assembly preflight launched

The41-case source expansion completed five outcomes before a retained timeout:
four prepared projections (1CC6,2R96,34FS,6J9V), all independently matched to
deposited axes, and a3FIM refusal for source chirality outside automatic-repair
scope. No chirality repair was silently added to this contract.

8V3Z reached1800.437s, peak RSS1374519296 bytes, during bounded **observed geometry
repair**, not gap generation or missing-heavy completion. The report retains two
severe input findings and266 warnings; the attempted locality was four observed
residues, with other source records fixed. The case is time-blocked, not accepted
or diagnosed as an ML failure. Original report and partial repair evidence remain.

`source-preparation-expansion-continuation-attempt01/` binds/reuses five completed
outcomes, retains8V3Z without retry or budget increase and declares only35 unchanged
first executions. `source-expansion-ledger-attempt01/` partitions the41 tasks into
four ready, one refusal, one timeout and35 pending. Across the original118 and this
41-task expansion, totals at this checkpoint are105 privately ready/axis-matched,
14 refusals, five time blocks and35 unexecuted; no combined final cohort is implied.

After the source executor released its lock, the three declared140-source hard
assembly preflight batches launched. The35-source continuation waits; there is
still only one guarded CPU study active. Neither the deposited metadata warnings
nor partial workdir creation are used as scientific completion evidence.

### Hard assembly inventory completed;39 source preparations declared

All140 hard-source assembly preflights completed:120 inventory successes and20
retained refusals (16 normalization, two altloc, two assembly). Among successes,
60 have explicit protein/nonprotein links needing review,24 have noncanonical
chemistry and21 retain partial/zero occupancy; sets overlap.39 have none of these
initial obstacles and are declared, **not launched**, under
`source-preparation-hard-attempt01/`. No source/geometry/backend pass is inferred
from this declaration.

The35 unchanged source first executions from the shorter-span expansion launched
after hard assembly preflight released the shared CPU lock. The39 hard-source
preparations wait for that executor. No completed task or timeout was repeated.

A recruitment feasibility check shows why the current program cannot be called
the200-entry benchmark: initial101 ready + at most39 additional successes from
the41-case expansion after its known refusal/timeout +39 newly declared hard
preparations gives an **absolute upper bound179**, even if every remaining task
passes. Four small controls and subsequent context/group admission reduce it
further. The union census still has141 unacquired25–49 candidates in122 sequence
components. Further bounded recruitment and context reviews are required, not a
target reduction, reclassification of refusals or declaration of completion.

### Source expansion finished; hard preparation running; medium recruitment launched

The35-case continuation completed all outcomes:28 ready, seven refusals and no
new time block. All28 ready requests independently match deposited axes; three
additional small assemblies are marked controls. The full41-task expansion now
has32 ready, eight refusals and one retained8V3Z timeout, with zero unexecuted tasks
(`source-expansion-ledger-final-attempt01/`). Across the two settled programs
(118+41 tasks):**133 ready/axis-matched,21 refused, five time-blocked**, including
seven small assembly controls among the133 ready. Longest-gap bins of the ready
projection union are26/50/23/13/21, not final cohort assignments.

New refusals are1YH3/3ZDN (D-Cα detected after completion attempts),4G4Y (source
chain break),5C8W/5IRV/8CNI (observed H in affected residues) and7FBO (damage exceeds
the bounded local repair mask). These outcomes do not establish D-Cα causality,
authorize global relaxation or permit silent hydrogen stripping. Existing failures
and8V3Z are not retried. The39 declared hard-source preparations started after
this continuation released the serial CPU executor.

Because25–49 readiness remains23 versus its50-entry quota, and the current
program cannot reach200 even ideally, `medium-original-sources-attempt01/` now
declares and acquires all141 unused25–49 metadata candidates in122 union
components. All733 prior original acquisitions are excluded. This is an explicit
recruitment stratum, not a method length ceiling or replacement for either long
quota. No coordinate edits, new neural task or final split is declared. Source
acquisition retains the existing serial60s/three-attempt/2GiB evidence guards and
16MiB compressed/128MiB decompressed per-source read ceilings. Source inventory,
assembly, chemistry and geometry admission are still required after download.

### All141 medium originals inventoried; assembly preflight queued

The141 new25–49 originals completed download and label-axis inventory, with zero
inventory refusals, source/feature disagreements or multiple coordinate models.
Total original-source custody is874 unique entries; this is not a readiness count.
Before residue-coherent altloc selection,59/141 medium sources contain partial or
zero-occupancy protein observations. Their classification remains explicit.

Connection annotations were independently retained for all141:48 entries have
protein–metal coordination and21 other cross-boundary connection types (overlap;
976 `metalc`,176 `covale` records), without force-field/covalency inference.
Evidence: `context-annotations-medium-attempt01/`.

`assembly-preflight-medium-attempt01/` declares three guarded assembly/altloc
inventory batches. They are **not launched** while the39 hard-source preparation
jobs own the serial CPU executor. No current hard-preparation result is inferred
from directory existence or polled progress. The medium acquisition itself did
not generate gap coordinates, rewrite source observations or freeze a final split.

### Ready-cohort capacity audit and correction of textual100+ overcount

`scripts/audit_real_gap_cohort_capacity.py` independently rereads the159 settled
tasks, both immutable ledgers, all seven deposited-axis audits and source-bound
requests. It verifies current source/request hashes, derives gap kinds/lengths
from target intervals rather than author-number adjacency, and requires longest
bins to partition all ready projections exactly once. No running output is read.
Evidence: `cohort-capacity-readback-attempt02/` (the initial audit is also retained).

This found an error in the **hand-maintained text/status**, not the original
artifact ledgers:100+ ready projections were stated as32 instead of21. The initial
118-task ledger correctly had21;32 was the total ready count of a different
shorter-span expansion. The mistaken32 was propagated into a later textual union
bin count. Text and status are corrected; original requests, outcomes, ledgers,
total133 readiness and all scientific thresholds are unchanged.

After excluding seven explicitly small assembly controls,126 ready projections
remain, spanning112 union components. Two lack complete30% group memberships;
93 are multichain,115 multigap and40 have at least1000 **observed** residues in the
selected biological assembly. These are descriptive input strata, not method
qualification. There is no rule silently limiting selection to one entry per
component; component custody must instead prevent leakage across final splits.

| Longest gap | Target quota | Non-small ready projections | Minimum shortfall before further admission |
|---|---:|---:|---:|
|1–8|30|23|7|
|9–24|50|47|3|
|25–49|50|22|28|
|50–99|40|13|27|
|100+|30|21|9|

118/126 have at least one terminal missing span, so construct/expressed-sequence
review remains necessary. This does **not** mean118 are invalid, or that every
longest gap is terminal; a short unresolved terminus also triggers review.
Twelve have a genuine internal gap of at least50 by deposited-axis placement,
including7N8X812,7CRH141,9EGO/7T10127,9JHJ126,5N9J109,7JHG104 and6ZBG86.
Internal missing spans also require scientific/construct admission; no sampling
or favorable crystal-pose assumption was introduced by this readback.

The combined continuation controls (`test_real_gap_*`, `test_1he8_*`, passive
minimizer and additional-metadata tests) pass **45 tests**. The new capacity
auditor/tests pass Ruff, and `git diff --check` is clean. These focused checks
are not a substitute for the still-pending full-suite rerun after all additions.

## Current boundaries

Current user scope amendment takes priority over the historical length quotas:
**do not include7N8X or analogous enormous/domain-scale omissions in the active
loop-completion test cohort**. The source-readiness records above are preserved,
not converted into cohort admission. Keep unclassified large spans pending and
retain huge cases only as a separate research reserve; no new proposal/refinement/
evaluation for those cases is authorized by their prior source preparation.
No universal numeric length cutoff or replacement quota table is silently
invented. See the current amendment at the top of the
[continuation handoff](../plans/mac-reconstruction-real-gaps.md) and
`cohort-policy-amendment-oversized-gaps-attempt01.json` in the retained evidence
root. Do not crop a frozen whole-entry request to hide a deferred gap.

The user further directs **work with what is already retained**. Recruitment is
paused: no new metadata/source acquisition or replacement structures are declared.
There are874 retained originals. The39 hard-source preparation tasks are existing
input inventories, not test admissions; the141 medium assembly inventory batches
remain declared and waiting for the serial CPU executor.

### Existing1HE8 endpoints: per-gap movement and junction readback

`scripts/analyze_1he8_joint_gap_endpoints.py` rereads all12 audited endpoints without
optimization/sampling and binds original requests, own starting proposals, final
PDBs and actual temporary FF bond graphs. Every generated atom is accounted for
exactly once by gap; all19 observed/generated outer peptide junctions are verified
against real FF C–N bonds. Evidence: `1HE8-joint-gap-endpoints-attempt01/`.

For the two accepted richer-context buffered arms:

| Metric | Sequence±2 | Spatial5Å+padding1 |
|---|---:|---:|
| Generated-heavy RMS motion from own proposal (Å) |0.502526|0.521569|
| Mean motion (Å) |0.365703|0.391343|
|95th percentile (Å) |1.065235|1.084959|
| Maximum atom motion (Å) |3.802355|3.835124|
| Final outer-junction C–N range (Å) |1.3330–1.3740|1.3324–1.3709|
| Maximum outer-junction deviation from nearest cis/trans ω (degrees) |20.818|29.179|

These are repair displacements, **not crystal RMSD, accuracy estimates or donor
ranking**. The own starting proposal has a maximum outer-junction ω deviation of
57.397 degrees. Geometry passed the unchanged gross-defect gate, not a claim that
every amide is ideal: the spatial endpoint retains a larger worst junction twist
despite its smaller worst observed-backbone displacement. Preserve this tradeoff
instead of automatically ranking either endpoint as universally better. All ten
gap motion summaries and individual before/after junction measurements are saved.
No continuous topology/no-passage guarantee is inferred from clean endpoints.

### Current loop scope and retained-input correctness tests

`current-loop-scope-attempt01/` applies the user amendment without changing
requests: of133 retained ready projections,7N8X is explicitly excluded, seven are
small controls and125 remain pending loop/context review. None is automatically
test-admitted, even for short gaps. Historical compulsory length quotas are not
current execution requirements; no new recruitment or numeric cutoff is added.

Two existing water-only protein projections,1NON (four18-residue internal gaps;
716 targets) and2BE1 (six8/10/20-residue internal gaps;686 targets), were used only
for weightless context-contract inspection (`retained-loop-context-contracts-attempt01/`).
Neither reads a hidden structural reference. Current Protenix full-input routes
refuse these sizes at the **adapter**512-residue guard, before weights or memory
measurement. The auto-selected local axes have354/365 residues. cc94 forms4/6
separate local axes with sizes114/100/96/83 and43/86/79/48/85/86, respectively.
These are different coupling/context methods, not matched inference results.

The shipped multichain Protenix template encoder was then exercised without
loading a model/dataloader (`retained-loop-protenix-encoder-attempt01/`). Its
14/22 engine entities are **fragments of physical source chains**, not14/22 native
biological chains. A separate reader verifies every entity sequence, placement,
case-sensitive original residue/insertion-code identity and all selected observed
atoms/coordinates, including terminal OXT. Source and derived-request staging
are byte-exact; only observed coordinates enter templates. Actual Protenix feature
construction and numerical checkpoint execution remain separate unexecuted gates.

Eight encoder regression cases pass, including negative checks for lost OXT,
coordinate changes, repeated atoms/templates, incorrect engine-chain identity,
reordered alignments and attempts to fill generated target slots with template
coordinates. Scope regressions prohibit implicit short-gap admission, promotion
of small controls, duplicate source counting and mutation of retained requests.

The user suggested model/checkpoint/platform-specific budgets, then explicitly
deferred limits work in favor of these other tests. No larger-input experiment,
production cap change or new sampling was started. Current work remains input/
mapping correctness on retained structures, plus the already-running inventory.

The complete focused continuation set now passes **61 tests**; Ruff on the new
scope/context/encoder code and `git diff --check` pass. The earlier full2474-test
lane is historical evidence, not a full rerun covering these additions.

### Native data-only checks declared at unchanged input sizes

The user requests continued autonomous correctness work and has deferred residue
limit experiments. `retained-loop-protenix-native-preflight-attempt01/` therefore
declares exactly **two native data-loader checks**, at existing354/365-token local
inputs, with zero checkpoint/sampler calls. This is not a larger-context study.
They wait for the shared executor; current declaration is not an execution pass.

`scripts/preflight_retained_loop_protenix_native.py` uses the owned pinned v1
source, existing assets and Kalign, not the different newer `engines/protenix`
checkout. Before upstream data imports it prohibits torch checkpoint loads,
module forward calls and socket connections. Data-loader preprocessing on CPU is
explicitly not a numerical inference fallback. Native feature tokens are mapped
from actual atom-to-token indices with source chain/insertion-code identity;
generated residues remain generated even if the native axis adds disposable OXT.
Partial-residue completion ownership is refused by this whole-gap-only check.

Template masks must exclude generated token rows/columns and cover observed
context tokens; a completely empty template cannot masquerade as safe zero
leakage. Four regression cases cover both mask axes, lost source coverage, token
aliasing, incomplete token indices, chain case/insertion codes and generated OXT.
Worker-module import was checked in the isolated native Python with explicit
project/src paths: neither OpenMM nor Torch is imported by module import alone.
Earlier probes without the executor's explicit PYTHONPATH failed import and are
not counted as native data execution. No native loader/checkpoint has run yet.

The shared serial CPU guard retains1800s/job,RSS8GiB, host-pressure/headroom and
10GiB evidence ceilings. No environment asset is downloaded or patched to make
these checks pass, and failure stops remaining launches without implicit retry.

### Autonomous serial continuation queue launched

To avoid requiring a user turn after every resource-stop notification,
`retained-correctness-queue-attempt01/queue-plan.json` declares and launches a
finite serial scheduler. It waits on the shared executor with a **blocking file
lock**, not progress polling. The stages are:

1. Account for the existing39 hard-source preparation jobs and continue only
   unchanged unexecuted cases after retained wall-time stops. No timed-out or
   completed case is repeated; technical/other-resource failures stop the queue.
2. Run the already-declared141-source medium assembly/altloc inventory.
3. Run the two declared native Protenix data-only correctness checks at existing
   local sizes. Any failure is retained and stops the remaining queue stages.
4. Declare canonical source preparations only from that retained medium inventory,
   bounded above by141 tasks, and continue their unexecuted cases under the same
   per-child guards. This is input inventory, not cohort admission or sampling.
5. Hold the shared lock for a full `pytest -m 'not slow' -q` rerun, retaining its
   full log and exit status with an1800s ceiling and child-process cleanup.

The scheduler has an exclusive owned reservation and does not adopt crash-left
directories. Verified completed non-source stages can be reused without rerunning;
partial native stages need an explicit recovery amendment. Result files must match
their completed execution receipts. A competing executor reservation fails closed.
New source acquisition, numerical model calls, model assets and limit experiments
remain **zero**. Completion is not claimed from scheduler launch: use its final
`completed.json` or retained `queue-blocked.json`, plus individual stage receipts.

New scheduler regressions cover completed-outcome reuse, retained time blocks,
technical stops and uncommitted directory refusal. Four additional tests exercise
the existing pinned Protenix context-axis helper: chain case/insertion codes,
coalesced same-sequence entity rejection, missing OXT and wrong native sequence.
These are negative correctness controls, not proof of native feature execution.
All new focused tests pass; the queue's full fast lane remains pending.

### Autonomous queue completed;194 source-ready projections independently checked

`retained-correctness-queue-attempt01/completed.json` confirms every declared
queue stage completed. The39 hard-source tasks are30 ready, seven refused and
two retained time blocks (9DJU/9SFN); only previously unexecuted cases continued.
The141 medium assembly inventories yielded129 successes and12 refusals (six
normalization, one altloc, five assembly). Their35 canonical preparation tasks
are31 ready, three refused and one retained7WMP time block. No timeout was retried.

All61 new ready projections independently match deposited positions, using five
new per-series readbacks under `source-placement-queue-readback-attempt01/`.
Its two final ledgers account for all74 new tasks. Across all233 settled source
tasks the result is **194 ready,31 refused, eight time-blocked**, zero unexecuted.
`cohort-capacity-readback-attempt03/` binds all four ledgers and12 placement audits.
There are12 small controls;182 non-small projections span156 union components,
four lack complete group memberships,171 have at least one terminal span,133
are multichain,170 multigap and59 contain at least1000 observed residues.

Applying the current user amendment (`current-loop-scope-attempt02/`) explicitly
excludes7N8X, retains12 small controls separately and leaves181 pending loop/
construct/context review. **194 ready is not194 scientifically admitted tests**.
The old length quota table printed by the historical capacity reader is descriptive
only; it does not reinstate superseded recruitment requirements or admit huge gaps.

Both native Protenix data checks completed on354/365 tokens and2763/2857 atoms.
Their72/76 generated tokens have zero structural template masks; every observed
context token remains covered. No checkpoint/model-forward call occurred in those
data checks. The full fast lane then passed **2546 passed, three optional skips,
33 slow deselected**,298.79s. Its exact log/exit receipt remain under the queue.
This precedes the later callback/sampling-adapter tests, checked separately.

### Synthetic MPS callback controls completed; legacy wiring limitation isolated

The legacy private MPS runner passes all source observations to its reinjection
constructor. For local contexts that correctly fails with "fixed coordinate
identity is absent from sampled atoms". It must not be confused with a failure of
the model, residue budget or native input features.

`retained-context-callback-controls-attempt01/` executed two explicit synthetic
MPS callback tests, not inference. The represented/preserved-only/rigid-alignment
counts are2193/2847/930 for1NON and2263/2595/1180 for2BE1. Passing only represented
observations and using the shared context alignment evidence restores fixed atoms
exactly at float32 precision; synthetic global-frame error was0 in both controls.
All source fixed PDB records remain byte-exact in an in-memory materialization,
including preserved-only atoms. No synthetic structure is saved as a candidate,
no diffusion step occurred, and the frozen MPS runner remains unchanged.

### Two private checkpoint-backed development controls declared

The next frozen study is limited to **two** actual native calls on the same1NON/
2BE1 internal-loop tasks: seed7, one candidate, one cycle,200 updates, no MSA,
unchanged354/365-token contexts, zero sampling restarts and zero repair calls.
All four/six missing intervals are generated jointly within each Protenix call;
this is not cc94's independent-local method. Input source readiness, deposited
axes, water-only omitted context, native masks and callback conformance are already
retained. These are prospective development controls, not a final cohort, the
ten-entry tranche, general model qualification or authority for domain-sized gaps.

`scripts/run_retained_loop_context_sampling.py` provides a new private context
adapter only: selected observations are projected during sampling; preserved-only
observations remain immutable at final materialization. It uses the existing
pinned v1 model/checkpoint, actual MPS, declared NumPy3×3 alignment and the existing
serial pressure-v2 RSS8GiB/MPS6GiB/1800s guards. No production/frozen runner or
residue limit is edited. Original input requests are unchanged.

Attempt01 failed **before worker/model loading** because the upstream checkout's
concrete `scripts` package shadowed the project namespace on PYTHONPATH. It records
zero model calls and no prediction directory;2BE1 did not launch. The original
entry script and failure are retained. The separately declared attempt02 removes
only upstream paths from startup PYTHONPATH (the worker adds them after loading),
binds the failed attempt's receipts and permits recovery only for this exact proven
zero-model startup failure. Seed7 has not been resampled. Attempt02 is launched;
actual model updates and outcomes still require completed receipts.

An independent raw-endpoint auditor is prepared. It checks complete source-record
custody, atom ownership,200 reinjection updates, actual model MPS calls, geometry
and global L before any repair is declared. It explicitly retains the limitation
that the legacy private runner does not separately dump the full historical sampled
atom-identity axis. Clean input features and endpoint geometry are not a continuous
trajectory/topology guarantee.

### Attempt02 import failure diagnosed; established MPS LayerNorm restored

Attempt02 failed in6.665867s (peak aggregate RSS476135424 bytes), without a
resource stop. Its retained stack reaches `runner.inference` import, then the CUDA
`fast_layer_norm_cuda_v2` import/JIT path and missing Ninja. This occurs before
runner construction, checkpoint loading and the later `seed_everything` call.
`device-final.json` records zero model MPS calls; no predictions directory exists,
and2BE1 was not launched. This is not a sampled failure or model-size measurement.

The existing successful Mac development/fixed-validation launchers explicitly set
`LAYERNORM_TYPE=torch`; the new private wrapper had omitted it. Preserve attempt02
and its declaration-bound entry script, then restore exactly this established
environment setting. No Ninja/CUDA installation, upstream patch, CPU fallback,
timeout increase or scientific policy change is involved. A launcher regression
asserts torch LayerNorm even if the parent requests CUDA, no generic fallback and
project-only startup PYTHONPATH. Four focused adapter/readback tests pass; Ruff
passes. Attempt03 is separately declared and launched with the same two actual
model-call budget. Recovery admission accepts only the two demonstrated pre-runner
import failures, zero model calls, absent prediction directory and no second launch.
Actual attempt03 execution/seed consumption remains pending its final receipts.

The canonical plan and durable status now reflect the completed correctness queue,
194 independent source readbacks, current181-review scope and synthetic/native
preflight passes. Historical source-series records remain preserved; the aggregate
ledger and current-scope overlay, not old per-series snapshots, define current totals.

### Attempt03 consumed1NON seed7; retained MPS refusal, no resampling

The restored environment loaded the pinned checkpoint and executed actual MPS
model tensors. `device-limit-exceeded.json` records179234 module hook events,
including200 `DiffusionModule` forwards. Peak MPS driver allocation was6740672512
bytes, exceeding the unchanged6442450944-byte (6GiB) limit. The wrapper self-sent
SIGTERM. No candidate/native arrays or completed reinjection trace were written;
200 module forwards do not independently prove200 completed callback projections.
No scientific geometry outcome can be assigned.1NON seed7 is consumed and must not
be sampled again under the two-call declaration.2BE1 had not launched.

The parent concurrently attempted `killpg` and got `PermissionError`, losing its
usual execution/host telemetry receipt. A separate process inventory found no
remaining exact child command. Preserve this missing-receipt limitation rather
than synthesize a successful guard receipt. The declaration-bound guard/entry
scripts are snapshotted beside the original plan; independent
`resource-failure-readback.json` hashes the device/log/plan/block evidence and
records the resource refusal. The model hook event count is not a count of
independent proposals or diffusion updates.

The demonstrated guard race is fixed: after group-signal permission denial, wait
at most10s for confirmed child exit before retaining the error in the normal
report. A still-live unmanageable child raises and forbids continuation. Normal
SIGTERM/SIGKILL cleanup remains bounded; all resource thresholds stay unchanged.
Three negative/cleanup tests pass (seven focused adapter/readback/guard tests total).

Attempt04 is a **one-call first-execution continuation for2BE1 only**, not recovery
of1NON. Its declaration binds the retained failure/readback and refuses if2BE1
already has a workdir or the original1NON child is alive. The original two actual
calls, seed7/200 updates/one cycle, no repair/restarts and RSS8GiB/MPS6GiB/1800s
remain unchanged. It is launched; its actual outcome awaits final receipts.

The expanded retained-data correctness set now passes82 focused tests, including
the three guard cleanup regressions. A fresh complete fast-lane rerun is declared
under `final-correctness-fast-attempt01/`; it blocks on the shared executor lock
before starting, avoiding competition with2BE1. It has the same1800s ceiling,
retained full log and process-group timeout cleanup. No new scientific call is
introduced by that test run. Its outcome is pending, not inferred from launch.

### Attempt04 resource refusal; original two-call tranche fully accounted

2BE1 executed its first and only declared seed7 call. The guard retained returncode
−15, `mps-driver-limit`,127.915571s, peak aggregate RSS2551414784 bytes and peak MPS
driver6874939392 bytes against the unchanged6442450944-byte limit. The child again
self-stopped and the group signal returned permission denied; this time the fixed
guard confirmed/reaped the exited child and retained the full execution/telemetry
receipt rather than crashing. The module trace includes200 diffusion forwards,
but no completed callback trace, native arrays or candidate endpoint was written.
No generic CPU fallback occurred. Swap growth remains telemetry, not a substituted
stop reason or justification for raising a guard.

`retained-loop-context-sampling-final-readback-attempt01/readback.json` independently
accounts for both actual calls and binds225 evidence files. Attempt03's two later
edited executables are checked against byte/hash-identical executed snapshots;
this substitution is explicit, not pretending current files match the old plan.
Both native configs retain200 updates, one candidate and one cycle. The results
are **2 resource refusals,0 endpoints,0 assessable raw geometry outcomes**—not0/2
geometry passes and not a measured universal model/input-size limit.1NON's missing
parent receipt remains visible;2BE1's complete resource receipt is retained.

The two actual calls are exhausted. No new seed, repeat, repair, smaller-context
resampling, alternative model or higher-resource run is declared. Native data and
synthetic callback success remain separate positive transport evidence; they do
not rescue the absent raw endpoints. Recruitment and input-budget calibration
remain paused as directed. The full final fast lane is still running; the overall
scientific handoff remains incomplete rather than being marked closed by these
technical controls.

### Final fast attempt01 retained environmental failure; PATH-only rerun

The standalone final test launcher omitted the selected environment's `bin` from
PATH (the earlier serial queue correctly prepended it). Attempt01 therefore reports
six `FileNotFoundError: dvbfixer` subprocess failures,2550 passed, three optional
skips and33 slow deselected in299.82s. Preserve the failed full log/exit receipt;
this is not a green suite and does not justify modifying test assertions or skips.
All six exact failing node IDs pass with the interpreter bin restored to PATH.
`final-correctness-fast-attempt02/plan.json` binds the failed log/result and declares
a full rerun with only that child-environment correction, unchanged1800s budget and
shared lock. No scientific sampler is involved. Attempt02 is launched, pending its
completion receipt. Agent documentation maps and license inventory checks also pass
(the inventory check is not legal clearance).

### Final fast attempt02: isolated inherited-pipe timing failure retained

With corrected PATH, attempt02 reports2555 passed, one failed, three optional
skips and33 slow deselected in304.44s. The remaining failure is the existing
`test_runner_preflight_timeout_cleans_up_inherited_output_pipes`: it returns
`runner-timeout` rather than success. No production preflight code changed in this
continuation. Do not report this full run as passed.

A separately retained no-model diagnostic (`preflight-cleanup-diagnostic-attempt01/`)
wraps process wait/group signalling without changing the tested timeout/assertions.
It passes: the parent responds in~0.399s against its existing0.7s effective wait
(0.2s request plus existing0.5s startup grace), then `killpg(SIGKILL)` returns
success to close inherited sleeper pipes. This establishes a working isolated path,
not the exact cause of the earlier full-suite timeout. A single isolated pass does
not erase the failed full result.

The pipe-cleanup test fixture now uses a minimal `/bin/sh` parent, background
`sleep10` retaining stdout/stderr and a builtin `printf` of the same protocol JSON,
instead of spending the short budget loading Python's subprocess module. The
scenario and both original assertions remain: timeout0.2s, elapsed<2s, passed
response. Production preflight, its startup grace/drain bounds, protocol validation,
separate timeout tests and all scientific guards are unchanged. All17 preflight
tests pass. This is a fixture-isolation change, not an increased timeout or broad skip.
`final-correctness-fast-attempt03/` binds the failed prior log and the changed test
file and launches another full fast lane under the unchanged1800s guard. Outcome
is pending; original failures remain retained.

### Final fast attempt03 passed; bounded tranche settled, scientific gates open

`final-correctness-fast-attempt03/result.json` records exit0 and the retained log
digest: **2556 passed, three optional Torch/Biotite skips,33 slow deselected**, two
SWIG deprecation warnings,295.44s. No test was skipped to hide a failure, timeout
raised or production preflight changed. Both earlier final failed logs remain
available; the fixture-only repair does not establish the exact cause of the
historical full-suite timing failure. Ruff, agent maps, release version sync,
generated CLI reference and license-inventory checks pass; legal clearance is not
inferred. These are technical correctness checks, not scientific qualification.

All current background jobs are settled. The source inventory remains233 tasks,
194 independently checked ready projections,31 refusals and eight retained time
blocks. Current loop scope is12 small controls, one explicit7N8X exclusion and181
pending reviews, not194 admitted tests. The two actual native calls remain resource
refusals without endpoints; their zero assessable geometry outcomes must not be
reported as geometry failures or successes. No declared model calls remain.

The canonical plan/status preserve unfinished cohort/split, endpoint transport,
matched comparisons/conditioning, topology and manual donor/public-integration
gates. A later model-executing stage requires a distinct declaration and budget,
not an implicit retry or increased resource envelope. Recruitment and limit
calibration remain paused. This closes the **bounded technical tranche only**;
WP1–WP7 are not marked complete and the overall scientific handoff is not closed.

### Renewed autonomous continuation: passive sampler-return study, fresh seed11

The user again explicitly requested continuing until the whole handoff is complete.
The prior two-call budget was only a development tranche, not a user-imposed stop
on the broader plan. Preserve its immutable225-file outcome readback and boundary
receipt; live status and both receipt-bound documents are snapshotted with exact
hashes before continuing. No historical refusal becomes a pass.

The pinned model source calls `sample_diffusion`, then distogram and confidence
heads. Old instrumentation observed200 diffusion forwards and a later driver limit,
but did not retain the sampler return. The new **passive** private interception
saves that actual return and passes the same tensor object onward. It never
substitutes coordinates, reruns a step, bypasses heads or changes weights/config.
The snapshot commits all payloads before `capture.json`: full native atom identity
axis (chain case/insertion codes), names/elements, full coordinate array,200 callback
records and source-restored raw PDB. Native fixed observations and complete source
custody receive independent readback; raw geometry failures remain failures. A
later resource refusal still prevents full-inference success even if raw survives.

Three synthetic no-model tests cover tensor identity, callback/method restoration
after downstream failure, rejection of incomplete traces, exact axis serialization,
payload hashes and exclusive output. Ruff passes. This tests transport code, not
checkpoint conformance.

`sampler-return-study-seed11-attempt01/plan.json` prospectively binds a **new**
two-call development control:1NON/2BE1 seed11, one candidate, one cycle,200 updates,
unchanged354/365-token contexts and RSS8GiB/MPS6GiB/1800s pressure-v2. Both seed7
requests remain untouched; separate derivative requests differ only in `seeds` and
bind their original lineage. No source recruitment, size experiment, MSA, repair,
asset download or public/frozen runner edit occurs. Resource failures are retained
and permit only the other prospectively declared first execution; technical errors
stop the rest, with no implicit retry. The study is launched. Its completion and
independent raw audit remain pending; no sampling success is inferred from launch.

### Retained topology-frame readback prepared, declaration blockers retained

While the seed11 study runs, a separate no-model/no-optimizer script prepares a
readback over all12 historical1HE8 arms. It uses N/CA/C segments **only from actual
serialized force-field bonds**, not residue-number/proximity interpolation, and
retains implicated atom identities, source/generated ownership, phase force groups/
parameters and callback metadata. Both source/proposal/final and all retained
callback/returned states are in its proposed scope. The existing experimental0.75Å
segment diagnostic stays nonfatal and cannot certify no passage between sparse
frames or native knot validity. Two synthetic controls pass: absent/case-distinct
bonds cannot become inferred bridges, and an exact generated/observed intersection
is reported without a topological-defect or publication-gate claim.

Its declaration attempt01 found the old normalized workspace path unavailable;
the historical mobility worker explicitly reads `source.pdb` beside `request.json`,
and that source is checked against the frozen digest instead. Attempt02 then hit
the outer120s command limit before writing `plan.json`; neither attempt executed
the topology readback. Preserve both owned directories; do not call them completed
or silently rerun under a larger budget. All retained frame files remain untouched.
The script and focused tests are available for a separately accounted continuation;
there is no new optimization and no change to the12 existing scientific outcomes.

### Seed11 native study completed: two raw returns independently materialized

Both new calls retain full-inference MPS-driver refusals (211.071114s and788.042486s)
under unchanged limits; neither is a completed native success. Passive capture
nevertheless committed both actual sampler returns before subsequent heads. The
independent audit under `sampler-return-study-seed11-readback-attempt01/` verifies
payload hashes, full case-sensitive native atom axes,200 ordered zero-post-overwrite
callbacks, represented observations directly in the native arrays, independently
reproduced materialization, full source records and exact complete endpoint masks.
Both are globally all-L and preserve observed records byte-for-byte.

**Raw geometry passes0/2**.1NON retains junction C–N distances up to3.9705Å, angle/
amide defects and clashes;2BE1 retains junction distances up to5.2609Å, angle/amide
defects and clashes. These are now actual assessable raw failures, unlike seed7's
absent endpoints. Raw/native/restored bytes and full-inference resource failures
remain retained; no successful geometry or native execution is inferred from the
transport success. No new donor approval or public publication occurred.

### Bounded joint CPU repair prospectively declared and launched

`sampler-return-seed11-repair-attempt01/plan.json` binds both audited saved raw
proposals and declares two CPU jobs, with no model calls or solvent. Each first
keeps every observed atom fixed and jointly moves **all generated atoms in all
gaps**. Only a remaining geometry failure (never an observation/admission failure)
permits one subsequent shared sequence±2 observed-anchor stage. Both stages retain
the original source references; stage2 starts from the saved generated-only result,
not a new NN proposal. There are no optimizer restarts, at most1000 full-phase
iterations per stage, and the same serial1800s/job/RSS8GiB/host-pressure guard.

Stage2 reuses the already tested0.05Å inward finite-wall buffer but independently
validates the unchanged outer backbone/Cα/heavy0.5/0.5/1.25Å limits. No gate widening,
error-only atom mask, source ratcheting or fixed-observation relabeling occurs.
The necessity-gate test passes; raw failures remain raw failures even if repair
later succeeds. The two jobs are launched; their actual CPU outcomes and separate
independent final readback remain pending.

### Retained topology subset and native covalent transport checks queued

The full12-arm topology declaration attempt03 again hit the same outer120s limit,
after retaining eight per-run inventory records but before its final bound plan.
No topology analysis executed from either partial all-arm declaration. A **separate
narrow task**, `1HE8-accepted-endpoint-segment-frames-attempt01/`, declares the two
previously joint-accepted endpoints and all32 saved callback/returned frames each,
plus source/proposal/final. It does not replace the uncompleted12-arm scope or
qualify its other ten arms. Its declaration completed under the unchanged120s
outer limit; analysis has its originally declared1800s ceiling and blocking shared
executor lock. It is queued behind the CPU repair. No model or optimizer runs in
this diagnostic, and clean sparse frames cannot prove continuous no-passage.

`retained-native-peptide-transport-attempt01/` separately declares two pinned native
data-only checks, queued on the shared lock. Existing template/identity preflight
is reused with checkpoint/network/model forward forbidden before upstream imports.
The added readback inspects the actual atom-array bond list and native chain labels
for every intended source C–N junction around generated gaps. It preserves the full
native graph/identity axis and reports missing bonds/cross-native-chain seams as
measured transport facts, without changing encoding, repairing sources or asserting
causality of the raw geometric defects. A synthetic disconnected/connected graph
control passes; source mapping alone is not mistaken for covalent continuity.
Actual native results are pending; no absent-junction finding is assumed in advance.

### Private review export prepared, without selecting or approving donors

While the guarded CPU jobs and queued readbacks execute, a separate exporter is
prepared for the saved raw/repaired natural-gap stages. It requires the independent
repair audit and its digests, copies source/raw/all available repaired stages into
one source frame, and emits a manual PyMOL overlay plus explicit pass/fail metadata.
It does not invent a complete reference or align away movement. Every candidate
remains unselected/unapproved; the all-gap joint native proposal is one lineage,
not independent per-gap donors. Full native resource refusal and raw failure labels
remain visible even if a repaired stage passes.

The export uses canonical `DatRecord.save/load` and independently checks that its
sidecars contain exactly `request.generated_atoms`, never fixed observations,
prepared-source additions or released anchors. A case-distinct/insertion-coded
roundtrip test passes. All nine new focused tests (capture, topology helper,
repair necessity/readback, native graph and sidecar ownership) pass. The export has
not yet run because its independent repair-audit prerequisite is still pending;
this is private workflow preparation, not donor approval or public integration.

### Shared CPU repair and two data-only diagnostics completed

Both seed11 CPU jobs completed and were independently re-assessed:
**0/2 accepted**. Generated-only repair is insufficient in both cases. Sequence±2
release removes1NON's amide/clash findings but leaves D74 Ramachandran outlier
(φ48.892°,ψ−7.693°); there are no observation/movement violations.2BE1 retains
bond/angle/amide/clash findings and five outer movement violations after release.
Both endpoints remain globally L and preserve records outside their declared
mobile masks. The independent audit explicitly uses saved PDB precision because
solver coordinates/trajectory were not separately retained in these jobs.

The native covalent-graph preflight confirms **all20/20 intended C–N junctions**
(1NON8/8,2BE112/12) in the actual native atom-array graph, with no cross-native-chain
junctions. The proposed missing-native-bond explanation is **not supported**.
Feature/identity correctness still does not imply a usable sampled loop.

The narrow1HE8 topology readback completed on70 source/proposal/final/saved frames,
3391 actual covalent backbone segments per complete endpoint, zero segment alarms.
No minimizer was restarted. This leaves sparse unsaved trajectories, sidechain
ring threading, native knots and the other ten arms unqualified; it is not a fatal
topological acceptance policy.

`sampler-return-seed11-private-review-attempt01/` was exported from the independent
failed-repair audit. It includes source/raw/both repair stages, manual overlays,
unapproved recipes and six canonical generated-only `.dat` sidecars. All manual
approvals/selections remain unset; raw/final failures are plainly labelled.

### Residual atom mapping: distortion appears during repair, no causal claim

`sampler-return-seed11-residual-readback-attempt01/` checks six saved endpoints.
1NON A74/D74 start with negativeφ; different generated-only/released endpoints
introduce or retain the reported positiveφ defects. This bounds when a defect
appeared at saved-stage resolution, not the exact optimizer transition or cause.
2BE1 B/PHE256 CA–CB goes1.535726→3.725882→3.566750Å across raw/generated-only/
released stages. The raw phenyl cycle is planar (max deviation0.003846Å), with no
piercing among the bounded tested paths. Both repaired cycles are **undefined**
for the declared0.10Å planar-surface diagnostic (deviations0.359220/0.313468Å),
not clean no-threading results. Graphs here are the explicit OpenMM endpoint PDB
standard-residue topology, not an unsaved historical temporary FF graph. No force-
balance/local-trapping/collision cause is asserted from those measurements.

The renewed complete fast lane passed **2565 tests**, three optional skips,
33 slow deselected in306.99s. This predates the next two MODELLER transport tests,
which pass separately.

### MODELLER comparator: multichain transport failure retained and preflighted

A new source-matched two-case MODELLER comparator was declared using the existing
one-initial/one-fast-loop protocol, all gaps selected jointly, no hidden reference,
native chirality sweep, failed-loop fallback or source-restoration repair. It uses
the full observed template (716/686 target residues), unlike Protenix's local
354/365-token contexts. It is not a matched-context or matched-RNG-seed superiority
comparison and does not calibrate a residue limit.

Attempt01 entered1NON `model.make` but failed its alignment check: the inherited
single-chain `FIRST:@/LAST:@` header read only161 residues versus644 in the4-chain
alignment. No native proposal was written, and2BE1 was not launched. Retain the
failure/readback/executed source; seed7 is not replayed or declared unconsumed.
A private explicit first/last-chain PIR header then passes actual no-generation
MODELLER alignment checks for both sources (644/610 observed,716/686 target).
Case-sensitive endpoint-chain regression passes. No shared/frozen helper is edited.

Attempt02 separately declared fresh MODELLER seed42 for both cases, binding the
native PIR preflight and old failure, unchanged1800s/RSS8GiB guards. Both new calls
completed. Total `model.make` attempts remain three including the retained failed
seed7 call, not two uniformly successful attempts.

### Independent MODELLER audit caught restoration adapter error; no resampling

The first audit refused worker `source-restored.pdb`: passing the full native axis
to the canonical materializer updates represented fixed coordinates too. The new
worker had labelled an un-restored structure as restored. Retain those files and
the audit failure; neither is accepted or silently overwritten.

The explicit **transport-only** recovery independently re-normalizes the unchanged
native output, recomputes the observed-only global fit, and materializes generated
coordinates only while copying all original observations. A drifted-observation
regression passes. Corrected derivatives and readback are separate under
`retained-loop-modeller-controls-readback-attempt02/`; no MODELLER call is repeated.
Native2BE1 geometry passes, native1NON does not; both are L. Observed drift after
global observed fit is RMS/max1.473957/18.452627Å and1.042662/9.779388Å respectively.
Correct source restoration preserves all observed records but passes geometry0/2.
The invalid worker restoration is explicitly recorded, not reclassified as an
original successful stage. Later import-only formatting keeps an exact executed
entry snapshot so current source need not falsely match the old declaration hash.

### Same common repair declared for correctly restored MODELLER proposals

`retained-modeller-common-repair-attempt01/` prospectively declares two CPU jobs
on those corrected derivatives. The native seed remains42; **optimizer seed11**,
generated-only→necessity-gated sequence±2,1000 iterations, zero restarts/solvent,
original source references, inward0.05Å buffer and outer0.5/0.5/1.25Å limits match
the Protenix repair exactly. There are zero new proposal calls. The shared private
worker is parameterized only for input paths and native provenance; its executed
seed11 source is snapshotted with the old declared hash, and audits permit only
that exact retained snapshot, never an arbitrary changed implementation. The new
jobs are launched; final outcomes/readback remain pending.

### MODELLER common repair completed: one bounded full-task pass

The two jobs completed; independent `retained-modeller-common-repair-readback-attempt01/`
accepts **2BE1** after the necessity-gated sequence±2 stage and refuses1NON.
All six2BE1 gaps are repaired together; generated-only leaves one amide defect,
and the anchor stage removes all reported geometry/admission findings, preserves
global L and remains inside all original outer movement bounds. Nonmobile source
records and observed metadata are exact.1NON retains bond/angle/amide/clash findings
and one movement violation. This is one case/seed/workflow success, not general
MODELLER superiority, native-source preservation, uncertainty calibration or public
qualification. Protenix and MODELLER context/seeds differ as already declared.
No raw, native-execution or original restoration failure is retroactively passed.

`retained-modeller-private-review-attempt01/` contains both cases, including failures,
and generated-only canonical sidecars. Native MODELLER seed42 and CPU optimizer
seed11 are distinct metadata fields; recipes select/approve nothing. The exporter
now validates the exact audit-owned repair root/files before copying, uses uniquely
backend/seed-labelled overlay objects and explicit local file paths, and distinguishes
full-inference resource refusal from completed classical proposals. The earlier
Protenix exporter source is retained with its recorded hash; its old review is not
overwritten or silently relabelled. A native-vs-optimizer seed regression passes.

### Prospective native MODELLER fixation control, fresh seed101

The native `AutoModel.select_atoms` implementation explicitly documents selection
of atoms to optimize and defaults to all atoms. The new private arm delegates this
initial-optimization selection to the same existing gap-only `select_loop_atoms`.
Every invocation must map to **exactly all generated atom identities** in the full
target axis; missing atoms, extra observations, duplicate IDs or empty fallback
fail closed. The final native drift is still measured: supported selection does
not itself prove no observation movement during earlier initialization.

`retained-modeller-fixation-seed101-attempt01/` declares one new candidate per case,
fresh native seed101 (not a replay of7/42), one initial and one fast loop, same
whole-source template/PIR and unchanged1800s/RSS8GiB/pressure guard. No repair,
neural call, post-hoc native coordinate clamp, hidden reference or public capability
is introduced. Source restoration now correctly passes generated coordinates only;
the earlier invalid full-axis restoration is retained as an explicit adapter error.
Three MODELLER transport/mask regressions pass. The two new calls are queued/launched;
their native mask receipts, observed drift and independent geometry readback are
pending and not assumed successful.

### Settled seed101 results and all-input consolidation

Both native seed101 MODELLER calls completed. Independent readback verified the
exact generated-only initial-selection masks, but **did not qualify observation
fixation**: observed-fit RMS/max drift was 0.291952/2.418977 Å for 1NON and
0.369726/2.468794 Å for 2BE1. Both corrected source-restored endpoints failed
geometry while preserving observations and global L. The same common CPU repair
then completed both stages for both cases and independently accepted **0/2**.
No thresholds were relaxed and no consumed seeds were replayed.

Twenty saved MODELLER stage measurements localize drift after the source-exact
`target.ini`. For 1NON, initial-to-loop observed coordinates are identical. For
seed101 2BE1, the 16 changed observed labels form an exact within-residue coordinate
permutation preserving the incident covalent graph. No renaming was applied and
source-record custody was not reclassified. The earlier initialization cause
remains unproven; denied vendor-source reads were not bypassed.

The private library now contains two entries × three arms × three stages:
**18 stage records, not 18 independent structures**, one geometry-passing stage
and zero manual donor approvals. It binds the same source/sequence/masks across
arms, canonical generated-only `.dat` sidecars, filters and manual overlays.
Native context and seed differences prevent a backend-superiority conclusion.

Read-only review covered **194 retained structures and 1280 declared gaps**:
181 require loop/construct/context review, 12 are small controls, and 7N8X remains
explicitly excluded. There are 112 context-omission review flags, 95 prepared-source
provenance flags, four incomplete-family flags and 26 known development-component
exposures. These flags overlap and are not automatic admissions or refusals.
Older inspected controls absent from the census make the exposure overlay a lower
bound, not proof that other components are untouched. No length-based admission,
new recruitment, replacement quotas or final split was introduced.

The [consolidated execution ledger](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/ledger-attempt01/execution-ledger.json)
and [case-arm table](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/ledger-attempt01/case-arm-outcomes.tsv)
verify 1253 directly bound payloads. All 233 source tasks are accounted for:
194 ready, 31 refused and eight time-blocked. Recent execution comprises four
actual Protenix calls (four full-native resource refusals, two assessable raw
endpoints), five MODELLER `model.make` attempts (one alignment failure, four
completed proposals), six CPU repair jobs and 12 stages. One of the six case-arms
passed final acceptance: MODELLER seed42 2BE1 after common repair with optimizer
seed11. Historical results and adapter failures remain unchanged.

This supersedes the earlier trailing claim that no checkpoint-backed proposals
executed: they did execute, with outcomes distinguished above. All calculations
in this series are now settled. Following the user's prioritization correction,
no further 1NON/2BE1 restoration runs are queued. Current work is reproducibility,
retained-input review and documentation. Fixed known sequences are unchanged;
restoration concerns missing coordinates, not new protein sequences. Final cohort,
group-safe split, broad scientific qualification and public integration remain open.

The independent [all-input handoff receipt](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/review-handoff-attempt01/receipt.json)
rechecks every per-entry report against its request/source, case-sensitive sequence
hashes, declared gap intervals and human-readable TSV. All 194 entries/1280 gaps
agree. A [blank adjudication worksheet](../../tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/review-handoff-attempt01/adjudication-worksheet.tsv)
retains controls and exclusions, groups rows by existing conservative component,
and leaves every scientific decision pending/unassigned. It introduces no admission
rule or outcome/length ranking. Future decisions require a separately versioned
review record; do not edit the bound facts or treat blank reviewer fields as approval.

### Consolidation software verification

The fresh complete fast suite passed **2573 tests**, three optional Torch/Biotite
skips, 33 slow deselected, two SWIG deprecation warnings in 300.37 seconds. The
selected interpreter's `bin` was prepended to `PATH`; the full log is retained in
`tmp/mac-real-gap-continuation-2026-11-12/handoff-consolidation/software-verification/pytest.log`.
The 47 focused retained-data/guard tests also pass. Ruff, local file-link checks,
agent documentation maps, CLI reference and license-inventory currency pass;
inventory currency is not legal clearance. Earlier failed/passing runs are retained.
These software results close the consolidation verification, not scientific WPs,
manual donor approval or final cohort admission.

### 2026-11-13: concrete context and construct reviews across all retained inputs

The [actual-context audit](../../tmp/mac-real-gap-continuation-2026-11-12/retained-context-omission-audit-attempt01/context-audit.json)
reads the original full assembly and original protein projection for each of the
194 ready inputs. Every retained protein observation has the same identity and
coordinates; no ATOM record was silently omitted. Actual omitted HETATM records
match the previous declarations exactly. There are 48 entries with no omitted
atoms, 34 with HOH/WAT-only omissions and 112 with other omitted components.
Reports identify exact omitted residue instances and nearest original observed
anchor atoms. These distances neither estimate contacts with unknown coordinates
nor establish that distant context is irrelevant. Prepared-source coordinates are
not substituted for the original context frame.

The [construct audit](../../tmp/mac-real-gap-continuation-2026-11-12/retained-construct-annotation-audit-attempt01/construct-audit.json)
rechecks all 1280 requested gaps against retained deposited absence axes, including
every compatible source-label candidate for assembly copies. All agree. Of these,
811 gaps in 182 entries are terminal. Deposited entity-fragment annotations are
available in 52 entries, database alignments in 194 and sequence-difference
annotations in 128. Raw annotation ranges are retained, not silently interpreted
as author numbering or proof of a local-loop task.

The [difference audit](../../tmp/mac-real-gap-continuation-2026-11-12/retained-construct-difference-audit-attempt02/difference-audit.json)
uses the dictionary-defined `struct_ref_seq_dif.seq_num` pointer to
`entity_poly_seq.num`, independently checked against the deposited sequence scheme,
entity, label chain, case-sensitive author chain and monomer. It never substitutes
an author residue number for a label index. See the
[wwPDB item definition](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_ref_seq_dif.seq_num.html).
No additional experimental structure or sequence was downloaded.

Explicit depositor descriptions `expression tag` or `cloning artifact` overlap
**229 gaps in 106 entries**: 196 gaps have partial coverage, and 33 gaps in 23
entries have complete coverage. Coverage requires agreement across all compatible
source-label candidates. The remaining 1051 gaps lack mapped artifact annotations;
this is not proof of native-sequence identity or absence of tags. No motif-based
tag inference, mask trimming or automatic exclusion is applied.

The first difference report retained 56 records without a coordinate-label match.
Inspection showed that all 56 are explicit database-relative deletions: deposited
monomer/sequence position are absent, while database monomer/position are present.
The second report classifies them as not applicable to a deposited coordinate
target, rather than inventing residues or hiding them as mapping failures. They
occur in five entries; zero other annotation records remain unresolved under the
declared exact mapping. Counts include assembly-candidate reuse and are not numbers
of independent structures. The first report and its executed source remain retained.

The [case-review table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-case-admission-reviews-attempt02/case-review.tsv)
and [bound per-case review data](../../tmp/mac-real-gap-continuation-2026-11-12/retained-case-admission-reviews-attempt02/case-reviews.json)
combine these findings for all 194 inputs. Three factual subchecks are verified for
every entry: deposited gap axes, original projection observation custody and actual
context omissions. Concrete remaining questions include artifact-task semantics,
terminal-construct interpretation, complete context, prepared-source ownership,
family grouping and previous development exposure. The 181 pending active-scope
reviews, 12 controls and explicit exclusion remain separate. No scientific
admission, final split, new restoration calculation or change to known sequences
follows automatically from these factual checks.

### 2026-11-13: gap-level scaffold and construct separation

Autonomous continuation resumes the current retained-input review, not the settled
two-case sampling chain. `scripts/build_retained_gap_reviews.py` joins the bound
case reviews, independently checked source reports, deposited construct/difference
reports and unchanged requests. It validates exact case-sensitive gap axes and
the observed/absent partition of each target chain before publishing a fresh table.
Duplicate identities, inconsistent boundaries, changed bindings or denominators
fail closed; an existing output directory is never overwritten.

The [gap-level table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-reviews-attempt01/gap-review.tsv)
contains all 194 entries/1280 gaps: 469 internal, 419 N-terminal and 392 C-terminal.
All 229 explicitly annotated artifact overlaps are terminal: 153 N-terminal and
43 C-terminal gaps have partial coverage, while 22 N-terminal and 11 C-terminal
gaps have complete coverage. The other 1051 gaps lack mapped artifact annotations,
not necessarily construct changes. No sequence-motif inference or mask trimming
was introduced.

There are 120 entries containing internal gaps: 12 internal-only and 108 mixed;
74 entries are terminal-only. Within the 181 active-scope pending reviews these
counts are ten internal-only, 101 mixed and 70 terminal-only. Controls and the
explicit exclusion remain retained. Shortcuts that admit internal-only cases or
drop the terminal gaps of mixed entries would be new scientific decisions, not
consequences of this table. Contiguous observed run lengths are descriptive
sequence-axis evidence, not spatial/domain boundaries or calibrated uncertainty.

The [independent stdlib readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-reviews-attempt01/readback-receipt.json)
verified 973 bound payloads, all TSV cells, all request axes, observed runs and
annotation counts/coverage without importing the builder. All scientific admission
and split fields remain unset. The 39 focused gap/case/construct/context tests
pass, including changed-evidence refusal and output no-overwrite controls. Initial
test collection with the unrelated base Python failed because the package was not
installed there; the established scientific interpreter passes. No environment was
modified. No new source acquisition, sampling, refinement or public integration
was performed. The subsequent complete fast suite passed **2610 tests**, three
optional Torch/Biotite skips, 33 slow deselected and two SWIG warnings in 296.12s.
The full log is retained at
`tmp/mac-real-gap-continuation-2026-11-12/retained-gap-reviews-attempt01/pytest-fast.log`.
Ruff, agent documentation maps, CLI reference currency, license inventory currency
and all 27 local file links in the two updated state documents also pass. These
checks establish software consistency, not legal clearance or scientific admission.

### 2026-11-13: deposited assembly membership and explicit connection custody

`scripts/audit_retained_context_connections.py` performs a read-only check of all
194 retained ready sources. It separates deposited polymer labels, declared
assembly1 membership and actual atom pairs in the original full/projection PDBs.
It preserves raw connection types/symmetry fields and never infers covalency from
distance. A sequence-compatible source-chain candidate is not evidence that the
specific deposited chain or its partners belong to the selected assembly.

The [connection audit](../../tmp/mac-real-gap-continuation-2026-11-12/retained-context-connections-attempt02/connection-audit.json)
finds three source protein/nonprotein annotations: one metal coordination in
3U8X and two in 3UQ5. All are outside the chosen assembly1 member-label set.
6ZH0's source interchain disulfide is between two polymer chains, not a
protein/nonprotein connection merely because one chain is outside assembly1.
Label co-membership alone still does not demonstrate actual operator-pair transport.

Canonical PDB parsing resolves 370 unique full-context explicit pairs: 315 between
retained atoms and 55 between omitted context atoms. There are zero resolved pairs
crossing the projection boundary. All 315 retained pairs in 54 entries survive
unchanged by exact atom identity; no new resolved projection pair appears. The
original full coordinates and every raw connection record remain retained. The
first audit is preserved separately; attempt02 adds the retained-pair custody check.

The [separate readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-context-connections-attempt02/readback-receipt.json)
checks 972 bound inputs, 194 reports, all TSV cells, exact pair sets and independently
recomputed deposited endpoint memberships without importing the research auditor.
Its graph readback deliberately reuses the canonical PDB reader; it is not an
independent graph parser or proof that missing annotations are absent interactions.
Among 371 raw records (316 SSBOND, 55 LINK), one resolves to no atom pair: 5XZ0's
SSBOND connects A/CYS93 to its crystallographic image. 5XZ0 remains a small control,
not an active primary-cohort entry. The deposited annotation
also names the same residue with distinct symmetry operations, while assembly1
uses an identity assembly operator. No image was synthesized, no self-bond added
and no correspondence between differing CIF/PDB symmetry-number spellings assumed.
This is a concrete pending crystal-context decision, not a silently cleared case.

The 49 focused connection/context/gap/admission tests pass, including polymer-vs-
assembly classification, case/insertion-code identity, stale CONECT retention,
lost/new retained-pair detection, serial-independent identity and changed-evidence
refusal. No scientific admission, mask or sequence change, new source, model call,
refinement or public capability follows. Complete nonbonded context/force-field
treatment and the existing 112 nonwater-context reviews remain open.

The subsequent full fast suite passed **2627 tests**, three optional Torch/Biotite
skips, 33 slow deselected and two SWIG warnings in 298.53s. The log is retained at
`tmp/mac-real-gap-continuation-2026-11-12/retained-context-connections-attempt02/pytest-fast.log`.
Ruff, agent documentation maps, CLI reference currency, license inventory currency
and all 31 local file links in the two state documents also pass. No scientific
admission or qualification is inferred from this software verification.

### 2026-11-13: actual gap lengths and first focused-submodel scope decisions

The user's concern about enormous gaps prompted separating lengths from counts.
The1280 earlier gap instances are summed across194 experimental entries, including
assembly chain copies; they are not1280 residues, one enormous single structure
or1280 independent biological examples. `scripts/summarize_retained_gap_lengths.py`
now binds exact source requests and reports scope-separated distributions and a
[per-entry length table](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/entry-gap-lengths.tsv).
The retained inventory includes36616 missing residue instances. Its876 descriptive
within-entry exact-sequence/interval classes do not establish independence or
identical context; copies remain separate records and original reconstruction tasks.

The original181-pending inventory contains455 internal gaps (median9, maximum141),
407 N-terminal gaps (median12, maximum672) and381 C-terminal gaps (median13,
maximum1069). Thus the size concern is real, not merely a misunderstanding of
the aggregate counts. Controls and the explicitly excluded7N8X/812-residue span
are kept separate. No universal maximum length or revised quota was introduced.

Four individual whole-entry requests are now deferred under the already approved
focused-submodel amendment, with [versioned source-only reviews](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/case-scope-reviews.json):

| Entry | Deposited study and reason the original whole-entry task is deferred |
|---|---|
|5MSC|A-domain structure; original target also asks for a535-residue terminal extension outside that observed domain.|
|7QL0|Explicit filament fragment266/297–391; six copies each have69 observed residues of441, with308 N-terminal and64 C-terminal residues absent.|
|8YNM|Death-effector-domain assembly; three chain copies each require297 residues beyond their observed182-residue segment.|
|8ZJK|Explicit core assembly; two ELMO1 copies require535 N-terminal residues outside the observed core.|

These are case-specific scope judgments from retained deposited titles and exact
observed/absent axes, not automatic decisions from length alone. They defer the
original complete requests, not the source structures permanently. No gap was
trimmed, sequence changed or new narrower task created, and no generated outcome
was used to select these decisions. Source readiness remains a different status.

The [effective review overlay](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/effective-scope-review.tsv)
now contains177 pending entries, four deferred focused-submodel entries,12 controls
and the explicit user exclusion. No final cohort or split is frozen. For the177
remaining pending entries, the descriptive distribution is:

| Gap kind | Instances including copies | Median residues | Range |1–8|9–24|25–49|50–99|100+|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|Internal|450|9|1–141|211|156|61|14|8|
|N-terminal|393|12|1–672|167|106|65|32|23|
|C-terminal|363|12|1–1069|143|94|41|57|28|

The remaining large spans are unadmitted review cases, not queued loop-modeling
jobs. Length bins are reporting labels only. Old scope tables are preserved;
the new overlay explains rather than silently rewrites their181-pending denominator.

The [independent stdlib readback](../../tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/readback-receipt.json)
checks all194 entries/1280 instances and399 bound inputs, recomputes distributions
and sequence-axis equivalences, and checks custody of the four decisions without
importing either reporting/review materializer. This is factual verification,
not a second independent scientific adjudication. The44 focused length/gap tests
pass. New sources, sampler/optimizer calls, public capabilities, request-mask
changes and final scientific admissions remain zero.

The subsequent full fast suite passed **2648 tests**, three optional Torch/Biotite
skips,33 slow deselected and two SWIG warnings in303.89s. The retained log is
`tmp/mac-real-gap-continuation-2026-11-12/retained-gap-lengths-attempt01/pytest-fast.log`.
This closes the bounded software verification, not the remaining scientific reviews.

### 2026-11-13: all retained1HE8 segment frames checked without new restoration

The user requested describing the task as restoration of missing coordinates at
known sequences. No sequence-design task or text prompt for a new protein was
introduced. Frozen technical ownership fields retain their existing names for
compatibility. This continuation sends no new request to a structural model.

`scripts/continue_retained_1he8_segment_frames.py` selects the exact complement of
the two previously checked arms within the immutable twelve-arm history. It leaves
the old auditor unchanged, verifies prior payloads and source/request custody,
and uses the existing serial pressure-v2 CPU executor:1800s/arm,8GiB process RSS,
10GiB output cap. The longer600s outer manifest-preparation allowance addresses
the historical120s orchestration interruption, not a larger molecular input or
weaker execution/scientific limit. Seed0 is an output reservation, not sampling.

All [ten declared readback jobs](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/completed.json)
completed in7.24–8.28s each, with peak measured process-tree RSS377–428MiB.
They screened367 optimizer frames plus30 source/proposal/final stage instances.
Each complete temporary topology has3391 backbone segments from actual serialized
harmonic covalent bonds; no numbering/proximity-based peptide bridge was invented.
All397 new instances yielded zero segment alarms at the unchanged0.75Å diagnostic
distance and1e-5Å interior-intersection tolerance.

The [combined twelve-arm table](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/combined-arm-summary.tsv)
now covers467 instances:431 optimizer frames and36 source/proposal/final instances.
These are saved-frame/stage counts within twelve arms of one structure, not467
independent structures or continuously observed trajectories. A
[separate readback](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/readback-receipt.json)
reconstructs the graph, checks exact particle identities/units and callback metadata,
and recomputes all467 screens without importing either audit worker. It shares
the canonical segment primitive, so this is independent graph/frame recomputation,
not an independent geometry algorithm or knot/ring classifier.

The ten historical arm rejections remain unchanged, as do the two prior combined
passes. Zero segment alarms do not resolve ring/surface warnings, close an open
path, prove absence of passage between saved frames, or qualify whole-entry cofactor
context. No structural-model call, optimizer call, new source, observed-coordinate
change, policy relaxation, donor approval or public capability was introduced.
The16 focused continuation/frame tests pass. An initial synthetic test incorrectly
expected an unsorted bond tuple; the fixture now uses canonical identity ordering,
without changing the actual topology reader or diagnostic method.

The complete fast suite passed **2662 tests**, three optional Torch/Biotite skips,
33 slow deselected and two SWIG warnings in308.63s. The retained log is
`tmp/mac-real-gap-continuation-2026-11-12/1HE8-remaining-segment-frames-attempt01/pytest-fast.log`.
The combined focused segment/continuation set passes23 tests, including the
existing rigid-transform, projected-crossing, bonded-neighbor and budget controls.
Ruff, agent documentation maps, CLI reference currency, license inventory currency
and46 local file links across the two state documents also pass. Software success
does not qualify rings, unsaved trajectories, physical accuracy or public support.

### 2026-11-13: focused TYR972 ring across saved1HE8 frames and historical controls

`scripts/audit_retained_1he8_ring_frames.py` extends the nonfatal readback to one
already implicated covalent ring, without new coordinate restoration or optimizer
calls. The ordered TYR972 cycle and SER1044 CA–CB witness identities come from the
retained historical investigation, not a new residue template or proximity graph.
Every ring edge and witness bond must exist in the supplied covalent topology.

All [16 declared CPU readback jobs](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/completed.json)
completed under the existing serial pressure-v2 bounds: twelve saved arms plus
four retrospective control structures. Saved frames use actual serialized harmonic
covalent bonds restricted to input heavy atoms. Original incomplete sources are
explicitly partial: available-coordinate subgraphs are labelled and missing ring
or witness observations are not silently treated as clean findings. Complete
frames fail closed if required topology coordinates are absent. The0.10Å ring
planarity scope and existing acceptance/movement policies remain unchanged.

The [arm summary](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/ring-frame-summary.tsv)
contains467 frame/stage instances from the twelve arms:455 defined ring surfaces
with zero piercings, plus12 original observed-source instances where the ring is
absent. There are no nonplanar undefined surfaces among these455 assessable arm
states. This does not mean every aromatic ring or all unsaved motion was tested.

Four controls reproduce the historical distinction:

| Control | Focused ring result | SER1044 CA–CB length |
|---|---|---:|
|Selected raw assembly|Undefined: nonplanar beyond the declared surface scope|1.535Å|
|Retained seed7 endpoint|Defined, no piercing|1.532Å|
|Earlier rejected seed42 endpoint|Defined, witness bond pierces the surface|3.039Å|
|Earlier rejected seed101 endpoint|Defined, witness bond pierces the surface|3.041Å|

The positive endpoints are from different historical experiments, not newly failed
members of the twelve-arm cohort or newly repaired structures. Their geometry
rejections remain intact. Control topology is labelled OpenMM canonical endpoint
PDB topology, not a recovered historical temporary force graph. Controls omit a
restoration ownership mask; raw probe ownership strings therefore have no
provenance interpretation. The [per-frame table](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/ring-frame-findings.tsv)
marks control ownership unqualified rather than relabelling restored atoms as
observations.

The [separate readback](../../tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/readback-receipt.json)
recomputes every arm/control probe and checks exact frame coverage, missing-source
status, counts and historical controls without importing the ring audit worker.
It shares the canonical graph helpers/ring primitive; this is not an independent
topology algorithm, a causal timing result or a proof that no path passage occurred.
Nine new wrapper tests plus existing ring/continuation controls pass28 focused
tests. Source/request masks, sequences, manual donor approvals, scientific
admission, movement limits and public capabilities are unchanged.

The full fast suite passed **2671 tests**, three optional Torch/Biotite skips,
33 slow deselected and two SWIG warnings in304.57s. The complete log is retained at
`tmp/mac-real-gap-continuation-2026-11-12/1HE8-ring-frames-attempt01/pytest-fast.log`.
These are software regression results, not broader ring/trajectory qualification.

### 2026-11-13: deposited construct targets are not full-native restoration requests

Source review returns to the existing development pool, without restarting the
settled1NON/2BE1 series or making new structural-model requests.
`scripts/review_retained_construct_targets.py` performs three explicitly bounded
case reviews:1AHW,1JPS and1KTZ. It checks each original request against every
compatible deposited source-label candidate, including exact canonical construct
sequence, single-coordinate-model heavy-observation indices and all absent-coordinate
intervals. It does not run another alignment, infer missingness from database
numbering or expand a deposited construct to a full native reference protein.

The [case reviews](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/construct-target-reviews.json)
support the source-level interpretation that these requests concern coordinate
restoration within the deposited experimental construct. This is distinct from
the earlier four deferred requests extending far beyond the observed focused core.
A title containing “domain” alone neither admits nor excludes an entry.

| Entry | Target/observed residue instances | Internal gap lengths | Terminal gap lengths | Omitted context |
|---|---:|---|---|---|
|1AHW|1294/1256|8,8|3,8,3,8|None|
|1JPS|658/626|4,7|1,8,4,8|Water only; treatment pending|
|1KTZ|468/376|18,18|12,10,6,12,10,6|Water only; treatment pending|

Counts include assembly copies, not independent biological replicates. Together
these are20 gap instances in three already-inspected entries, not a newly admitted
three-entry evaluation cohort. Every terminal interval stays in its original
request; short terminal completion is not silently called two-sided loop repair.
The absence of mapped artifact annotations establishes neither native-sequence
identity nor absence of construct changes; the deposited sequence is the exact target.
Database-only deletion annotations stay outside that coordinate target.

The [human-readable summary](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/construct-target-summary.tsv)
keeps remaining gates explicit. The retained preparation scope refuses terminal
gaps for all three cases; this historical profile result is preserved, not treated
as a newly tested admission verdict for every backend. Water treatment/complete-
context endpoint assessment, backend-specific full-request admission and a frozen
prospective method/context/seed/resource manifest remain prerequisites for a new
execution. No public profile was modified or bypassed.

The [separate readback](../../tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/readback-receipt.json)
recomputes all candidate source sequences/absence axes,20 original gap instances
and source-label ambiguity without importing the reviewer. It checks factual
custody, not a second independent scientific adjudication. All three source inputs
are unchanged preparation originals and already belong to known development-exposed
components; no untouched final component or split is inferred. The effective scope
stays177 pending, four deferred focused submodels,12 controls and one user exclusion.
The23 focused construct/annotation/difference/admission tests pass; structural-
model calls, optimizer calls, new sources, masks and final scientific admissions
remain zero for this continuation.

The full fast suite passed **2683 tests**, three optional Torch/Biotite skips,
33 slow deselected and two SWIG warnings in301.67s. The retained log is
`tmp/mac-real-gap-continuation-2026-11-12/construct-target-reviews-attempt01/pytest-fast.log`.
Ruff, agent documentation maps, CLI reference currency, license inventory currency
and58 local file links in the two state documents also pass. These are software
checks, not backend admission, legal clearance or final scientific qualification.

### 2026-11-13: weightless full-target RFD3 input check on1KTZ

`scripts/preflight_retained_rfd3_full_target.py` now materializes a separate private
conditioning artifact for the exact retained1KTZ request:468 residues/four chains,
all eight gaps,2986 observed atoms and710 atoms needing coordinates. It includes
the entire target, not a local crop or substituted sequence. Canonical native
conditioning represents2984 observed atoms plus710 zero-coordinate placeholders;
the two already-observed true-terminal OXT atoms remain preserved-only observations
outside the native axis, explicitly not native constraints or new coordinates.
Original18 SSBOND annotations are retained in the conditioning artifact, and any
CONECT serial update uses the shared canonical helper rather than guessed serials.
The original source, request, full water context and metadata remain bound.

Attempt01 failed before upstream input construction because the prototype eagerly
imported the parent scientific executor, pulling Gemmi into an isolated RFD3
environment that does not contain it. Its [failure log](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt01/execution/1KTZ-rfd3-native-data-seed-0/output.log)
and exact source snapshot remain intact. Attempt02 lazily imports that dependency
only in controller mode and hardens duplicate-source-serial refusal. Conditioning
bytes and manifest are identical between attempts; input specifications differ
only by their owned artifact path. No package was installed or engine patched.

The [native stage summary](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/native-stage-summary.tsv)
comes from actual `DesignInputSpecification.build()` and
`prepare_pipeline_input_from_atom_array()` execution under the existing serial
pressure-v2 CPU guard. The data job finished in8.29s with measured process-tree
RSS about1.91GiB. Owned engine revision/source diff/hashes and the interpreter are
bound. Torch/JIT checkpoint loading, module forwards and Python socket connections
are blocked before upstream imports; offline settings are enabled. Safetensors is
not installed in this isolated runtime, and its guard is explicitly reported false
rather than claimed installed. No checkpoint, model forward, MPS qualification,
sampling, input-size calibration or scientific acceptance occurred.

| Native stage | Atoms/fixed/unfixed | Peptide edges present | Junction edges present | Declared disulfide edges present |
|---|---:|---:|---:|---:|
|Builder (intermediate)|3694/2984/710|0/464|0/10|0/18|
|Pre-model pipeline|3694/2984/710|464/464|10/10|0/18|

The missing builder peptide graph is intermediate: its completion in the next
stage is measured, not hidden as464 broken output bonds. Both stages preserve all
468 source components/four chains, every fixed sequence/atom mask and relative
observation coordinates after a common COM translation; maximum residual is
7.97e-6Å under the declared1e-4Å float32 representation tolerance. Missing-coordinate
placeholders remain exactly zero, not reconstructed positions. Native arrays and
per-stage diagnostics are retained before subsequent validation can fail.

The [disulfide readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/native-disulfide-readback.tsv)
shows all18 observed SG–SG endpoint pairs retain their fixed coordinates and source
distances, despite absent explicit edges in these native arrays. Thus adjacency
loss is not physical disulfide rupture, a model-output defect or a proven final
feature-consumer defect. It does leave the declared strict graph-transport check
unresolved/failed. No distance-based rebonding, mask relaxation or model execution
was used to force a pass. Further investigation must establish the actual downstream
bond-feature contract before interpreting the finding or declaring readiness.

The [separate array readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/readback-receipt.json)
recomputes identity/mask/sequence/coordinate and bond counts without importing the
worker or RFD3. Recorded terminal-profile refusals, pending water treatment,
source-scope decisions and all177 pending complete reviews remain unchanged.
The full fast suite passes2700 tests, three optional Torch/Biotite skips,33 slow
deselected and two SWIG warnings in296.93s;21 focused transport/conditioning tests
pass. The full log is retained at
`tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-full-target-data-attempt02/pytest-fast.log`.

### 2026-11-13: isolated AF3 token-bond semantics, not full native featurization

Read-through of the pinned RFD3 pipeline finds `AddAF3TokenBondFeatures()` after
virtual-atom padding. Installed AtomWorks `ml/transforms/bonds.py:246–271` reduces
atom adjacency to token adjacency, then explicitly zeros pairs for which both
token-start atoms have `atomize=False`. Its documented scope is polymer–ligand
and ligand–ligand bonds, not the complete covalent graph. RFD3's feature aggregation
can additionally clear bonds involving unindexed tokens; the ordinary and batched
encoders consume `token_bonds`. These are source observations, not executed neural
paths. Covalent-modification annotation selects polymer/non-polymer pairs, not
ordinary protein–protein links; standard amino acids are excluded from default
atomization. Actual full-input atomization and other conditioning channels have
not been measured here.

A separately declared [five-control native transform audit](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-token-bond-semantics-attempt01/result.json)
executed in3.05s on five two-node synthetic arrays. A polymer pair with an explicit
edge has raw token adjacency but zero AF3 bond features; the same pair without an
edge produces the identical feature matrix. Positive polymer–atomized and
atomized–atomized controls retain their explicit edge; an unbonded atomized pair
remains zero. All matrices are symmetric2×2 with zero diagonals. The
[independent NPZ readback](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-token-bond-semantics-attempt01/readback-receipt.json)
recomputes these results without importing the worker, AtomWorks or RFD3.

Sources, interpreter command, declared timeout, arrays and guards are retained
under `tmp/mac-real-gap-continuation-2026-11-12/rfd3-token-bond-semantics-attempt01/`.
Checkpoint/JIT loading, Torch module forwards and socket connections were blocked;
no RFD3 engine, real structure, request or environment was modified. These are
generic feature semantics, not a completed1KTZ feature pass or evidence about
restored coordinates. The prior strict declared graph-transport failure remains
unchanged. In particular, absence of an S–S bit in this intentionally restricted
channel neither establishes physical rupture nor qualifies source-link custody or
coordinate fixation through sampling. Next distinguish actual native token axes
and fixed-coordinate conditioning from explicit topology custody, without patching
the engine, inferring bonds by proximity or loading a model to force a pass.

### 2026-11-13: native source-config prefix through feature aggregation

The separate declaration in
`tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-feature-prefix-attempt01/`
executes the real native inference transform prefix on unchanged complete1KTZ
conditioning. It resolves the pinned source-default validation transform with
batch size1 and no residue cache, **not checkpoint-authenticated hyperparameters**.
No checkpoint metadata is read. The prefix ends after
`AggregateFeaturesLikeAF3WithoutMSA`, before `BatchStructuresForDiffusionNoising`,
EDM noise and random coordinate augmentation. Module-forward/checkpoint/socket
guards remain active. The CPU job completed in7.74s at about1.60GiB process-tree RSS.

The [six-stage readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-feature-prefix-attempt01/stage-summary.tsv)
tracks input, native atomization, virtual-atom handling, design-reference features,
token bonds and final feature aggregation. All stages retain3694 canonical atoms;
the final axes are468 non-atomized sequence-fixed tokens/four chains with no
virtual atoms and no unindexed tokens. Fixed known sequence prevents virtual
padding here. Dense association naming replaces1368 displayed atom names with
internal names such as `V0`, but `gt_atom_name` retains original identities: these
internal aliases are not new virtual atoms or missing observations.

All2984 represented observations preserve exact identity/masks and coordinates
after the existing common centering translation (maximum residual7.97e-6Å under
the unchanged1e-4Å tolerance). All710 missing-coordinate atoms remain zero. Final
`motif_pos` matches fixed native coordinates exactly and masks unfixed coordinates
to zero. `ref_pos` is all zero on this source-config protein-only route, not the
observation channel. All468 restype one-hot rows decode to the known residues;
`atom_to_token_map` matches the source-bound axes. The468×468 `token_bonds` matrix
is all zero as expected for this non-atomized protein input. All18 declared S–S
edges remain absent in the actual atom graph; the prior strict failure is retained.

The [independent verifier](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-feature-prefix-attempt01/readback-receipt.json)
imports neither the native worker nor engine and refuses six deliberate in-memory
mutations: fixed mask, atom-token map, observed motif coordinate, unfixed motif
leakage, virtual flag and encoded sequence. Original saved arrays are not mutated.
RFD3 encoder source consumes motif coordinates through a fixed-pair distance
embedding; this is a source observation, not a neural-forward or fixation test.
No sampler, noise, restored coordinate, optimizer, engine patch or input/profile
change occurred. Feature transport does not establish sampling-time hard fixation,
the output coordinate frame, checkpoint-specific behavior or whole-context
acceptance. Those remain separate from source-link custody, preserved-only OXT,
pending water treatment and the unchanged scientific admission gates.

Read-only next-boundary inspection finds that ordinary sampler reinsertion and
final frame alignment are conditional on `allow_realignment`; its source default
is false. `insert_motif_at_end=True` alone is therefore not a hard output guarantee.
The saved `coord_to_be_noised` annotation matches native coordinates, but no sampler
path was executed. Preserve these source observations separately from effective
checkpoint settings and actual fixation evidence; alternate/batched paths remain
to be reviewed before declaring the output contract.

### 2026-11-13: autonomous layered fixation and stored-configuration audit

The separately declared `rfd3-fixation-contract-attempt01` executes eight synthetic
numerical-driver controls: ordinary/padded paths × two realignment flags × two
end-insertion flags. Inputs are six invented Cartesian points, not a protein; the
denoiser is a deliberately non-native analytic displacement, not an `nn.Module`.
The driver permits deformation when that arbitrary denoiser changes fixed points;
end reinsertion restores the test shape only under both flags. The
[independent saved-point audit](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-fixation-contract-attempt01/readback-receipt.json)
recomputes movement and pair distances without Torch. The later interpretation
record explicitly limits this counterexample to the numerical driver, not native
RFD3 as a whole. No production flag was changed to force a pass.

Continuing into the actual denoiser finds the essential additional mechanism:
ordinary `RFD3_diffusion_module.py` multiplies per-atom time by the complement of
the fixed mask; the padded denoiser does the same under valid-atom masks. Both
pass those times into `scale_positions_out`. For `f_pred=edm` and finite residual
updates, time zero makes that formula return its input coordinate exactly. This
is a mathematical identity, not a learned preference or final-realignment policy.
The separate [formula/composition audit](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-edm-fixation-formula-attempt01/readback-receipt.json)
checks six unbound scaling-method controls and six three-update analytic driver
compositions over EDM/noise-prediction/unconditioned modes and both tensor shapes.
Only the EDM compositions preserve fixed points bit-for-bit with realignment off.
Independent NumPy algebra reproduces the formula and update recurrences. There
are18 analytic non-neural calls, zero structural-model calls and zero real-input
sampling jobs. No module is initialized or forwarded. This finite-update evidence
does not assert that time-zero arithmetic sanitizes NaN/Inf, or qualify symmetry,
CFG, checkpoint-specific neural execution or actual1KTZ output.

The current local checkpoint's configuration was then read statically, not loaded.
`rfd3-checkpoint-metadata-attempt01` reads the ZIP directory and815458-byte
`data.pkl` member, parses167426 opcodes, and represents GLOBAL/REDUCE/NEWOBJ/BUILD
and persistent references as inert local records. It never invokes their
callables, unpickles, imports Torch/RFD3/OmegaConf or reads tensor archive members.
The symbolic parser fails closed on unsupported opcodes;19 controls include an
inert malicious-callable example and unsupported-protocol refusals. Literal
OmegaConf content, not neighbouring string matches, establishes the stored paths:
`train_cfg.model.net.diffusion_module.f_pred=edm`, sigma16, and sampler
`allow_realignment=False`. Stored end-insertion/jitter/fix-fraction fields are
absent and are reported as absent, not reconstructed as checkpoint settings.
Metadata hash, archive stat and selected values are retained in the
[static verification](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-checkpoint-metadata-attempt01/static-verification.json).
Tensor weights, historical-weight identity and arbitrary runtime overrides are
not authenticated by this metadata-only read.

The first stored validation transform differs from the earlier source-config
prefix in three fields: `provide_reference_conformer_when_unmasked=False`,
`max_binder_length=170`, and an original Linux residue-cache path. A separately
declared [stored-config prefix](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-rfd3-checkpoint-prefix-attempt01/readback-receipt.json)
runs the unchanged full1KTZ input using those values, overriding only batch size1.
The absent cache path is not replaced or created; no cache residues are loaded.
Checkpoint/JIT/pickle loading, module forwards and sockets are blocked, and the
prefix stops before noising. All96 saved annotation/feature arrays are identical
to the earlier source-config prefix. The independent auditor is reused with its
previous six mutation controls, not presented as a second independent algorithm.
The known18-link graph failure, OXT custody, sequence/observation masks and all
admission/profile/context gates remain unchanged.

A separate [retained-result readback](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-retained-fixation-readback-attempt01/readback-receipt.json)
uses only the historical1UBQA1–72/gap10–13 source, raw CIF and saved source-restored
PDB; no held-out coordinates or new materialization are used. It recomputes the
proper rigid fit and reproduces native observation RMS3.36e-14Å/max4.25e-14Å for546
observations; their source-restored records are byte-identical. Metadata records
the old legacy execution path but not effective prediction-mode hyperparameters;
do not infer them from missing fields or qualify1KTZ from this single smoke.
The first readback script failed at insertion-code byte indexing (`int.strip`);
its exact source and error explanation are retained. Explicit technical
continuation replaces integer indexing with a one-byte slice, without changing
inputs, criteria or outcomes. Finally, trainer output source explicitly sets
`atom_array.bonds=None` before formatting; that is not an executed final-output
custody check or proof of physical bond rupture. Original explicit source links
still require independently retained final-materialization evidence.

### 2026-11-13: approved development of a separate source-link custody plan

After the explicit decision boundary, the user approved developing the
[experimental RFD3 plan](../plans/rfd3-source-link-custody.md). This authorizes
planning, not a new adapter, terminal-profile bypass or model-backed run. The old
native18-link failure remains failed, separately from proposed final source-link
custody and eventual scientific acceptance. The plan preserves original request
identities, observed SG endpoints, OXT pass-through and all water/context gates.

Source review of `pdb_materialize.py` and `pdb_metadata.py` identifies a material
integration boundary: the general materializer can update fixed coordinates when
they are included on the candidate axis, and can drop invalid CONECT records.
The proposed wrapper must therefore refuse unresolved/invalid source links before
materialization, independently audit raw fidelity and frame transport, supply only
requested restored atoms, and independently compare final observation records and
resolved source links. Neither the helper name nor copied SSBOND count guarantees
the new contract. Raw fidelity, source-link custody and geometry remain distinct.

The plan specifies an offline negative-control matrix and staged gates, without
inventing a new sampler tolerance or reusing the prefix's float32 tolerance for
serialized output. No implementation, new model/optimizer execution, source
acquisition, seed reservation, public capability or scientific admission occurred.
Next is separately scoped offline adapter/readback qualification; terminal/water
decisions and a prospective precision/runtime declaration still precede any run.

### 2026-11-13: authorized offline custody implementation and independent controls

The user approved the private adapter/offline work and reiterated autonomous
continuation until a genuine decision is required. Four private modules now live
under `src/dvbfixer/model/diffusion/`: `rfd3_custody_records.py` binds original
source/context/request and strictly resolves canonical ownership and explicit
links; `rfd3_custody.py` checks raw sequence/masks/inventory and one proper global
frame before generated-only materialization; `rfd3_custody_validation.py` audits
final source bytes/links and exposes existing generated/junction/global-L geometry;
`rfd3_custody_geometry.py` adds declared fixed-endpoint length diagnostics using
canonical bond/disulfide metrics on diagnostic-only declared-pair topologies.
No source/native topology is patched and no pair is inferred by proximity.
Warnings remain canonical warnings in evidence but hold immutable-link review.
These are private pure/offline services,not a public runner or case admission.

The source ledger rejects invalid/dangling/duplicate link records,unsupported
image operators,source identity/element/occupancy mismatches and unsupported
partial-residue ownership. Raw audits require explicit precision,exact boolean
masks,complete canonical inventory/sequence and a nondegenerate proper fit;drift
is refused before restoring observations. Independent final readback re-resolves
identities and links and checks emitted missing-atom coordinates at PDB precision.
The module auditor openly reuses canonical PDB/link/frame primitives;the full1KTZ
artifact reader separately uses Gemmi/NumPy and its own bounded SSBOND-only decoder.
No report claims a second independent whole-topology algorithm.

Initial tests and the first SS-only1KTZ control passed. A supplemental invented
positive LINK/reciprocal-CONECT control then exposed a real integration defect:
per-record `_parse_explicit_links` needs coordinate lines to construct its serial
map,not only its parsed-atom argument. Exact pre-correction source/test snapshots,
the error note,initial control and initial fast-suite log remain retained in
`rfd3-custody-offline-attempt01/`. Explicit attempt02 supplies those coordinate
lines for CONECT and adds a positive reciprocal-record regression. No source
request,mask,geometry threshold or invalid-link rule was changed. Supplemental
valid LINK plus reciprocal-CONECT records now survive exactly,and three mutated
variants fail. A later test collection-only import spelling correction precedes
the final geometry-test runs;no model seed or scientific outcome was involved.

The [complete1KTZ control](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-custody-offline-attempt02/readback-receipt.json)
uses unchanged retained pre-model arrays with710zero placeholders. It writes
raw/source/context/request inputs before adapter calls under the existing guarded
serial CPU executor,with Torch/RFD3/safetensors imports and sockets prohibited.
All2,986 observations including OXT remain byte-identical,all18 source S–S pairs
resolve independently,and the common-frame residual is7.96e-6Å under the existing
prefix-only1e-4Å budget. The710 inserted entries are **not restored coordinates**;
the output filename says `NOT-A-PREDICTION`,whole-candidate geometry is unassessed
and scientific acceptance is false. Separate declared source-endpoint length
metrics emit zero findings;that does not qualify the placeholder structure.
The old native18-edge graph failure remains failed in separate evidence.

A separately declared [retained1UBQ replay](../../tmp/mac-real-gap-continuation-2026-11-12/rfd3-custody-geometry-replay-attempt01/readback-receipt.json)
uses the existing raw CIF and evaluation mask for1UBQA1–72/gap10–13,without any
held-out reference coordinates,new sampling or repair. Proposal masks are built
from the frozen evaluation request for this offline service test,not fabricated
as a native sampling trace. The new adapter reproduces the entire old PDB hash
`9fe73c9cf8fbd518c50aa0aa2e4b1b22326dfdcca41a9038a8525162e9a865ee`.
Independent final-byte readback and reused canonical geometry/global-L diagnostics
pass for546 source observations and28 old proposed atoms;combined declared-link
review also passes for this no-S–S historical control. This cropped protein-only
context does not decide whole1UBQ or1KTZ water policy or enlarge the test cohort.

The new tests include case-sensitive chains,insertion codes,OXT,ANISOU,common
frames,reflections/degeneracy,nonfinite values,mask/sequence/element changes,
raw/final drift,invalid annotations and declared-link distances. Four fixed-link
length/review controls demonstrate that exact source custody can coexist with
bad physical geometry;no source atom is moved to force a pass. There are41 custody
and four declared-geometry controls,64 focused tests with the existing19
materialization/metadata tests. Scientific execution,public profiles,frozen
requests,source recruitment,size calibration and terminal/water admission remain
unchanged. Offline-only precisions do not become future sampler defaults.
The final complete fast suite passes2,745 tests,three optional Torch/Biotite skips,
33 slow deselected and two SWIG warnings in301.67s. The earlier2,740/2,741 logs
remain retained at their actual implementation stages;they do not replace the
final corrected-code/declared-geometry result.

### 2026-11-13: authorized complete-1KTZ terminal pilot, retained resource stop

The user's `приступай` authorizes the separate case-specific experimental
full-1KTZ terminal study, conditional on source/water preflight and prospective
execution settings. It does not alter the old terminal-profile refusal, native
18-edge S–S graph failure or public/frozen request schemas. The request remains
468 residues, four chains and eight gap instances; no crop or new source is used.

The new private water checker uses canonical vdW radii, 0.6/0.9Å overlap thresholds
and the existing heavy H-bond envelope on noncovalent protein–water cross pairs.
It performs no bond inference/exclusions across sets and preserves original water
records, occupancy and altloc evidence. Eleven synthetic tests pass, including
very short contacts that the generic proximity-derived bond graph could hide.
The guarded source preflight runs without Torch/RFD3/safetensors imports, network,
model execution or repairs; source/request/full-context hashes bind unchanged
originals. Observed global L and all 18 declared source-link length metrics pass.

The [independent preflight receipt](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-preflight-attempt01/readback-receipt.json)
confirms 2,986 source observations and all 326 original oxygen waters, with unit
occupancies and no altlocs. There are **eight baseline warnings, zero errors**.
An initial independent comparison exposed a worker-summary defect: it compared
uppercase `Severity.value` to lowercase keys and reported zero warnings. The
detailed original findings already contained all eight warnings; independent
Gemmi/NumPy all-cross-pair recomputation agrees. Original worker, incorrect
summary and detailed findings remain retained; the corrected receipt explicitly
records the defect. The actual severe-contact gate uses enum comparisons and
was unaffected. No source atom, water, threshold or admission policy changed.

The [prospective execution declaration](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-run-attempt01/declaration.json)
reserves one seed101 under the previously explicit MPS seed-at-model amendment.
It binds the original owned checkpoint, pinned engine/source patches and 364
input/implementation paths; checkpoint SHA256 is
`9b3f85923e0d51e9453e15cdd2f8c666e7ce096a60577f57d11bbc54ae6d67c1`.
One sample/batch, 200 configured timesteps, two recycles, fast/Torch kernels,
compile disabled, workers0 and MPS fallback disabled remain unchanged. Actual
native atom identity/`gt_atom_name`, source components, masks, reference tensor
and lossless final float32 NPZ capture are declared before launch. Maximum
all-observation common-frame residual is1e-4Å before private generated-only
materialization; no reflection, per-chain frames, private observation reinjection,
refinement or tolerance escalation is allowed. Native sampler settings, including
end-of-sampler motif insertion, remain unchanged and are reported as native behavior.

Source conditioning is protein-only with 2,984 native fixed atoms and 710
missing-coordinate atoms. Two observed OXT remain preserved-only; all 326 waters
are immutable final-assessment obstacles, not native conditioning. No full-native
context claim is made. The full software suite passes **2,756 tests**, three
optional Torch/Biotite skips, 33 slow deselected and two SWIG warnings in300.92s;
75 focused tests, Ruff and agent-doc checks pass before the scientific attempt.

The only attempt stops after **19.20s** for `mps-driver-limit`, under the unchanged
serial pressure-v2 guard. Sampled peak MPS driver allocation is15,411,675,136 bytes
(14.35GiB), above the6GiB stop threshold. A periodic stop check is not a hard
allocation cap; retain the overshoot rather than claiming memory stayed below
6GiB. Peak aggregate RSS is3,153,903,616 bytes; swap growth5,847,449,600 bytes is
retained telemetry, not the declared stop reason. The process exits -15 and is
reaped; the guard records a subsequent group-signal permission error after
confirmed child exit. No automatic retry, different setting, limit increase,
source repair, refinement or public publication occurs.

The [independent partial readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-run-attempt01/readback-receipt.json)
confirms actual forward masks/identity for 3,694 atoms (2,984 fixed,710 missing),
known sequence throughout, and7.96e-6Å common-frame residual on the **input**.
Effective runtime metadata confirms float32 MPS, EDM,200 configured timesteps,
two recycles,`allow_realignment=false` and`insert_motif_at_end=true`. The explicit
seed/RNG amendment and input evidence survive. The device monitor observed58
MPS-bearing submodule outputs; these are not58samples, completed timesteps or a
completed seed. SIGTERM left no final native coordinate tensor or candidate PDB;
no raw-output fidelity, final geometry/global-L or restored-water verdict is
possible. The audit retains this as a resource-blocked attempt, **not** a geometry
rejection or a successful structure restoration. Further complete-target execution
needs a separate resource/host decision; no new structure or cropped task is
silently substituted. Broader cohort review and scientific qualification remain open.

### 2026-11-13: user correction — Mac memory engineering, not host transfer

The user clarified: the objective is being able to execute on Mac. The proposed
Linux/CUDA transfer was not the task. Continue engineering the complete retained
Mac request with the original limits, not shorter inputs, different checkpoints,
new structures, model-size sweeps or relaxed scientific gates.

Source inspection finds a large temporary allocation in native
`SinusoidalDistEmbed`: full pair distances are expanded to32 sine/cosine
frequencies before projection. A private query-row64 expression evaluator keeps
all3694x3694 dense output pairs but bounds those temporary arrays. It retains
native CPU frequency construction/transfer, distance norms, projections, masks
and float32. Native `low_memory_mode` is not used: its initializer path bypasses
token-feedback pooling and the optional atom transformer. Shared parameters alone
cannot establish equivalence to the normal path. No engine files are edited.

The [isolated controls](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-mac-memory-controls-attempt01)
retain each prospective plan, guard log and result. The first test command
dereferenced the venv Python symlink to the base interpreter, where Torch is
absent; it executed zero tests. A separate technical amendment preserves that
failure and uses the existing venv entry without dereferencing, with no install
or environment change. Seventeen CPU/MPS checks then pass. A fixed3694x3694x16
operator using random tensors, no checkpoint, returns every pair, all finite,
with a sampled native submatrix difference of0.0 and sampled driver peak
1,976,287,232bytes (1.84GiB). This is not a protein proposal or full-model proof.

The separately declared sinusoidal-row full-target arm retains the same seed101,
known sequence, all chains/gaps/source links/waters and budgets. Its complete
sinusoidal output is retained in telemetry, but a subsequent allocation triggers
the child MPS sentinel at6,538,625,024bytes (6.09GiB), after11.94s. The parent
records exit -15 with a null reason because the child exits between polling
iterations. Independent readback resolves the cause from the immutable child
sentinel without rewriting `execution.json`, and rechecks actual native input
identity/masks/common frame. There is no final tensor or candidate PDB.

Reference inverse-distance embedding and pointwise pair-MLP/projection operators
are then tiled with the same full pair axes and no skipped layers. Forty-two
CPU/MPS controls pass; a single fixed3694atom468token operator composition,
including unchanged native atom-to-token pooling, peaks at3,897,245,696bytes
(3.63GiB) and yields finite pooled values with matching sampled native embeddings.
This random-weight operator composition is explicitly not full initializer or
checkpoint qualification. The full-pair technical arm completes all four tiled
operations, native TokenInitializer and the first denoiser atom encoder, but
stops at7,658,504,192bytes (7.13GiB) after12.02s. Independent partial readback
passes actual input identity/masks/frame; no final coordinates exist.

Explicit synchronize/empty-cache hooks at selected native module boundaries
release only unused allocator cache, not live tensors, watermarks or model
features. Two additional controls verify bit-identical input/output, unchanged
parameter-key inventory and unchanged CPU/MPS RNG, and leave CPU execution
untouched. All44 isolated controls pass before the separately declared cache arm.
That arm completes three native diffusion-module calls and seven tokenencoder/
decoder returns; its next call remains partial. After21.79s the child sentinel
records6,601,523,200bytes (6.15GiB), again with parent null-reason/-15 race retained.
Boundary telemetry reaches `diffusion_token_encoder.pairformer_stack.0` before
the stop. First-noisy coordinates are saved immediately rather than relying on
a finally block that SIGTERM can bypass. No completed seed/step count or final
scientific verdict is inferred from submodule counts.

The next source-localized temporary is the pairwise gated native `Transition`:
last-feature RMSNorm, two expanded linear projections, SiLU(A)*B, then output
linear. Query-row tiling retains the complete feature normalization and gating
expression; no feature dimension, attention or reduction is split/approximated.
Twelve CPU/MPS dense-native-vs-row controls bring the isolated matrix to56passes.
The [Transition-row arm](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-mac-transition-attempt01/declaration.json)
is separately declared and launched with the same complete target, weights,
seed101,200 configured timesteps/two recycles and unchanged pressure-v2 caps.
It is a technical implementation arm, not an independent scientific replicate
or an unchanged retry. Actual completion/raw fidelity/final geometry/water
outcomes remain pending; helper controls are not a claim that full Mac execution
or structural restoration now works.

### 2026-11-13: complete Mac execution finally passes; candidate quality does not

The Transition-row arm completes on the same Mac/MPS under the unchanged
pressure-v2 budget: returncode0, no guard stop,443.090s wall time, sampled peak
MPS driver5,152,391,168bytes (4.80GiB), aggregate RSS2,906,308,608bytes. No CPU
fallback, watermark override, dtype change, source edit, chain/gap crop or resource
increase occurs. The full468-residue/four-chain/eight-gap request remains intact.
This is one successful complete-target engineering smoke, not a claim of broad
Mac scalability or scientific terminal-gap qualification. The normal200-point
schedule makes199 native diffusion-module calls and398 tokenencoder/decoder
returns with two recycles. Earlier same-seed technical arms are failures, not
additional completed scientific seeds or bit-identical baseline comparisons.

The [independent retained-byte readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-terminal-study-mac-transition-attempt01/readback-receipt.json)
rechecks prospective source/weights/implementation bindings, actual native
`src_component`/`gt_atom_name`, masks and known sequence. Final float32 NPZ exists;
raw common-frame maximum observed residual is7.9570714e-6Å, below the declared
1e-4Å limit, before private generated-only materialization. The worker's raw
audit agrees. All2,986 original observations remain byte-identical, all710
missing-coordinate atoms are present, all18 original S–S annotation records are
unchanged, and independently emitted missing-atom coordinates match the native
common frame at PDB precision. Shared canonical declared-link metrics emit no
findings and global L-chirality passes. The native18-edge graph failure is still
a separate historical failure; this final source-link pass does not rewrite it.

The candidate nonetheless fails geometry. Generated/junction assessment reports
one canonical bond-length finding,12steric findings,four peptide-junction reports
and aggregate bond/overlap gate reports; these categories are not all independent
defects. No admission-policy failure is needed to reject bad coordinates.
Independent all-water cross-set assessment agrees with the worker on13ERROR and
10WARNING contacts involving restored atoms, while the original source baseline
remains8WARNING/0ERROR. Water records/positions are unchanged and were assessment
context, not native conditioning. No proximity-derived covalent exclusion can
hide these contacts.

`candidate-UNQUALIFIED.pdb`, actual native arrays, source/context/request,
raw/custody/geometry/water audits, operator/cache telemetry and independent
receipt remain retained in the project. There is no public publication, automatic
loading, coordinate refinement, water movement, criteria change or further
sampling. The engineering objective **full execution on Mac within the original
caps** is demonstrated; acceptable missing-coordinate restoration is **not**.
Further quality/conditioning/refinement work requires its own bounded declaration
and cannot discard these rejected raw coordinates or their findings.

### 2026-11-13: fixed-context derivative blocked by unresolved disulfide target

The user renewed autonomous continuation after the Mac execution/candidate-quality
split. A separate prospective derivative study binds only the existing full1KTZ
candidate/source/context/request, with no new checkpoint sampling. Exactly710
missing-coordinate heavy atoms in92 residues may move;2,986 observed protein
atoms and326 water oxygens remain fixed. The existing v8 CPU/one-thread refiner
uses AMBERff19SB/TIP3P, NoCutoff and no solvent box, with geometry-first bonded/
softcore/full-FF phases capped100/75/75/250 iterations, zero restarts, seed101 for
FF preparation and a separate3Å missing-atom displacement bound. This is not
model CPU fallback, source repair or a rewrite of raw geometry acceptance.

The new private `rfd3_refinement.py` builds an FF-only input with original waters
and unique temporary serials, and audits a derivative without falsely demanding
refined missing atoms equal native output. All source and non-mobile records,
metadata, auxiliary/link records and physical layout remain unchanged; only the
declared missing-coordinate columns may differ. New FF hydrogens/temporary heavy
completions are not observations or publication atoms. Fifteen unit controls pass,
including an actual small-kernel water-O mass-zero/disposable-H control using an
existing tracked fixture. This is software evidence, not a new cohort structure.

The [guarded attempt](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-fixed-context-refinement-attempt01/readback-receipt.json)
finishes its retained failure report in5.74s. FF/H preparation had run, but the
worker stops before constructing the main simulation/minimizing missing heavy
coordinates: `CustodyError: temporary FF disulfide pairs differ from original
custody`. Its exact18-pair FF guard was too narrow. No refined PDB was produced,
no source/water coordinate was rewritten and no corrected retry was started.
The original plan/code/failure and rejected raw candidate remain retained.

An offline source-only readback finds19 OpenMM-inferred disulfides, not18: the
extra is observed A77–C77, already present at2.233Å in the immutable source.
Raw-candidate and source-only inferred pair sets are identical; no new bond was
invented from missing-coordinate atoms. Original ASU PDB/CIF declare only9
intra-chain pairs, rendering to18 copies in this biological assembly. Neither
original format declared the77interchain pair. Therefore do not diagnose this as
an assembly renderer silently dropping an original explicit annotation. Canonical
full-source diagnostics also flag the77pair's ~153° twist: both source length and
twist warnings are retained, not repaired by moving observations or dropped.

The original DBREF identifies mature TGF-beta3 A1–112 as UniProtP10600:301–412.
A separately declared single read-only [UniProt metadata receipt](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-disulfide-semantics-attempt01/P10600-receipt.json)
retains exact response bytes/HTTP headers/digest. This acquires biochemical
annotations, not another structure, model, checkpoint or benchmark case. Frozen
A/C sequences match the112-residue mature sequence exactly. UniProt disulfide
features307–316,315–378,344–409,348–411 and interchain377–377 cite experimental
PubMed8819159 evidence. Mapping through the exact DBREF supports a **native
oxidized hypothesis**:21 total pairs after retaining receptor/source pairs,
including A77–C77 and CYS7–CYS16 for each A/C copy. All42 CYS SG endpoints then
have one partner. Sequence match/native annotations do not prove the redox state
of the unresolved termini in this particular crystal.

The [independent saved-annotation/Gemmi readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-disulfide-semantics-attempt01/readback-receipt.json)
recomputes actual raw SG7–SG16 distances9.263Å and17.549Å. With SG16 immovable,
even reaching a deliberately generous2.5Å closing envelope requires SG7 movement
of at least6.763Å and15.049Å. It is geometrically impossible to close either
native terminal bond inside the declared3Å local-cleanup bound. This is a target
topology/terminal-reconstruction decision, not permission to tune minimization,
drop a disulfide, raise the old bound or modify source/water coordinates.

No native21-link hypothesis has been applied, no old request or annotation was
edited, and the source18-link custody pass, native18-edge graph failure and raw
geometry/water rejection remain at their original meanings. User confirmation
is needed before selecting the native oxidized target and a separately bounded
disulfide-aware restoration of the missing A/C termini. Other177 case reviews,
public integration and general qualification are not declared complete.
Software checks pass75 focused tests and the full fast suite:2,771 passed,
three optional Torch/Biotite skips,33 slow deselected,two SWIG warnings in297.22s.

Evidence housekeeping: the derivative readback's broad file inventory also
hashed a still-running, independent pytest log. Scientific inputs/results did
not change. Its historical517-byte prefix was recovered from the completed
append-only log by an exact original-digest match and archived explicitly as
`RECOVERED`; the original receipt and untruncated final log remain intact. The
[binding note](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-fixed-context-refinement-attempt01/live-software-log-binding-note.json)
distinguishes that administrative log progression from immutable scientific
evidence. No expected scientific hash or test verdict was rewritten.

### 2026-11-14: native oxidized target confirmed; topology passes, derivatives reject

The user answered “подтверждаю” to native oxidized TGF-beta3. Retain that
[authorization](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-topology-preflight-attempt01/authorization.json)
separately from the original18source annotations and frozen request. It resolves
the target-state decision; it does not prove the original crystal's redox state,
permit observed/water motion, rewrite old failures or enable a public capability.
Native21pairs, source-only inferred19pairs and explicit18pairs remain distinct.

The actual licensed MODELLER10.8 topology preflight completes in1.56s with no
resource stop. It loads full protein/water heavy topology, transports frozen
residue ordinals and applies21DISU patches. Actual native bonds contain exactly
21unique target pairs with42singly partnered CYS. The active native selection is
exactly202atoms/24missing A/C1–12residues. All3,494nonactive protein coordinates
and326water oxygens stay within7.97e-6Å; the written PDB has zero fixed-coordinate
residual. The [independent Gemmi/SSBOND readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-topology-preflight-attempt01/readback-receipt.json)
does not invoke MODELLER or optimize coordinates. This establishes native topology
and transport, not missing-coordinate geometry or repaired raw graph conditioning.

A separate prospective [terminal declaration](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-terminal-restoration-attempt01/declaration.json)
uses one deterministic native stereochemical procedure:250/750/1000CG-iteration
caps, all21DISU patches, canonical CA–N–C–CB chirality restraints, dynamic sphere
packing after the first bonded stage, no MD/randomization/neural sampling or
unrestricted initial automodel optimization. Only202terminal atoms are active.
It declares a new60Å active displacement ceiling from original raw coordinates,
not a changed/replayed3Å local-refinement experiment. Fixed native atoms include
the3,494other protein atoms,326waters and temporary terminal extras. Native stage
arrays/PDBs are saved and checked before final materialization; no posthoc clamp.

This procedure completes in4.69s; every stage has zero nonactive native motion.
Native global-L checks precede missing-atom-only materialization. Its
[independent readback](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-terminal-restoration-attempt01/readback-receipt.json)
confirms all2,986observed records,3,494inactive protein records and326water oxygens
unchanged,454L centres,21target bonds and native-output/final-coordinate agreement.
Maximum active motion is16.117Å. A/C SG7–SG16 distances improve from9.263/17.549Å
to2.064/2.321Å. All-generated water contacts fall13ERROR/10WARNING→8ERROR/7WARNING;
active-terminal contacts fall6ERROR/4WARNING→1ERROR/1WARNING. Baseline remains
0ERROR/8WARNING. The result is **not accepted**: terminal length/geometry and
whole generated/junction checks remain failures. Native fixation passes are not
scientific restoration passes.

One subsequent [fixed-context FF declaration](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-context-refinement-attempt01/declaration.json)
binds that exact rejected terminal parent, unchanged original source/fullwater/
request and approved21-pair graph. It is not a retry/rewrite of the historical
exact18guard failure. Existing v8 geometry conditioning runs100/75/75/250iteration
caps, AMBERff19SB/TIP3P NoCutoff, same-Mac CPU/Threads1, no solvent box, zero
minimization restarts. Exactly710missing heavy atoms are mobile and3,312observed
protein/water heavy atoms fixed. The3Å posthoc missing-atom movement bound is
explicitly relative to the terminal parent, **not** original RFD3 raw. Both
parent motion and raw lineage remain visible.

FF preparation confirms the exact approved21-pair topology,3,312fixed atoms,
326fixed oxygens and zero protein–water covalent edges. Execution completes in
45.20s with peak aggregate RSS204,587,008bytes, no guard stop/fallback or raised
limit. Independent readback confirms source records/global L and all21target
lengths below severe thresholds. All-generated water contacts improve to0ERROR/
3WARNING; the three warnings are terminal contacts, distinct from the8unchanged
source/water warnings. Yet **reject** the derivative:

- Fifteen missing atoms exceed the declared3Å movement scope; maximum5.237Å is
  B24CZ. This is a scientific custody/movement failure, not a technical crash.
  Preserve its failure report; do not raise the bound retrospectively or call
  the output an accepted parent for a later repair without retaining that failure.
- One non-planar C3–C4 amide (omega−113.4°) and11severe protein overlaps remain,
  concentrated at missing C3/C4/C7 versus fixed nearby residues. Two additional
  assessment issues summarize those failed gate categories, not extra atom pairs.
- Canonical S–S diagnostics still warn about A7–A16 twist~155°, C7–C16 length
  2.295Å and twist~59°. Fixed A77–C77 remains2.233Å/~153° with its two source
  warnings. “No severe length” is not canonical/full disulfide geometry acceptance.

The [FF independent receipt](../../tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-context-refinement-attempt01/readback-receipt.json)
binds full final bytes, kernel coordinates, actual FF mobility/graph and independent
Gemmi/NumPy chirality/water/movement readback; it does not rerun the optimizer.
All source/raw/native/derivative failures remain. No accepted complete restoration,
new structure, neural sample, changed source/water coordinate or public capability
is claimed. Next work is conformational escape in the missing termini and an
enforced movement policy, not blind minimization, outcome-driven bound changes or
another target-redox confirmation. No new user decision is presently established.

Software adds ten native-graph controls and ten pre-materialization transport/
terminal-derivative controls, including distinct uppercase/lowercase/insertion-code
identities, fixed-water drift rejection and protection of inactive missing atoms.
The focused set passes95tests. The full fast suite passes2,791tests, with three
optional Torch/Biotite skips,33slow deselected and two SWIG warnings in299.77s.
Its completed log is retained separately; independent scientific readbacks never
bind that log while it is still being written.
