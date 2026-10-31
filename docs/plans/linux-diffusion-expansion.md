# Linux Diffusion Expansion Implementation Plan

- Status (2026-10-20): partially implemented and experimentally evaluated; see
  the current evidence below and [plan index](README.md). Original work-package
  checklists are historical decomposition, not a claim that all their software is
  missing. New scientific scope and paused validation-policy work require separate approval.
- Target branch: `feature/dvbfixer-hardening`.
- Scope owner: Structure Preparation bounded context.
- Exclusive coordinating change group: `missing-atom-rebuild-chirality`.
- Production baseline retained: MODELLER remains the default and supported
  backend until a separate evidence-backed production decision changes it.
- Related plans and decisions:
  - [`diffusion-gap-reconstruction.md`](diffusion-gap-reconstruction.md)
  - [`diffusion-cli-integration.md`](diffusion-cli-integration.md)
  - [`diffusion-runner-environments.md`](diffusion-runner-environments.md)
  - [`diffusion-heterogen-gap-reconstruction.md`](diffusion-heterogen-gap-reconstruction.md)
  - [`small-diffusion-apple-silicon-experiments.md`](small-diffusion-apple-silicon-experiments.md)
  - [`../adr/0011-experimental-diffusion-cli.md`](../adr/0011-experimental-diffusion-cli.md)

## Outcomes

- [ ] Complete checkpoint-backed Linux/NVIDIA acceptance for
  `protenix-v1-cuda` in a pinned environment.
- [ ] Complete class-aware protein geometry validation before widening the
  public scientific scope.
- [ ] Add reference-dependent benchmark metrics without combining unrelated
  score scales into one synthetic score.
- [x] Replace the global 3-12 gap limit with profile-specific, evidence-backed
  capability strata.
- [ ] Complete causal acceptance and broader engine support for surrounding-aware
  gap reconstruction, building on the implemented Protpardelle local-partner
  conditioning slice.
- [x] Add backend-aware batch output and diffusion-enabled `zbs`, and expand the
  implemented one-chain mosaic-first `homology` slice without changing MODELLER
  defaults.
- [ ] Preserve fail-closed admission, independent validation, atomic bundle
  publication, exact identities, and zero D-C-alpha output throughout.

## Future scope: diffusion completion of partial residues (2026-10-29)

**Status: private masks/contracts, insertion, validation and generated-only
refinement implemented; narrow private cc94 sampling executed. Public integration
and broader scientific qualification remain open.** The user
requested that diffusion also reconstruct missing heavy atoms within otherwise
observed protein residues. This is distinct from missing whole residues, bounded
relaxation of existing source atoms, hydrogen placement, or template-based
PDBFixer repair. It is not a prerequisite for starting the existing Linux/CUDA
acceptance work and must not silently widen current admission.

The supplied 5J13 case motivates this work: besides its separate SO4 admission
issue, ordinary diagnose reports residues with incomplete heavy-atom sets.
Removing nonprotein molecules alone does not implement partial-residue diffusion.

### First implementation slice (2026-10-30, after `931f977`)

`model/diffusion/atom_completion.py` now inventories exact canonical missing atom
identities without modifying source bytes or an existing frozen request. Its
immutable masks keep observed atoms, whole-residue-gap atoms and partial-residue
atoms separately owned; chain case, insertion codes and target positions remain
exact. Independent audit recomputes the inventory from the digest-bound source.
Backbone defects are classified, not granted sampling capability. Unknown atoms,
nonterminal OXT, substitutions, alternate locations and contradictory ownership
fail rather than being silently normalized. Chemical atom names come from the
existing canonical catalog and sequence placement remains in the shared service.

The original supplied 5J13 contains **20 partial sidechains / 62 missing canonical
heavy atoms**. The inventory matches PDBFixer's independent missing-atom set
exactly after insertion-code whitespace normalization. The earlier prepared
protein-only projection retained **63 newly present heavy atoms**, including one
terminal OXT: raw prepare had 4,052 heavy atoms and the target-aware projection
4,050, versus 3,987 originally observed protein heavy atoms. Those preparation
counts must not be confused with 62 canonical atoms eligible for future
diffusion completion. The original sulfate-containing input is still not
admitted for sampling by this planning service.

The actual pinned cc94 engine reader at
`ee378400f25b801fa481028000f9060183d7fb4c` was exercised with the production motif
writer and a separate 1UBQ input containing both a whole-residue gap and a partial
LYS A27. Reader atom37 masks retained N/CA/CB/C/O as present and CG/CD/CE/NZ as
absent. Input/request/writer/reader digests and arrays are retained under
`tmp/atom-completion-931f977/`. This proves **input-mask transport only**: weights
were not sampled, and no structure was published as atom-completion output.

Verification: 34 atom-inventory regressions and a 100-test focused set passed;
the complete fast lane passed 1,840 tests. The portable evidence record is
[partial-residue masks and reader probe](../research/partial-residue-mask-foundation-2026-10-30.json).

- [x] Implement typed source-bound partial-residue mask inventory, exact ownership
  audit, gapless planning and real 5J13/PDBFixer inventory regression.
- [x] Verify missing-versus-observed atom37 masks through the pinned actual
  cc94 reader without altering the engine or runtime.
- [x] Define separate atom-completion v1 wire requests/capabilities, with strict
  JSON decoding, source/target/mask binding and explicit resource ceilings.
  Existing schema-4 `DiffusionRequest` and frozen semantics remain unchanged.
- [x] Implement private generated-only materialization for partial residues and
  mixed whole gaps, preserving original ATOM/HETATM/ANISOU bytes and valid links.
- [ ] Integrate the separate protocol with runner execution; capabilities deny
  all execution by default. No existing profile advertises completion support.
- [ ] Teach context selection, engine views and PDB materialization to insert only
  missing atoms while retaining known atoms. Private insertion is implemented;
  engine/context integration remains open. The legacy materializer remains
  whole-gap-only; the new adapter never supplies it fixed-coordinate updates.
- [ ] Integrate private generated-only refinement/validation with source-bound
  independent publication and `.dat` ownership for mixed observed/generated residues.
- [ ] Run predeclared checkpoint-backed tests before advertising capability or
  enabling CLI/model/Homology/batch/resume/GUI integration.

- [x] Define backend-neutral **atom-level** completion requests: exact observed
  and missing atom identities, expected residue chemistry and separate masks
  for whole-residue gaps versus partial-residue completion. Preserve chain case,
  residue numbers and insertion codes. Do not mark an entire observed residue
  as generated merely because one side-chain atom is absent.

### Separate contract and private materialization slice

`atom_completion_contract.py` adds a schema-1 `AtomCompletionRequest` and
`AtomCompletionCapabilities`, reusing strict JSON encoding without changing
schema-4. Requests bind source SHA, canonical partial masks, exact placements,
complete gap atom sets, unique bounded seeds, template ownership and an explicit
profile. Capability defaults deny execution; feature checks and runner resource
ceilings are separate from scientific admission. Contract import works with
`python -S`: scientific package re-exports are now lazy while retaining the
existing public names. Canonical atom tables and their ordering were independently
compared against `931f977` and are unchanged.

`atom_completion_materialize.py` accepts **generated atoms only**. It derives names
and elements from the request, validates finite/fixed-column coordinates, allocates
unused serials and inserts missing atoms after their observed residue/ANISOU block.
Mixed gaps reuse a legitimate gap-only adapter, not fake gaps for observed residues.
Original coordinate and ANISOU records must remain byte-identical and in order.
Ambiguous MODEL records, split residue blocks, duplicate serials, orphan ANISOU and
dangling CONECT are refused; stale MASTER counts are dropped and records after END
are rejected. CRLF and valid source links are retained. The 31 new
contract/materialization tests use synthetic coordinates to check data integrity,
**not scientific geometry or successful diffusion reconstruction**.

The continuation passed **147 focused tests** and **1,871 fast tests** (3 explicit
Torch/Biotite skips, 33 slow deselected) in 162.72 s. All 91 existing package
re-exports resolve; separate Python 3.11 stdlib-only contract roundtrip, Ruff,
targeted mypy and documentation checks passed. No new scientific sampling claim
is inferred from these implementation tests.

### Private validation and refinement slice (2026-10-31)

`atom_completion_validation.py` independently audits candidate bytes against the
source and request, without trusting runner or materializer declarations. It
requires the exact source-plus-generated atom set; source observations (including
hydrogens), ANISOU and explicit links remain verbatim and ordered. Generated names,
elements and residue chemistry are checked. Shared diagnose routines cover bond
lengths/angles, amide planarity, chirality, Ramachandran, chi12, coincidences and
sterics; explicit target-neighbour C-N checks catch disconnected bonds even if a
parser omits them. The changed-residue/junction/neighbour cohort is gated, and
D-Cα is rejected globally. All source missing atoms/gaps and geometry findings
remain in the report, including outside the changed cohort. A private pass does
not certify unrelated source geometry or grant publication.

`atom_completion_refinement.py` delegates exactly the generated atom mask to the
existing v8 OpenMM kernel and reconstructs the output from unchanged source bytes.
It never grants observed mobility or marks entire partial residues as generated.
This first adapter is **protein-only** and refuses raw D-Cα instead of reflecting
observed side chains. Kernel failures retain raw bytes and diagnostics; both raw
and refined reports survive, and refined rejection is not silently treated as
acceptance or hidden behind a raw fallback.

Three predeclared, actual Reference controls passed: LYS A27 four-atom completion
from 1UBQ reference coordinates, the same control with a +0.75 Å NZ x perturbation,
and a mixed control with A10–A13 deleted (32 missing-mask atoms total). Seed 7,
250 iterations, zero restarts, unchanged observed records. Both initial and final
geometry passed for all three; the mild perturbation was not claimed as a raw
failure. **These are reference-derived controls, not diffusion output.** All
artifacts and full source diagnostics remain under
`tmp/atom-completion-e2904a9/audited-controls/` (including module digests and runtime
versions; earlier controls remain retained); see the
[portable control record](../research/partial-residue-refinement-controls-2026-10-31.json).

- [x] Private generated-only refinement and independent source/candidate integrity
  plus geometry validation, including partial-only and mixed reference controls.
- [ ] Runner/context/engine transport, sampler trace/provenance, canonical `.dat`
  ownership and independent atomic bundle publication for completion candidates.
- [ ] Checkpoint-backed atom completion and predeclared scientific qualification;
  keep public options disabled until those gates pass.

Verification: **165 focused / 1,889 fast tests passed**, 3 explicit runtime skips,
33 slow deselected; Ruff, targeted mypy and documentation checks passed. The
fast lane completed in 165.29 s. No checkpoint-sampling claim is inferred.

### Private cc94 checkpoint pilot (2026-10-31, after `eacd0d5`)

`atom_completion_sampling.py` adds a separate full-monomer engine view and narrow
source admission: one <=128-residue heavy-atom chain, up to four partial
sidechains, and optionally one internal 3–12-residue gap. Backbone defects,
heterogens, explicit covalent context, partial occupancy, terminal gaps and
unsupported resource/context requests are denied before weights load. The view
keeps exact source identities (including chain case and insertion codes), uses
native relative target positions and transports actual observed/absent atom37
masks. Only declared, spatially disconnected gap-bridge bonds inferred by PDBFile
are excluded from the temporary diagnostic topology; original findings and exact
exclusions remain in evidence. A bonded observed peptide cannot be reclassified
as a gap.

The standalone **private** adapter
`deploy/protpardelle-1c/atom_completion_research.py` requires explicit experimental
acknowledgement. It reuses pinned cc94 preflights, strict weights/config identity,
native 500-update ODE conditioning, actual-device checks and the relative-position
probe. It writes no production `RunnerResult` or bundle and never registers a
public profile. Original inputs, native atom37 arrays, reader masks, engine/source
maps, trace, raw/rejected proposals and refinement diagnostics are retained in
exclusive output reservations, never adopted or overwritten.

Two predeclared seed-7 MPS cases **actually sampled 500 updates each**: partial
LYS A27 (4 missing atoms) and that partial residue plus A10–A13 (32 missing atoms).
The original 76-residue protein projection was refused for occupancies 0.45/0.25
at residues 73–76. An explicitly labelled **1–72 protein subset** was used instead;
occupancies and original fixture bytes were not rewritten. Missing coordinates
were withheld from the sampler. An initial adapter import-order failure is
retained separately; it sampled zero updates and was not a scientific rejection.

**Both raw diffusion proposals failed geometry.** Separate generated-only CPU
refinement produced passing candidates, independently revalidated from source
bytes in another process. The missing atoms' raw coordinates were independently
recomputed from saved native output and the declared rigid fit/rounding. No
canonical atoms were missing before refinement; temporary OXT A72 completion was
not published. Known observations were byte-identical in raw and refined outputs.
Generated-atom refinement RMSD was **3.92 Å** (partial) and **1.58 Å** (mixed), so
this demonstrates a diffusion-proposal-plus-refinement path, **not good native
sampling quality or broader qualification**. Training membership is unknown.

Evidence: `tmp/atom-completion-eacd0d5/pilot-v2/` and its separate independent
audit; [portable pilot record](../research/partial-residue-cc94-pilot-2026-10-31.json).
The new 11 tests cover scope refusal before weights, native transport, exact
identity mapping, occupancy refusal and preservation of known peptide evidence.

- [x] Narrow private source admission, engine transport and checkpoint-backed
  partial/mixed MPS pilots, with source/native/refinement evidence retained.
- [ ] Production runner/result/trace integration, canonical `.dat` ownership and
  independently audited atomic publication; no public profile is enabled.
- [ ] Prospective multi-seed/multi-structure qualification and non-monomer context;
  the two pilot passes do not close these gates.
- [ ] Compare against template-only atom rebuilding before claiming a scientific
  benefit: the substantial refinement required here establishes neither superiority
  nor good native proposal quality.

Final implementation verification: **1,900 fast tests passed**, 3 explicit
Torch/Biotite skips, 33 slow deselected, 171.33 s; Ruff, targeted mypy,
documentation and Python 3.11 CLI-reference checks passed.

### Separate trace, `.dat` and experimental publication

`atom_completion_trace.py` defines a strict stdlib-only completion trace v1;
existing schema-4 `SamplerTrace` and public gap requests are untouched. Native
cc94 conditioning is shared by execution and auditing. The research adapter
automatically supplements successful attempts with a bound trace and the exact
adapter source snapshot after its progress file is finalized. Historical/private
attempts can be supplemented only with source bytes matching their recorded
executed digest, without rewriting the original requests, native arrays or
progress. Existing attestations are never adopted/overwritten.

`atom_completion_evidence.py` independently audits all digest-bound artifacts,
preflight/policy, source admission, engine mapping, motif and reader masks. It
reconstructs raw generated coordinates from native atom37 plus the declared fit
and PDB rounding, then independently rechecks raw/refined observations and final
geometry. NPZ member sets, expanded byte counts and NPY shape/dtype headers are
bounded before allocation; objects/pickle and nonfinite tensors are rejected.
This is chain-of-custody/integrity evidence, not a cryptographic proof of
execution on an adversarial host.

`atom_completion_publication.py` builds a canonical `DatRecord` containing exactly
the requested missing atoms, never existing observations or temporary OXT/H.
Chain case and insertion codes are preserved in ownership. It publishes only
independently passing candidates to a distinct
`experimental-atom-completion-bundle`, with original source/native/raw evidence
and both independent reports. Payloads are durable before an exclusive manifest-last
commit via the existing CephFS-compatible helper. Missing index, corrupted files,
undeclared entries and stale/crash-left destinations fail closed. Failed staging
is retained; original run artifacts are not modified or deleted. The ordinary
public gap-bundle reader explicitly rejects this different experimental kind.

Both retained pilot candidates were independently re-audited and published in
`tmp/atom-completion-661a8b2/published/`, without resampling: **4 and 32 `.dat`
atoms**, zero known observations marked as new. Their original raw failures and
large refinement movement remain visible; publication is not scientific
qualification. See the [portable publication record](../research/partial-residue-publication-2026-10-31.json).

- [x] Separate typed completion trace, native/raw linkage audit, exact canonical
  `.dat` ownership and manifest-last experimental publication/readback.
- [ ] Production schema/API/CLI/model/Homology/GUI integration after qualification;
  no public profile is enabled by experimental bundle readiness.

A fresh **prospectively declared** 2-case × 5-seed cohort (7/19/42/101/202)
executed all **10 × 500 MPS updates**, with unchanged source/conditioning/refinement
policy and no fallback. Automatic trace → audit → `.dat` → manifest-last publication
and a separate final readback passed for all ten bundles. **Raw geometry passed
0/10; refined geometry passed 10/10.** Generated-atom refinement RMSD ranged
1.21–3.92 Å. These are operational/case-specific results on the same cropped
protein, not broad qualification or superiority over template-only reconstruction.
The cohort and all raw/native/refined diagnostics are retained under
`tmp/atom-completion-661a8b2/prospective/`. A durability-only flush/fsync addition
to attestation writing occurred during the cohort; the original declared module
hash snapshot remains retained, scientific policies were unchanged, and every
bundle was independently re-audited with the final reader. Do not treat this as
a frozen-source hardware-repeat acceptance.

Verification: 22 new publication tests; **82 focused / 1,922 fast tests passed**,
3 explicit Torch/Biotite skips and 33 slow deselected, 192.47 s. Ruff, targeted
mypy, user/agent docs and Python 3.11 CLI-reference checks passed.

### Template comparison and additional structures (after `460cce9`)

`atom_completion_benchmark.py` adds a **separate template baseline**, not a fake
checkpoint trace. It calls the existing seeded PDBFixer rebuild helper, maps any
whole-gap insertion from the already-bound target placement and preserves known
coordinates. It never receives held-out coordinates. The same generated-only
CPU refinement, seeds and geometry gates are used in both arms. One deterministic
template heavy-atom proposal is reused across five refinement seeds; those are
not five independent template proposals. Exact-name held-out RMSD is measured
without fitting in the unchanged source frame, with all failures retained.

The old 1UBQ diffusion cohort was evaluated **retrospectively**. Before new
sampling, two additional cases and seeds 7/19/42/101/202 were fixed: 1FR2's 83
observed chain-A residues, deleting ASN24 CG/OD1/ND2; and an explicit protein-only
1A3N A1–59 subset, deleting GLN54 CG/CD/OE1/NE2. Partners/heme are excluded in these
labelled projections; original structures and occupancies are not rewritten.
All ten new cases actually sampled 500 MPS updates. Reference coordinates went
only to the independent scoring process, not either proposal/refinement arm.

| Partial sidechain | cc94 final geometry | Template final geometry | cc94 median held-out RMSD | Template median |
|---|---:|---:|---:|---:|
| 1UBQ subset / LYS27 | 5/5 | 5/5 | 1.458 Å | 1.235 Å |
| 1FR2 chain A / ASN24 | 5/5 | 5/5 | 1.766 Å | 0.992 Å |
| 1A3N A1–59 / GLN54 | 4/5 | 5/5 | 2.988 Å | 2.243 Å |

Medians include all five exported refined candidates, not a best-seed subset.
The rejected GLN54 seed 19 has chi12/clash failures and 5.944 Å RMSD; the accepted
four's median is still 2.984 Å. **No comparative advantage is demonstrated**:
the template baseline is closer to the withheld deposited coordinates at each
tested site. This supports retaining the partial-completion route as private
research, not enabling a public default. It is not a general checkpoint-family,
training-domain or physical-ground-truth conclusion.

The mixed 1UBQ restricted PDBFixer baseline failed one shared proposal after five
internal seeds with D-Cα at residues 11–13, so none of its five planned outer
refinements executed. That is **not a comparison against full legacy Modeller
loop reconstruction or its downstream chirality safeguards**. Keep the refusal
and the missing stronger comparator explicit rather than calling all template
methods inferior for gaps.

Evidence: `tmp/atom-completion-460cce9/comparison/` and the
[portable comparison record](../research/partial-residue-template-comparison-2026-10-31.json).
The first baseline controls also uncovered missing conda-owned PDBFixer assets
in three inspected general prefixes (15 files each, including `soft.xml`), while
the isolated CPU runtime has its complete PDBFixer package. No prefix was patched
or deleted. Doctor now reports this resource failure and selected-profile readiness
fails closed instead of equating package importability with working rebuilding.

- [x] Limited sidechain template comparison and two additional projected-structure
  cohorts, with negative results and all rejections retained.
- [ ] Demonstrate benefit on a prospectively defined broader cohort and compare
  whole gaps against a representative full legacy/Modeller baseline before public
  partial-completion qualification. Current sidechain evidence does not justify
  public enablement.
- [ ] Audit each backend/checkpoint's ability to condition on partial residues
  and generate the requested heavy atoms, including side-chain-only defects,
  backbone defects and mixed defects adjacent to sequence gaps. Add maintained
  adapters only for demonstrated capabilities; retain explicit refusals for
  unsupported chemistry or missing conditioning atoms. Do not disguise
  template-built atoms as diffusion-generated output.
- [ ] Keep source assessment before motif selection, but distinguish absent
  atoms from invalid existing geometry. Never invent missing source atoms just
  to satisfy context construction or a force-field preflight. Define how motif
  selection and optional refinement operate with incomplete conditioning.
- [ ] Keep known atoms fixed by default. Any necessary movement must use an
  explicit, bounded, independently audited relaxation policy, with original
  coordinates, atom provenance, masks and all attempts preserved. Frozen
  requests must not be silently repaired or reclassified.
- [ ] Validate output completeness, exact names/identities, bond connectivity,
  junction geometry, clashes, chirality and observed-atom preservation using
  shared diagnose/domain primitives. Publish generated-atom ownership and
  `.dat` evidence through the canonical schema; accept no unrequested atoms or
  duplicate identities. Handle hydrogens separately from heavy-atom completion.
- [ ] Add CLI/model and mosaic Homology tests, then batch/resume and GUI coverage,
  for partial residues alone and combined with multiple whole-residue gaps.
  Include case-sensitive chains, insertion codes, disulfides, unsupported inputs
  and failure-evidence retention; keep heterogen support a separate capability.
- [ ] Freeze checkpoint-backed multi-seed CUDA/MPS benchmarks before execution:
  artificial heavy-atom deletions from complete held-out sources plus authentic
  incomplete inputs such as 5J13. Preserve originals and label any protein-only
  derivative separately. Report admission refusals, executed failures and
  accepted candidates separately; never relax gates or select seeds by reference
  fit. Record backend-specific qualification before exposing broader support.

See the [current Linux continuation handoff](installation-profiles-apple-acceptance.md#latest-native-checkpoint-and-linux-continuation)
for work that can start now, without claiming this future capability.

## Current Expansion Baseline

**Public cc89 CUDA slice:** `protpardelle-1c-cuda` is now exposed through installed
runtime discovery, CLI and GUI. `auto` context selects the corrected fixed local
policy; the admitted slice is one internal 5–10-residue gap, canonical proteins,
at most 512 whole-source residues. Sampling and CUDA refinement use separate
environments/processes. Four public model/Homology workflows accept 5/5 candidates,
including multi-seed publication and verified resume; independent audit confirms
zero D-Cα and exact outside-gap records. This does not expose full-complex neural
sampling or joint gap-bearing chains, or change the historical MPS profile/MODELLER default.
See the [public acceptance](../research/cc89-public-cuda-acceptance.md) and
[deployment instructions](../../deploy/protpardelle-1c/README.md#public-linuxnvidia-cuda-profile).

**Input-representation correction:** actual attention/RoPE audits retain chain
position breaks, but the full Fab/ASU denoising axis and proxy-fragment position
features were inappropriate for the local repair experiment. Explicit local
conditioning and physical source positions reduce raw CDR1 RMSD from 4–32 Å to
0.83–1.35 Å; the ordinary adapter's unchanged refinement then passes all 16
original CDR mask/seed pairs, median 1.152 Å. This same-source diagnostic repair
does not replace the original holdout or establish independent acceptance.
The complete general-protein recheck now passes 44/44 after unchanged refinement,
retaining all 37 original passes and recovering seven failures. Combined corrected
acceptance is 60/60, with zero detectable D-Cα and exact outside-gap atom-record
preservation. Geometry/source preservation is the primary outcome; reference RMSD
is a secondary conformation descriptor, not an acceptance gate.
See [the input audit](../research/cc89-input-representation-audit.md).
The [source-disjoint prospective follow-up](../research/cc89-prospective-holdout.md)
passes 4/4 on its first new general-protein source, 9WW8. A separate frozen
8EK6/9Z0U/9LDN cohort now passes 20/20 after unchanged refinement. Combined
prospective acceptance is 24/24 on four new sources, zero detectable D-Cα and
exact outside-gap preservation; this remains a bounded source-disjoint check.

**Expanded known-sequence holdout:** all 154 planned jobs on 13 X-ray source
structures are complete and independently revalidated. Under the explicit
atom37 terminal-preservation policy (observed OXT retained as preserved-only),
cc89 passes 16/60 raw and 37/60 after unchanged CUDA refinement: 37/44 general
protein outcomes, 0/16 antibody CDR outcomes, with three guarded chirality
failures retained. MODELLER's unconstrained comparator geometry passes 30/30;
four public Protenix/Homology controls pass. This is evidence against promoting
the current full-context cc89 mode broadly, not completion of that acceptance
gate. See [the known-sequence study](../research/cc89-known-sequence-holdout.md).

**Subsequent cc89 chain-phase investigation:** 42 frozen CUDA samplings isolate
the denoiser's extra rotary chain-phase component as the cause of the previous
11–20 Å reconstruction failures. An explicit CUDA-only research repair preserves
all input atoms/gaps and sampler chain/residue features and reproduces causal
controls byte-for-byte. The unchanged generated-only refinement passes 6/6 at
0.17–0.25 Å backbone RMSD. This repaired pilot does not complete wider source/seed
acceptance; see [the causal report](../research/cc89-chain-phase-root-cause.md).

**2026-10-19 CUDA completion audit:** shared profile-specific length admission,
adaptive-source hardening, and Homology deletion-window/failed-request regression
coverage are implemented. All 50 full-input clean-cohort outcomes were audited for
withheld accuracy: 49/50 refined geometry passes include a length-50 result at
8.69575 Å backbone RMSD. Rich cc89 CUDA refinement accepts 0/6; two cases stop at
the chirality guard with no refined output. These are open scientific gates,
not missing adapter wiring. See [completion evidence](../research/cuda-expansion-completion.md).
The subsequent corrected Apple port/public integration is now complete in its
bounded local-v1 scope; see [Apple acceptance](../research/cc89-public-apple-acceptance.md).
Its larger operational cohort remains pending and does not replace NVIDIA evidence.

**Input-scope correction after manual review:** the
[chain-truncation audit](../research/a100-input-truncation-audit.md) reproduces
substantial removal of observed residues/chains before artificial gap creation.
The A100 studies below concern cropped fragments; their fixed-environment checks
do not establish preservation of the complete deposited structural context.

The corrective [intact-input recheck](../research/intact-input-recheck.md) is
complete: 47 original masks retained, 25 admitted cases attempted, 13 accepted
after rigid local refinement, 11 returned geometry failures and one guarded
chirality rejection with no refined output. The 22 incompatible masks remain
explicitly unsupported with their full CIF inputs retained. This recheck uses
local engine conditioning while preserving the complete observed source in the
scientific files and force-field context; it does not widen public scope.

The next [prospective A100 calibration and cc89 control study](../research/cuda-calibration.md)
has frozen source membership and calibration/validation splits. The 15 cc89
controls are complete and independently revalidated (0/15 accepted). All 92
prospective fragment outcomes have completed and passed the offline evidence
audit. A separate peptide-junction audit found no measurement or pairing defect;
it distinguishes bond-length acceptance from coarse chain-break detection.
The subsequent [paired local-relaxation study](../research/a100-failure-and-relaxation.md)
is complete: 33/92 raw versus 86/92 relaxed accepted, 53 rescues and zero
pass→fail pairs, with six long-gap failures retained. Reference accuracy can
worsen despite geometry acceptance; scientific/public scope remains unchanged.

The completed [relaxation follow-up](../research/relaxation-followup.md) audits
reference accuracy on 64 nonduplicated mask/seed cases. Bounded-anchor refinement
does not outperform rigid (57/64 versus 58/64), so the frozen rule retains rigid.
On three disjoint September–October sources, geometry acceptance improves from
6/15 to 12/15; three accepted long-gap candidates still have reference RMSD near
6 Å. Public scope remains unchanged. New-source outputs did not tune the method.

Further implemented research slices and 31 verified frontier outcomes are recorded
in [reconstruction frontiers](../research/reconstruction-frontiers.md): private
joint/discontiguous cc89 views, two-round sampled-trajectory context expansion,
explicit Homology edit flanks, and 13–25/26–50 research capability strata. Rich cc89
geometry and adaptive improvement are not generally accepted. Public admission
is profile-specific: 3–12 for Protenix/legacy MPS and 5–10 CYS-free internal gaps
for corrected cc89 CUDA/MPS; research strata do not inherit public acceptance.

Blocks WP5–WP8 and WP12 now have an additional implemented slice and frozen
expanded studies: [context/Homology expansion protocol](../research/protein-context-expansion.md).
The shared selector handles joint Protenix gap-bearing chains; Protpardelle uses
the shared compatibility/explicit selectors with stated adapter boundaries.
Homology uses authoritative coverage requests and opt-in bounded substitutions
and deletion repair. Expanded outcomes, sensitivity, causal audit, and paired
Homology comparisons require completed verified archives before acceptance claims.

Workflow implementation update: the current scoped slice of WP9–WP11 and WP12
orchestration is implemented and verified. This includes backend-dependent output
metadata, sequential diffusion batch/resume, validated ZBS PDB/`.dat` handoff,
Protenix/Protpardelle mosaic selection, and GUI native publication/workspace
archival. See [workflow integration evidence](../research/diffusion-workflow-integration.md).
The remaining unchecked scientific expansion tasks below (covered substitutions,
broader context contracts, interface cohorts, and new gap-length strata) are not
implied by that workflow slice. WP3 remains paused/cancelled pending a separate
decision to resume validation-policy work.

- [x] Multichain input, multiple ordinary gaps, and independent per-gap-bearing-
  chain Protpardelle invocations are implemented.
- [x] Protpardelle receives deterministic contiguous crops from neighboring
  canonical protein chains when fixed heavy atoms are within `12 A` of a gap
  anchor. Partner atoms enter the real multichain motif axis and are recorded in
  sampler trace context mappings.
- [x] The current partner selection expands each contact crop by two sequence
  neighbors, retains contiguous observed runs, and obeys the 512-residue axis
  budget.
- [x] N-terminal and C-terminal one-anchor generation is implemented for the
  accepted Protpardelle profile, with broader frozen hardware cohorts pending.
- [x] An initial one-chain mosaic-first Homology insertion slice is implemented
  in CLI and GUI with authoritative coverage metadata and fixed covered atoms.
- [x] An isolated authoritative-SMILES fixed-ligand geometric-repulsion slice is
  implemented but remains separate from learned surrounding-protein context and
  from general heterogen support.
- [ ] Existing local-partner evidence is a checkpoint smoke, not the required
  frozen interface cohort or causal partner-on/off ablation.
- [ ] Multiple gap-bearing chains are not jointly sampled and the current crop
  does not adapt to the predicted gap trajectory in the accepted Protpardelle MPS
  profile. Private CUDA joint/discontiguous views and between-round adaptive
  selection are implemented and tested; their broader scientific acceptance is
  still pending.

## Non-Goals

- [ ] Do not remove the current gap-length limit before each replacement stratum
  passes preregistered acceptance gates.
- [ ] Do not call restored or re-spliced neighboring coordinates
  `surrounding-aware`; the engine must consume verified context features.
- [ ] Do not concatenate separate chains into a fake continuous protein to work
  around a model that lacks multichain conditioning.
- [ ] Do not silently downgrade context-conditioned requests to gap-only,
  preservation-only, another profile, or MODELLER.
- [ ] Do not infer ligand, glycan, PTM, cofactor, metal, or covalent-link
  chemistry in this protein-context expansion. The heterogen plan owns those
  classes.
- [ ] Do not add ML frameworks, engine packages, or checkpoints to the core
  wheel or `environment.yml`.
- [ ] Do not expose a new public scope before its fixtures, thresholds,
  comparison baseline, and failure policy are frozen.

## Why Gap Length Remains A Capability Stratum

- The current 3-12 range is an evidence boundary, not an architectural claim
  about the maximum sequence length accepted by Protenix.
- Longer gaps increase conformational uncertainty, candidate-count
  requirements, closure failures, context sensitivity, and ranking error.
- A long missing segment can transition from local repair into de novo or
  full-chain prediction; those tasks require different evidence and output
  claims.
- Replace one global constant with explicit profile capabilities only after the
  following strata are evaluated independently:
  - [ ] 3-12 residues: retain as the accepted baseline.
  - [ ] 13-25 residues: first expansion stratum.
  - [ ] 26-50 residues: separate high-cost research stratum.
  - [ ] More than 50 residues: require a separate architecture and product
    decision before implementation.
- Record candidate count, runtime, RAM, VRAM, closure rate, hard-gate pass rate,
  diversity, and ranking enrichment by gap-length stratum.

## Execution Rules For Subagents

- [ ] Assign one coordinating agent as the sole owner of shared contract,
  validation-policy, CLI-dispatch, and DDD-map changes.
- [ ] Give each implementation agent a disjoint work package and explicit file
  ownership; agents must not opportunistically edit another package's files.
- [ ] Require every agent handoff to report:
  - changed files;
  - contract or schema changes;
  - commands run and exact outcomes;
  - unresolved assumptions;
  - hardware, checkpoint, and environment identity when applicable;
  - whether generated documentation or fixture manifests need regeneration.
- [ ] Land shared contract changes before adapters and consumers branch from
  them; do not reconcile parallel incompatible schemas after implementation.
- [ ] Keep outcome-bearing benchmark execution separate from threshold design.
- [ ] Commit logical change groups independently so a failed scientific scope
  can be reverted without removing validated infrastructure.
- [ ] Do not commit checkpoints, credentials, private URLs, raw host
  environments, or ignored benchmark workspaces.

## Dependency Graph

- [ ] `WP0 -> WP1` establishes the reproducible Linux baseline.
- [ ] `WP2 -> WP3` establishes reusable measurements before hard-gate policy.
- [ ] `WP2 + WP3 -> WP4` enables broader scope evaluation.
- [ ] `WP5 -> WP6 -> WP7` establishes context contracts, selectors, and engine
  adapters before context-conditioned claims.
- [ ] `WP3 + WP4 + WP7 -> WP8` enables surrounding-aware acceptance studies.
- [ ] `WP1 + WP3 -> WP9` enables backend-neutral orchestration hardening.
- [ ] `WP9 -> WP10 -> WP11` enables batch and ZBS integration.
- [ ] `WP5 + WP7 + WP9 -> WP12` enables mosaic-first Homology.
- [ ] `WP13` performs final repository and hardware verification after all
  approved packages complete.

## WP0: Baseline And Host Inventory

- Suggested owner: integration coordinator.
- Owned surfaces:
  - branch state and release metadata;
  - generated-reference status;
  - hardware/environment evidence under `docs/research/`.
- Tasks:
  - [ ] Start from the current remote `feature/dvbfixer-hardening` tip without
    including unrelated local or untracked files.
  - [ ] Record GPU model, compute capability, VRAM, driver, CUDA runtime, CPU,
    RAM, local storage, and available OCI runtime.
  - [ ] Record the operator-provided checkpoint path only where operationally
    required and store only its SHA-256 in tracked evidence.
  - [ ] Run the current CPU diffusion suite and preserve its exact result.
  - [ ] Run CLI-reference, GUI-spec, DDD-map, version, lint, and type checks;
    resolve pre-existing generated-reference drift as a separate change group.
  - [ ] Freeze the current profile metadata, protocol version, validation
    thresholds, and 231-case confirmatory aggregate before modifications.
- Exit gate:
  - [ ] The baseline is reproducible and any pre-existing failures are separated
    from failures introduced by later packages.

## WP1: Linux/NVIDIA Protenix Acceptance

- Suggested owner: platform runner agent with access to the NVIDIA host.
- Owned surfaces:
  - `deploy/protenix-v1/`;
  - Linux runner environment and image definitions;
  - profile-specific acceptance scripts and evidence.
- Tasks:
  - [ ] Pin Python, PyTorch, CUDA runtime, Protenix source revision, callback
    patch digest, OpenMM stack, auxiliary artifacts, and package lock.
  - [ ] Pin the base image by digest and record the final image digest when an
    OCI runtime is available.
  - [ ] Verify checkpoint SHA-256 before model loading; do not download or
    redistribute the checkpoint from core DVBFixer.
  - [ ] Assert CUDA availability, expected device class, model/input/output
    placement, dtype, callback count, and disabled fallback.
  - [ ] Bound timeout, process groups, stdout, stderr, artifact sizes, workspace
    use, RAM, and VRAM; clean all private workspaces after failure.
  - [ ] Run protocol preflight, one checkpoint-backed reconstruction, same-seed
    independent repeats, different-seed candidates, and failure injection.
  - [ ] Repeat the frozen 231-case cohort and compare hard-gate pass rate,
    selected geometry, fixed-coordinate adherence, runtime, RAM, VRAM, crash,
    timeout, and OOM rates with the accepted confirmatory baseline.
  - [ ] Add a manual or self-hosted NVIDIA acceptance lane while keeping core CI
    CUDA-independent.
- Stop gates:
  - [ ] Stop on silent CPU fallback, wrong-device tensors, source/patch/checkpoint
    mismatch, incomplete callback evidence, fixed-coordinate drift, or partial
    publication.
  - [ ] Do not widen public scope until the production wrapper reproduces the
    accepted baseline.

## WP2: Pure Geometry Measurement Layer

- Suggested owner: geometry agent.
- Owned surfaces:
  - `src/dvbfixer/diagnose/geometry.py`;
  - new pure measurement modules approved by the coordinator;
  - focused geometry tests.
- Tasks:
  - [x] Keep measurements free of severity and diffusion-specific pass/fail
    policy.
  - [x] Preserve full case-sensitive residue and atom identities, including
    insertion codes.
  - [x] Add class-aware Ramachandran measurements for general, GLY, PRO, and
    pre-PRO residues.
  - [x] Add residue-aware chi1-only, chi1/chi2, chi3-chi5, and
    backbone-dependent rotamer measurements where reference data support them.
  - [x] Add residue-specific bond-length, bond-angle, peptide-planarity,
    cis/trans, junction, and terminal measurements.
  - [x] Add expected intra-residue, peptide, and disulfide connectivity evidence
    without inferring unsupported chemistry.
  - [x] Provide deterministic heavy-atom and hydrogen-aware steric measurements
    for generated/generated, generated/fixed, and junction neighborhoods.
  - [x] Return explicit undefined results for missing, duplicate, degenerate, or
    unsupported geometry rather than a passing value.
- Exit gate:
  - [x] Existing unscoped `diagnose` behavior remains compatible unless a
    separately approved diagnostic policy change is documented and tested.

## WP3: Diffusion Validation Policy

- Suggested owner: exclusive validation-policy coordinator.
- Owned surfaces:
  - `src/dvbfixer/model/diffusion/validate.py`;
  - `src/dvbfixer/model/diffusion/quality.py`;
  - versioned validation manifests and focused tests.
- Tasks:
  - [ ] Map approved WP2 measurements to versioned diffusion hard gates.
  - [ ] Add GLY/PRO/pre-PRO Ramachandran policy.
  - [ ] Add chi1-only, chi3-chi5, and full residue/backbone-dependent rotamer
    policy where calibrated reference data exist.
  - [ ] Validate canonical atom completeness, duplicate atoms, identity,
    connectivity, peptide links, disulfides, planarity, clashes, and final
    chirality after the last heavy-coordinate change.
  - [ ] Preserve fixed-heavy RMSD at no more than `0.01 A` and maximum
    displacement at no more than `0.03 A` unless a separately versioned study
    changes those limits before outcomes are viewed.
  - [ ] Record every skipped or undefined metric with a stable reason; never
    interpret missing evidence as passing.
- Exit gate:
  - [ ] Every published candidate passes the same independent validator
    regardless of engine-native confidence or score.

## WP4: Reference-Dependent Benchmark Metrics

Implemented local/interface measurement slice and frozen-cohort descriptive
follow-up: [reference metrics evidence](../research/diffusion-reference-metrics.md).
GDT-HA/TM-score remain conditional on meaningful declared large-region inputs and
a reviewed global-alignment engine; they are not reported for short-loop cases.

- Suggested owner: benchmark agent.
- Owned surfaces:
  - `src/dvbfixer/model/diffusion/benchmark.py`;
  - benchmark schemas, aggregation scripts, and tests.
- Tasks:
  - [x] Add generated-region and junction-neighborhood lDDT.
  - [x] Add interface lDDT and contact precision/recall for declared interface
    cases.
  - [ ] Add GDT-HA only where the evaluated region is large enough for the
    metric to be interpretable.
  - [ ] Add TM-score only for declared full-domain or sufficiently long-region
    comparisons; do not report it for short loops.
  - [x] Preserve gap backbone/all-heavy RMSD, fixed/anchor RMSD, top-1, oracle
    top-k, pairwise diversity, ranking enrichment, and failure classes.
  - [x] Report paired MODELLER comparisons and confidence intervals by stratum.
  - [x] Keep engine confidence, MODELLER `molpdf`, physical energies,
    independent geometry, and reference metrics as separate named scales.
- Exit gate:
  - [x] No aggregate score hides a hard-gate failure or compares unlike
    objectives as if they shared a physical scale.

## WP5: Backend-Neutral Context Contract Expansion

Implemented Protenix context-v1 slice and bounded causal pilot:
[implementation and evidence](../research/protein-context-pilot.md). Schema-4
optional fields preserve legacy serialization. The tasks below remain the full
cross-backend acceptance target, not a claim that the pilot completes WP5–WP8.

- Suggested owner: contract coordinator; this package lands before any adapter.
- Owned surfaces:
  - `src/dvbfixer/model/diffusion/contract.py`;
  - protocol serialization and compatibility tests;
  - provenance and sampler-trace schemas.
- Existing baseline:
  - [x] Schema version 4 records exact partner-chain/atom conditioning contexts
    and verifies represented fixed atoms in sampler traces.
  - [x] Existing provenance distinguishes local protein-partner context and
    fixed heterogen context.
- Tasks:
  - [ ] Extend the protocol to represent four disjoint atom roles consistently
    across protein, mosaic, and future supported chemistry context:
    - `generated`;
    - `movable-junction`;
    - `fixed-conditioning`;
    - `preserved-only`.
  - [ ] Represent context residues and atoms with exact chain, residue number,
    insertion code, atom name, coordinates, and source identity.
  - [ ] Record context-selection policy, shell thresholds, complete selected
    identity digest, represented-context digest, and engine mapping digest.
  - [ ] Require the runner result and sampler trace to state which context atoms
    reached the model, which were projected after each update, and which were
    only restored after sampling.
  - [ ] Add stable capability labels:
    - `gap-only`;
    - `same-chain-conditioned`;
    - `multichain-conditioned`;
    - `preservation-only`.
  - [ ] Reject old protocol versions when they cannot prove context semantics;
    do not infer capability from absent fields.
- Exit gate:
  - [ ] A result cannot claim surrounding-aware behavior unless the trace proves
    that the declared fixed context was represented and consumed by the sampler.

## WP6: Deterministic Context Selection Expansion

The shared source-only selector now supports sequence flanks, 6/10 Å
observed-anchor shells, complete selected residues, discontiguous spans, and
fail-closed budgets for Protenix context-v1. The existing Protpardelle selector
and trajectory-adaptive selection remain separate follow-ups.

- Suggested owner: context-selection agent.
- Owned surfaces:
  - a new backend-neutral context-selection module under
    `src/dvbfixer/model/diffusion/`;
  - scope/admission integration and focused tests.
- Existing baseline:
  - [x] Protpardelle selects contiguous partner-chain crops within `12 A` of the
    original anchors, expands by two sequence neighbors, and enforces a shared
    512-residue axis budget.
- Tasks:
  - [ ] Always include both anchors and complete selected residues.
  - [ ] Compare the current selector with configurable same-chain sequence
    flanks, initially evaluated at 8 and 16 residues on each side.
  - [ ] Compare the current `12 A` anchor selector with complete-residue spatial
    shells preregistered before outcome-bearing runs:
    - direct context: at most `6 A` from a generated or junction atom;
    - secondary context: greater than `6 A` and at most `10 A`;
    - distant context: ablation-only negative control.
  - [x] Expand selected Protpardelle contacts into contiguous chain fragments
    and record the expansion.
  - [ ] Generalize the expansion policy into a backend-neutral selector rather
    than duplicating it in each engine adapter.
  - [ ] Enforce a deterministic context budget with priority order: anchors,
    direct contacts, interface residues, secondary shell, then sequence flanks.
  - [ ] Preserve case-distinct chains and insertion codes; never select through
    concatenated display strings or numeric residue ranges alone.
  - [ ] Return a stable unsupported reason when required context exceeds engine
    capability or budget; never silently drop neighboring chains.
- Exit gate:
  - [ ] Identical normalized input and policy produce byte-identical context
    identities, ordering, and digests.

## WP7: Engine Context Adapters

### WP7A: Protenix

The context-v1 adapter maps selected entities and exact atom identities,
projects represented fixed coordinates at every update, and limits alignment
to the gap-bearing primary span. Six frozen local/full pilot probes demonstrated
pre-projection denoiser dependence on partner coordinates; see the linked
WP5 evidence. This is a bounded implementation and mechanistic result.

- Suggested owner: Protenix adapter agent.
- Owned surfaces:
  - `deploy/protenix-v1/`;
  - profile adapter and trace tests.
- Tasks:
  - [ ] Audit the pinned engine's multichain and template representation rather
    than relying on upstream marketing or schema acceptance.
  - [ ] Map selected neighboring protein fragments onto a stable engine atom
    axis and prove a round trip to exact DVBFixer identities.
  - [ ] Pass neighboring chains as fixed conditioning entities.
  - [ ] Reinject every represented fixed coordinate after every denoising update.
  - [ ] Fail on feature loss, atom reordering without a verified map, chain
    collapse, callback mismatch, or context drift.
  - [ ] Keep refinement movement limited to generated and explicitly movable
    junction atoms.

### WP7B: Protpardelle Capability Audit

- Suggested owner: compact-model agent, isolated from the Protenix adapter.
- Owned surfaces:
  - `deploy/protpardelle-1c/` research/audit scripts;
  - compact-model evidence documents.
- Existing baseline:
  - [x] The pinned model accepts neighboring canonical protein crops as separate
    chains on the real multichain motif axis.
  - [x] Exact target-to-partner chain and atom mappings are recorded and checked
    against represented fixed sampler atoms.
- Tasks:
  - [ ] Prove causal denoiser use with controlled partner-on/off and perturbed-
    context ablations; motif-axis representation alone is insufficient for the
    final scientific claim.
  - [ ] Quantify which discontinuous same-chain, multiple-chain, chain-identity,
    and fixed-coordinate features affect generated coordinates.
  - [ ] Retain `local-partner-conditioned` rather than whole-complex terminology
    until joint-sampling and frozen interface gates pass.
  - [ ] Do not concatenate neighboring chains with the target sequence.

### WP7C: Purpose-Trained Compact Context Model, If Required

- Suggested owner: model-research agent; begins only after WP7B records a
  blocker and a separate implementation approval is granted.
- Tasks:
  - [ ] Design an invariant or equivariant fixed-context encoder.
  - [ ] Preserve chain-aware and insertion-code-aware residue embeddings.
  - [ ] Add cross-attention or an equivalent explicit interaction between the
    generated gap and fixed environment.
  - [ ] Initialize or diffuse generated atoms only and project fixed represented
    coordinates after every update.
  - [ ] Train with same-chain, interface, context-dropout, and negative-control
    masks using leakage-controlled structure splits.
  - [ ] Freeze dataset, split, source, configuration, and checkpoint provenance
    before confirmatory evaluation.
  - [ ] Keep weights and runtime dependencies outside core DVBFixer.

## WP8: Surrounding-Aware Ablation And Acceptance

Eight arms have completed on three eligible frozen pilot cases, with one
unsupported source retained in the denominator. Independent current-policy
revalidation reproduces all 24 outcomes. Broader strata, negative controls,
sensitivity studies, and independent scientific review below remain pending.

- Suggested owner: scientific evaluation agent who did not implement the engine
  adapter being evaluated.
- Owned surfaces:
  - frozen manifests and result aggregation under `docs/research/`;
  - no production policy code.
- Required ablation arms, using identical structures, masks, and seeds:
  - [ ] anchors and gap only;
  - [ ] same-chain sequence flanks;
  - [ ] same-chain spatial environment;
  - [ ] same-chain plus neighboring-chain fragments;
  - [ ] full supported protein environment;
  - [ ] context plus per-step fixed-coordinate projection;
  - [ ] context, projection, and localized boundary refinement.
- Required strata:
  - [ ] solvent-exposed loops without interchain contacts;
  - [ ] buried and protein-core loops;
  - [ ] chain-interface gaps;
  - [ ] antibody CDRs near the H/L interface;
  - [ ] oligomeric interfaces;
  - [ ] neighboring-chain contacts within 3-6 A;
  - [ ] distant neighboring chains as negative controls;
  - [ ] insertion-code and case-distinct-chain cases.
- Required context metrics:
  - [ ] interface lDDT;
  - [ ] contact precision and recall;
  - [ ] minimum interchain distance and severe interchain clashes;
  - [ ] gap orientation relative to neighboring chains;
  - [ ] gap RMSD after alignment on fixed environment only;
  - [ ] change from gap-only to context-conditioned inference;
  - [ ] sensitivity to shell size and context budget;
  - [ ] fixed-context RMSD and maximum displacement;
  - [ ] ranking enrichment with and without context.
- Hard gates:
  - [ ] The trace proves that claimed context reached the sampler.
  - [ ] All fixed represented context passes existing coordinate thresholds.
  - [ ] No chain or insertion-code identity is lost.
  - [ ] No false interchain bond is introduced.
  - [ ] Context does not create a severe interface clash.
  - [ ] Unsupported context fails before inference without automatic fallback.

## WP9: Backend-Neutral Model Orchestration

- Suggested owner: integration coordinator.
- Owned surfaces:
  - shared model orchestration and backend dispatch;
  - command metadata and output contracts.
- Tasks:
  - [ ] Introduce backend-neutral model request, options, candidate, outcome, and
    artifact-set types.
  - [ ] Keep CIF normalization, sequence placement, chain discovery, CONECT/link
    handling, and output naming in common preprocessing.
  - [ ] Keep MODELLER behavior and file output unchanged.
  - [ ] Keep diffusion directory bundles and validation/provenance semantics
    explicit rather than hiding them behind file-only assumptions.
  - [ ] Expose backend capability metadata for gap length, terminal gaps,
    multiple gaps, multichain context, batch, ZBS, and Homology.
  - [ ] Reject unsupported option combinations before preprocessing or runner
    launch.

## WP10: Batch Diffusion

- Suggested owner: batch agent.
- Owned surfaces:
  - `src/dvbfixer/batch.py`;
  - command registry metadata;
  - batch-specific tests.
- Tasks:
  - [ ] Make output mode and suffix backend-dependent.
  - [ ] Publish each diffusion result as `<stem>_model_diffusion/`.
  - [ ] Preflight every destination and collision before the first runner launch.
  - [ ] Isolate workspaces and cleanup per input.
  - [ ] Add bounded GPU concurrency, defaulting to one inference job per GPU.
  - [ ] Define explicit seed lists or deterministic per-input seed derivation and
    record it in provenance.
  - [ ] Produce a bounded aggregate report with success, unsupported, failed,
    and skipped outcomes.
  - [ ] Preserve no-partial-bundle and source/output containment guarantees.

## WP11: ZBS Diffusion Integration

- Suggested owner: ZBS workflow agent.
- Owned surfaces:
  - `src/dvbfixer/zbs.py` and ZBS-specific tests/docs.
- Tasks:
  - [ ] Add explicit `--model-backend {modeller,diffusion}` and compatible
    diffusion options without changing the default.
  - [ ] Prohibit automatic fallback to MODELLER or another diffusion profile.
  - [ ] Treat the validated selected candidate inside the diffusion bundle as the
    input to `prepare` while preserving the complete bundle as provenance.
  - [ ] Merge `.dat` only through `DatRecord`.
  - [ ] Keep all intermediates beside the final output under batch operation.
  - [ ] Preserve final numbering and postflight diagnose order.
  - [ ] Define restart behavior that never trusts an incomplete or unvalidated
    existing bundle.
  - [ ] Use `--no-solvent` for development and acceptance iterations unless a
    separately requested solvent evaluation is running.

## WP12: Mosaic-First Homology Diffusion

- Suggested owner: Homology agent; coordinate shared contract changes through
  the integration coordinator.
- Owned surfaces:
  - `src/dvbfixer/homology_plan.py`;
  - Homology backend orchestration;
  - Homology fixtures and focused tests.
- Existing baseline:
  - [x] `homology --backend diffusion` and the GUI Model workflow use
    `selected_template_mosaic.pdb` as the authoritative frame.
  - [x] Companion coverage records zero-based half-open masks, precedence,
    template ownership, and fixed covered atoms in provenance.
  - [x] A one-chain five-residue internal insertion passed a checkpoint-backed
    MPS smoke without fixed-coordinate drift.
- Contract tasks:
  - [x] Keep `selected_template_mosaic.pdb` as the single authoritative
    coordinate frame.
  - [x] Export versioned coverage metadata from `materialize_template_plan`:
    target chain, template ownership, zero-based half-open spans, fixed/movable/
    generated masks, source provenance, and frame transforms.
  - [x] Never reconstruct template ownership later from the serialized PDB.
  - [ ] Feed supported neighboring template chains through the expanded WP5
    fixed-conditioning contract and prove causal context use.
- Initial slice:
  - [x] Support one target chain and one internal uncovered insertion.
  - [x] Generate uncovered insertion and terminal residues while covered
    sequence mismatches remain fail closed.
  - [x] Keep every covered template atom fixed and validate mosaic adherence.
  - [ ] Compare diffusion and MODELLER from the identical target and template
    plan.
- Expansion sequence:
  - [ ] Generate template-covered substitutions through explicit bounded
    junction windows.
  - [ ] Multiple uncovered spans.
  - [ ] Substitutions and deletions.
  - [ ] Multiple templates already fitted into the shared reference frame.
  - [ ] Multichain targets and interfaces.
  - [ ] Antibody H/L targets with distinct PDB chains, insertion codes, and
    shared interface context; never collapse both logical chains to `V`.
- Promotion gate:
  - [ ] Do not expose public Homology diffusion until mosaic adherence,
    multichain identity, interface context, geometry, ranking, and atomic
    publication gates pass.

## WP13: Final Verification And Decision Records

- Suggested owner: integration coordinator with independent scientific review.
- Tasks:
  - [ ] Run focused tests for every completed package.
  - [ ] Run the full CPU diffusion contract, scope, runner, validation,
    publication, provenance, benchmark, batch, ZBS, and Homology suites.
  - [ ] Run `pytest -m 'not slow' -q`.
  - [ ] Run `ruff check src/dvbfixer`.
  - [ ] Run strict mypy targets plus new typed diffusion entry points.
  - [ ] Run `python scripts/check_agent_docs.py`.
  - [ ] Run `python scripts/gen_cli_reference.py --check`.
  - [ ] Run `python scripts/gen_gui_spec.py --check`.
  - [ ] Run GUI typecheck/tests if registry or generated schema changes.
  - [ ] Repeat checkpoint-backed NVIDIA smoke and every approved outcome-bearing
    cohort in the pinned environment.
  - [ ] Record separate decisions for Linux profile acceptance, each scientific
    stratum, surrounding-aware capability, batch, ZBS, Homology, and any future
    production promotion.

## Global Stop/Go Gates

- **STOP:** thresholds, context shells, candidate counts, or comparison rules
  change after outcome-bearing results are inspected without a versioned new
  study.
- **STOP:** an adapter claims context conditioning without trace evidence that
  the model consumed the declared context.
- **STOP:** any fixed atom, chain identity, insertion code, explicit supported
  link, or authoritative mosaic span is silently changed or dropped.
- **STOP:** a runner falls back across device, engine, profile, context mode, or
  MODELLER without an explicit new invocation.
- **STOP:** an unsupported input launches inference or leaves partial public
  output.
- **STOP:** expanding gap length materially lowers preregistered hard-gate pass
  rate or ranking performance without an explicit unsupported boundary.
- **GO:** accept only the exact engine, environment, scope stratum, context
  capability, and workflow that passed its frozen gates.
- **GO does not mean:** diffusion becomes default, MODELLER is deprecated,
  arbitrary chemistry is supported, checkpoints may be redistributed, or an
  unevaluated profile inherits another profile's evidence.
