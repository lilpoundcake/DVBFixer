# WP2O: source-bound state-invariant incoming carbonyl type authority

Date: 2026-11-20; baseline `502a8d9`; persistent
`feature/chemical-reconstruction`. Freeze before new type-bound port controls.

## Scope and documentation plan

Prove ONLY a real incoming carbonyl-C force-field type, without selecting an
incoming protonation/tautomer/disulfide state or assigning any partial/formal
charge, H inventory, probabilities or full-residue chemistry. Reuse the existing
AMBER variant owner, RTP parser, pinned CCD authorities, exact source/probe
transport, port constraints, projector and projected-H auditor. Do not add a
parallel histidine table or choose HID/HIE/HIP to bypass the old HIS refusal.

Documentation: this protocol, validation JSON, active handoff, task discovery and
local quality note. Pending links: none. User approved starting the incoming-type
authority work. No production/default/dependency/public/DatRecord change, new
licensed data or broad native/H quality claim. Altloc selection is NOT this slice.

## Pinned complete family authority, not a convenient subset

RTP source remains `FF/amber19sb.ff/aminoacids.rtp`, SHA256
`ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419`.
The canonical residue set is the existing CANONICAL_HEAVY_ATOMS owner; variant
relationships come ONLY from `ffutils.variants.AMBER_VARIANT_TO_PARENT` (namespace
is explicitly ff19SB, not CHARMM). Freeze its sorted mapping digest
`f92fc4f12e8997e14efbd407703e1df34460dbbaa78f34e1ac137891ad678806`.

Reuse parse_rtp and freeze all 27 present canonical/AMBER-variant nonterminal
templates, not just HIS or the C row. Snapshot each name, sorted full atom rows
(name/type/charge/group), sorted normalized bonds and sorted improper rows.
Tuple-payload digest is
`2f34f225fadd62b052ad1d8d74421af05cc94e7adace041c01fd9d8ba808b940`.
Unknown/terminal/other-namespace labels, changed catalog/registry, omitted state,
duplicate atom/bond, missing/malformed local motif or incompatible types refuse.

For an explicit AMBER variant label, inspect that exact template (do not rewrite
the label). For a canonical parent label, inspect ALL registry-mapped variant
templates plus the parent template if present in the PINNED complete catalog.
Missing required variant is refusal, never subset proof. Thus unresolved HIS
inspects HID/HIE/HIP without selecting one. Standard HIS has no exact RTP entry
in this frozen catalog; that absence is explicit data, not a guessed default.

Every inspected template must contain the internal carbonyl C with exactly local
C–CA and C–O bonds, and the outgoing-peptide improper `CA +N C O`. Require a single
common C atom type, and accept only exact type C for the current scalar port rows.
RTP lacks bond orders: obtain C–CA single / C=O double and the OXT/HXT leaving
branch from the pinned canonical CCD parent's local carbonyl footprint; verify
neutral nonaromatic retained C/O and the leaving flags. Use ONLY this footprint,
not CCD's whole-residue charge/H state. Actual outgoing C–N substitution comes
from the real observed directed peptide link, never a virtual-neighbour position.

Charges deliberately are NOT required to agree or transferred: HID/HIE C has
0.5973 while HIP C has 0.7341 in the pinned RTP. The proof records candidate rows,
but selected_state and assigned_charge are ALWAYS None. Sidechain/H differences
remain unassigned; donor Component.states/selected_state remain untouched/opaque.

## Source admission and proof object

Add a distinct source-binding record and pure type proof, not a new ChemicalState.
Bind source/prepared SHA, original request digest, complete observation digest,
exact incoming C and recipient N identities, literal source residue label, and
observed incoming C/O/CA plus recipient N. These must be fixed, unowned and unit
occupancy; retain chain case/insertion codes. Donor must be unselected opaque
context, with no source OXT/HXT and an exact order-1 covalent C–recipient-N link.
Refuse extra/unknown C attachments contradicting the projected carbonyl motif.

Source adapter reuses prepare_source_probe/audit_source_transport and reparses
through the EXISTING PDB adapter for the actual label. Rebind complete probe/
connection transport, including unassessed same-component CONECT endpoints;
do not implement an independent bond-perception/parser or infer from proximity.
Output only observed core records and metadata, never withheld target reference
coordinates or the full source text as construction input. A separate source-byte
audit is mandatory before science and on saved readback; proof shape/hash alone
does not establish source provenance. Catalog/CCD files are pinned/reloaded by
the evidence I/O adapter. No incoming coordinates/state/charges are changed.

## Explicit opt-in port v2 and unchanged physical stages

Use a separately tagged type-proven planar-port declaration carrying the proof.
Pure port binding must check its source request, incoming/recipient identities,
catalog/CCD/variant-family data and common type before using the SAME scalar
constraints. Existing v1 declaration fields/digests, driver HIS refusal, old parent
objects and default behavior remain unchanged. Do not silently upgrade old records.
Existing projected-parent H lineage/audits must also revalidate the v2 proof.
Extend the private strict known-DTO registry only; no public schema/API capability.

No change to source preparation, geometry objectives, starts, movement caps,
solver settings, masks/states, bond/angle/plane/clash/stereo policies, H compilers,
H names or conditional CIP scope. One unchanged old baseline, one new constraint-
bound heavy projection and (ONLY for passing projection) one H proposal. No H-
dependent heavy movement/retry/search. Native references/scores extracted ONLY
after both proposal phases finish, never in proof/solver/H inputs.

## Frozen narrow continuation and tests

New v2-policy declaration uses the SAME first file-order interior PRO chain A in
full retained 1UBQ, 1FR2 and 1A3N sources from WP2M. Whole CB/CG/CD mask and exact
explicit PRO state/leaving-H authority remain. 1UBQ is regression/known-case;
1FR2's formerly refused incoming-HIS context is prospective for this new type
policy, NOT a new independent protein family. 1A3N altloc must remain refused,
not silently selected. There is no required projection/H pass count or native-RMSD
improvement. All source/proof/baseline/projection/H refusals remain in denominator 3.
Old WP2M/N refusals and native-score worsenings stay immutable. A larger independent
structure/missing-mask cohort is a subsequent separately frozen qualification gate.

Tests use tracked/synthetic inputs: complete family/explicit-label behavior,
conflicting/missing types/topology, differing charges without assignment, changed
owner/catalog/CCD/source/label/identity/masks/links, donor carbonyl attachments,
forged proofs, no tautomer choice/whole-state inference, strict DTO roundtrip and
constructor-free audits. Independently read back saved source proofs/coordinates/
metrics/H/CIP without constructors/solver, replay old controls byte-identically,
and run focused/split boundaries/full fast/Ruff/scoped/source-bound typing/maps/
generated/fixture checks. Preserve ambient legacy typing failure and WP2F timeout.
