# WP2I: explicitly declared peptide leaving atoms and unchanged-graph HIP naming

Date: 2026-11-20; baseline `1d52019`, persistent
`feature/chemical-reconstruction`. Freeze before any new numerical controls.

## Scope, documentation and authority

Add a separate private declaration/derivation route; preserve all old defaults,
refusals, protocols and results. Reuse existing pinned CCD sources (PRO and HIS
in the WP2H manifest, ALA in the original manifest), canonical inventory transport,
the existing `ffutils.variants.VARIANT_TO_PARENT` naming owner, chemical-witness
assembler, source preparation and independent validators. No download or new
licensed data, production/public/default/dependency capability follows.

Documentation plan: this protocol, validation JSON, active handoff and task map.
Pending links: none. The user's continuation authorizes this bounded slice.

Declare exact target identity, requested state, parent and ALA-witness digests,
canonical inventory digest, and the complete two substitutions. The derivation
must verify every binding. All removals must already be marked source leaving
atoms and pass the existing assembler's branch, port, valence and CIP checks.
Require directed left-ALA -> target -> right-ALA connectivity. No optional default
leaving atom, source-neighbour substitution, chemistry inference or automatic
fallback. The virtual ALA units prove chemistry only and never enter output.

For PRO explicitly remove its source-marked N-H on the incoming bond, not absent
H2. Carboxyl removals remain source OXT/HXT and the outgoing ALA removes H2.
For explicitly requested HIP, use the existing variant owner to bind HIP -> HIS,
then PROVE exact HIP RTP atom/connectivity inventory against the unchanged
cationic CCD HIS graph after peptide leaving-atom removal. No atom additions,
formal-charge edits, bond-order edits, tautomer choice or pH inference are allowed.
Calling a neutral state HID/HIE cannot silently delete a sidechain H or move charge;
those declarations must refuse on exact inventory mismatch. Naming is not itself
chemical authority: graph sanitization, explicit H inventory, source stereo and
conditional witness proof all remain required.

Use a distinct declaration/authority digest and provenance schema, keeping
`LocalPeptideAuthority` compatibility only through explicit subclassing. Old
`derive_local_peptide` must still reject old PRO/HIS/HIP cases as before. Reuse
the existing full-complex request/source transport; do not truncate context or
rewrite observed residue names, coordinates, occupancy or records to say HIP.

## H construction boundary

Use the WP2H extended compiler plus WP2G joint-carbon groups unchanged. For a
declared peptide N with no local H and two observed heavy neighbours (PRO),
explicitly require no backbone-H typing, add a fixed heavy plane containing that
N, the real incoming C and both local heavy neighbours, at the existing 0.05-A
heavy-plane tolerance. Bind the incoming observed covalent peptide connection;
do not invent an N-H atom/type. If a backbone H exists (HIP), require the exact
existing N/H/CA parameter types and reuse pinned ff19SB amide geometry. This does
not establish all actual-neighbour/full-polymer heavy-angle chemistry.

All observations remain fixed. H construction is complete-selected-H only,
without measured target H/D, phase search, source repair, ring-heavy rebuilding,
reflection or numerical optimization. Retain WP2H's separate policy and bind the
state declaration additionally in the new compiler origin.

## Eight declared entries

1. Synthetic complete-heavy PRO -> explicitly requested PRO, all seven H missing.
2. Synthetic complete-heavy HIS -> explicitly requested HIP, complete H mask.
3. First interior PRO in chain A of original 1UBQ -> explicitly requested PRO.
4. First interior HIS in chain A of original 1UBQ -> explicitly requested HIP.
5. Chemical-proof-only HIS -> HID, unchanged source graph, exact HID inventory.
6. Chemical-proof-only HIS -> HIE, unchanged source graph, exact HIE inventory.
7. Chemical-proof-only PRO declaration removing absent H2 rather than leaving H.
8. Chemical-proof-only declaration with a mismatched source-authority digest.

Synthetic heavy target coordinates are pinned CCD ideals. For HIP use the WP2H
incoming carbon construction (1.335 A, Cpre-N-CA angle 121.96 degrees, heavy
N/CA/C frame); for PRO use 1.335 A along the negative normalized sum of its two
normalized N->heavy-neighbour directions. This places the new external anchor in
the observed-heavy N plane without consulting any generated H. Synthetic anchors
are declared inputs, not edits to observed PDB sources. Do not substitute another
observed target or alter the synthetic source after outcomes. No pass count is
required; retain every chemical, fixed-input, geometry or CIP refusal.

## Unchanged evaluation and gates

Fixed-heavy 0.10 A / 12 degrees, generated 0.02 A / 5 degrees, separation 0.7 A,
H planes 0.005 A, heavy planes 0.05 A, stereo minimum 0.05. Audit raw/rounded
independently with geometry plus conditional CIP. CIP alone does not validate H
positions/names. Independently read back saved sources/declarations/requests/plans/
candidates and remeasure without constructors; replay previous controls byte-for-
byte, not as new scientific evidence. Test exact identity/state/digest/role binding,
source leaving flags, neutral-state refusal, PRO no-NH inventory/typing and external
plane, declared HIP explicit H/count/charge, corrupt candidates, rigid equivariance,
context preservation and old-default rejection. Run focused tests, split slow
boundaries, Ruff, scoped/source-bound strict typing and full fast suite.

Broader all20 heavy/H coverage, neutral histidine transformations, measured/partial
H, actual-neighbour/polymer chemistry, native-H/name/orientation quality, linked/
ring-heavy reconstruction and broad retained-complex qualification stay open.
WP2F's unresolved combined slow timeout stays unresolved. No production or
PyRosetta redistribution/license permission is inferred from this slice.
