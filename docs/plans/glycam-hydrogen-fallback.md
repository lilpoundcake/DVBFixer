# Authorized GLYCAM fallback connectivity repair

## Authorization, baseline and documentation scope

The user approved a separate fallback fix after the
[real-input qualification](../research/oh-native-qualification.md). Branch:
`fix/glycam-charmm-hydrogen-fallback`, starting from `eb049aa`. No merge into main,
release, empirical-library integration or native FF retuning is authorized.

Documentation plan: this prospective implementation/validation contract, a
scalar-only qualification record, a new current-source inventory amendment, and
updates to the preparation task/contract/invariant maps and canonical handoff.
Keep all old sources, failed outputs, hashes, protocols and cohort status intact.
Add links only after their targets exist; validate links and public payloads.

## Source-bound observation (not a new native run)

Original case 11, GLYCAM hinge/CH3 with CHARMM request, remains a failed member of
the old 12-case cohort. Its saved `h-2-before` state contains unique H1 atoms in
`H/UYB232` and `K/UYB242`, each unbonded. At the H-chain site, C1–H1, C1–O5 and
C1–C2 distances are 1.090961/1.467420/1.566966 Å, respectively; all three bonds
are absent before fallback. The successful GLYCAM path supplies template bonds;
the CHARMM path lacks GLYCAM templates, then rebuilds from `fixer` for fallback.
OpenMM recognizes existing H by **bonded parent**, not merely by residue/name,
and adds a second H1 when the original is unbonded. Postflight makes no change.

The upstream algorithm is inspected at the runtime's exact OpenMM revision
`36a30cbca54e727b216b606f3c011b67201eb8b4`, not inferred from final naming alone.
This observation does not establish why the original CONECT/repair transport
omitted those bonds; qualify the bounded restoration independently.

## Proposed implementation contract (before production edits)

- Apply only in legacy prepare's existing whole-FF-error fallback, before its
  forcefield-free `Modeller.addHydrogens`. Successful whole-FF paths stay intact.
- Identify canonical GLYCAM residues using the shared registry, by native
  residue object identity. Only residues with existing unbonded H need repair.
  Reject duplicate atom names in a selected residue instead of collapsing them.
- Use a private, isolated topology view of affected residues and the existing
  `build_glycam_system` wrapper with the shipped GLYCAM XML. Use its canonical
  exact-name templates; do not guess by coordinates or hard-code UYB/H1.
- Harvest only missing **intra-residue** template bonds. Never add/remove an
  external bond, allocate/delete/rename an input atom, move any coordinate,
  reclassify protonation, change requested FF, or use this template FF for
  native minimization. Fallback remains explicitly forcefield-free.
- Validate all affected residues before applying any additions. Unknown
  templates, missing/ambiguous parents, inconsistent elements or incompatible
  existing H bonds must fail closed with a typed topology error. Do not strip
  or delete an H merely to hide a duplicate, and never invoke another fallback
  after that error.
- The O–H postflight contract, PROPKA/Reduce, disulfides, `keepIds=True`, original
  inputs, variant identity, `.dat` schema, native forces and minimizer settings
  remain unchanged. Unknown non-GLYCAM ligands retain the old route.

## Prospective validation and native budget

1. Keep the original failed cohort and its false complete-cohort status. Bind
   the exact saved case-11 pre-fallback state and native artifacts; do not
   overwrite or adopt any historical output directory.
2. Before implementation, add a deterministic minimal unbonded UYB-H1 regression
   and confirm it fails on baseline. Test already-bonded idempotence, two chain
   cases/insertion codes/occurrences, duplicate names, missing parent, wrong H
   parent/element, unknown template, unrelated ligand and all-or-nothing failure.
   Native-free tests must verify exact identity/coordinate/external-bond
   preservation and template provenance; an unbound H must not remain silently.
3. On the actual saved pre-fallback state, restore topology once without H
   sampling. Expect the two C1–H1 and four C1–ring links to be the only additions;
   verify immutable atoms/positions and all inherited links, then idempotence.
   Retain any different result as a blocker, not an invitation to broaden scope.
4. Freeze a new source/runtime/input-bound execution plan after engineering
   tests, before native execution. Enroll the same six tracked sources and two
   FF requests, same table order, pH, water/naming/defaults and per-H seed label
   1 as the old cohort: **12 new preparations, one per case**. This is a changed
   candidate's prospective cohort, not an unchanged retry or replacement of
   the old failure. Preserve native FF vs forcefield-free route accounting.
5. Keep **1,200 s/case, 7,200 s aggregate, eight H calls/eight native H calls
   per case**. Retain failures/timeouts/unexecuted cases; no replacements, new
   seeds or outcome-dependent repeats. Observer snapshots must now serialize
   chain/residue occurrence indices explicitly. Archive the original frozen
   observer unchanged; no silent source/hash rebinding of old evidence.
6. Require unique final atom identity, valid `DatRecord`, no new O–H contract
   violations, exact matched repair replay, immutable heavy/valid H and inherited
   external links, plus bounded actual-route coverage. Partial PDB existence is
   never success. Native CHARMM glycan and direct ligand-OH qualification remain
   separate claims; a forcefield-free success does not establish either.
7. Run focused tests, Ruff, strict-target mypy, generated CLI/GUI checks,
   agent-map checks and one complete fast suite after source/inventory changes.
   Preserve every failure; do not weaken historical assertions or add broad skips.

Stop on unexplained identity/movement/graph changes, a new refusal or missing
source evidence. Any scientific amendment must be recorded before further
execution. All private coordinates/native/RNG artifacts stay below the ignored
`tmp/hydrogen-postflight-contract-01/glycam-fallback-01/`; tracked records carry
only identities, hashes and scalar measurements.
