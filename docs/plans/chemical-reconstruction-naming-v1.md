# WP2C: explicit inventory naming transport

Date: 2026-11-19. Prospective protocol; baseline `7802240`, same persistent
`feature/chemical-reconstruction` branch. Freeze before executing the new controls.

## Scope and documentation plan

Implement a private, pure inventory adapter using the existing domain
`atom_name_policy` with explicit `AMBER19SB` / `STANDARD` selection. Do not patch
the canonical naming tables, add residue-specific fallbacks, search graph
isomorphisms, infer state, mutate CCD chemistry or rename source PDB observations.
The public naming service and WP2B runner/default behavior stay unchanged.

Documentation: this protocol, a separate naming validation JSON, the active
handoff and existing agent task map. Pending links: none. User continuation
authorizes this bounded implementation, not a public reconstruction route.

Use the already pinned Amber19SB RTP and five existing CCD amino acids. A raw
`TemplateInventory` must carry its original exact atom/bond names and external
tokens. Transport records a complete immutable simultaneous raw-to-canonical
bijection, rule IDs, source pin, explicit dialect/profile and both inventories.
Reject duplicate destinations including implicit identity collisions, invalid
names/metadata, unknown local tokens and malformed external references. Transport
`-C`/`+N`-style references without dropping the sign or treating them as local
atoms. Preserve order and all rows, including external impropers/CMAP terms.
No sequential chained renaming and no opportunistic second policy after failure.

The naming record alone proves no chemistry. Feed the transported inventory to
the existing explicit peptide substitution/witness proof. Exact retained atom and
bond inventory must agree. RTP remains connectivity/inventory evidence only:
it does not provide formal charge, bond order, atom stereo or H geometry. The
CCD source retains those chemical annotations unchanged. Equivalent-H name
transport is a deterministic naming convention, NOT independent qualification
of stereospecific proton coordinates or pro-R/pro-S assignments.

## Declared controls

Five positive hypotheses: raw RTP ALA/THR/ILE/LYS/PHE each transported with the
one declared canonical policy, then joined chemistry/inventory proof. ILE is
expected to pass now because this existing policy recognizes the CD/HD aliases;
retain any contrary result rather than modifying policy after inspection.
Sixth control: raw RTP LYN with CCD LYS must still refuse unsupported explicit
variant projection. No pH/charge edit or state-name substitution is allowed.

Record six entries, exact input pins, full transport/authority evidence, every
refusal and final source hashes. Source pins and the protocol hash must be
checked before output reservation; write a new non-overwriting project evidence
directory. Use original CCD coordinates only as chemical authority, not held-out
experimental coordinates; no coordinate proposals/sampling/geometry trials.
WP2B ILE refusal stays an unchanged regression and historical result. Do not rerun
or repair WP2B observed structures as scientific controls in this slice.

Adversarial regressions: incomplete/duplicate/colliding aliases, unknown local
tokens, malformed external tokens, simultaneous H-name transport, source-pin
binding, mapped connectivity mismatch, input immutability, wrong state, and
unchanged WP2B refusal. Do not accept negative controls by weakening gates.

## Acceptance and remaining gates

Run focused reconstruction/naming tests, Ruff and private scoped strict mypy;
run the full fast suite when available. Preserve transitive legacy typing failures
separately. Verify old fixture/force-field bytes and protocol immutability, then
commit coherent implementation/evidence metadata on the same persistent branch.
This slice does not qualify whole-polymer state, actual neighbour chemistry,
observed reconstruction, H placement, source preparation, default/dependency
replacement or public/production/redistribution permission. External-anchor
amide-H geometry and prospectively specified observed-source admission remain
the next independent work items.
