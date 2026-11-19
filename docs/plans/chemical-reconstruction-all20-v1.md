# WP2J: all20 explicit-state heavy/H coverage and exact-name heavy holdouts

Date: 2026-11-20; baseline `91ef103`, persistent
`feature/chemical-reconstruction`. Freeze before new source acquisition or controls.

## Scope and documentation plan

Exercise the existing shared reconstruction system over all twenty canonical
SOURCE residue types, in separate full-sidechain-heavy and complete-H arms.
This is bounded coverage/evidence, not a new geometry strategy, all20 capability
qualification, a combined heavy-then-H pipeline or a public/default replacement.
Do not conflate source residue coverage with all protonation/tautomer states.

Documentation: this protocol, new source fixture provenance, validation JSON,
active handoff and task map. Pending links: none. User's continuation authorizes
this bounded work on the persistent branch. No PyRosetta data or new native engine.

## Fixed source and state table

Source -> EXPLICIT requested state:

`ALA->ALA, ARG->ARG, ASN->ASN, ASP->ASH, CYS->CYS, GLN->GLN, GLU->GLH,
GLY->GLY, HIS->HIP, ILE->ILE, LEU->LEU, LYS->LYS, MET->MET, PHE->PHE,
PRO->PRO, SER->SER, THR->THR, TRP->TRP, TYR->TYR, VAL->VAL`.

ASH/GLH/HIP are caller-declared experimental states, NOT native-state/pH inference.
Use the existing variant owner, canonical RTP inventory transport and explicit
`PeptideStateDeclaration` route for EVERY state. All charge, H inventory and bond
orders come unchanged from the pinned CCD after explicit leaving-atom projection.
PRO removes source leaving H, all other targets remove source leaving H2;
carboxyl OXT/HXT and witness H2 are explicit as before. A mismatch refuses rather
than choosing another state, alias policy, source or target. Neutral HIS/HID/HIE
graph transformations remain outside this protocol.

Reuse original five amino-acid CCDs and WP2H's eight unchanged sources. Download
only ASP, GLU, LEU, MET, SER, TYR, VAL unmodified from
`https://files.rcsb.org/ligands/download/{ID}.cif` into separate
`tests/fixtures/chemical_reconstruction/ccd_all20_v1/`. Pin source URLs, hashes,
revision/atom/bond counts, CC0/usage policy and this pre-download protocol commit
in its manifest; commit source pins before numerical controls. Update the global
fixture manifest and explicit metadata allowlist. Do not edit old manifests.

## Eighty deterministic entries

Twenty source/state pairs x synthetic/observed x heavy/H = 80 declarations.
No random seeds, native engines, success-dependent selection or required pass count.

- Synthetic: pinned CCD ideal heavy target; same explicit virtual ALA chemistry
  and incoming-C construction as WP2I (PRO negative normalized two-heavy bisector;
  others heavy N/CA/C frame, Cpre-N-CA 121.96 degrees; N-Cpre 1.335 A).
- Observed: first file-order interior matching SOURCE residue in chain A of
  `tests/fixtures/protein_context/1UBQ.pdb`, except CYS/TRP use first interior
  matching source in chain A of `tests/fixtures/regressions/1FR2.pdb`.
  Never substitute a later residue or different file on refusal.
- Heavy: withhold ALL canonical target sidechain heavy atoms from the exact
  projected authority (every non-H/D atom except N/CA/C/O). Require all mask atoms
  present in observed source before withholding; no reduced mask to hide missing
  input or disorder. GLY has no sidechain-heavy mask: retain an explicit
  not-applicable entry, not an unchanged/reconstruction pass. Use existing CB
  stereo-cone plus deterministic chain compiler, including all ring closure and
  geometry audits; no ring solver/rotamer search or reference-driven torsions.
- H: complete selected H inventory with ALL selected heavy atoms observed, using
  the existing declared-state/extended/joint/source-metric compilers. This is NOT
  H generation conditioned on this protocol's reconstructed heavy output.

Reuse full source preparation, exact identity/links, opaque fractional-context
policy and independent source audit. Preserve original/prepared bytes and every
observation outside the explicit heavy holdout. Do not repair source planes,
strip context, choose altloc, rewrite occupancy or rename observed records.

## Distinct unchanged policies; no pooled acceptance statistic

Heavy arm uses existing WP2E `AuditPolicy(0.10 A, 12 degrees, 0.7 A, 0.02 A)` for
all measured bonds/angles, as its old plan does. H arm uses unchanged fixed-heavy
0.10 A / 12 degrees and generated-involving 0.02 A / 5 degrees, shared separation
0.7 A, H planes 0.005 A, heavy planes 0.05 A, stereo minimum 0.05. Do not describe
these arms as sharing an identical acceptance gate or pool them into a scientific
success rate. Audit raw and PDB-rounded additions; retain construction/geometry/
rounding/source/chemistry failures. Only the full-H arm has complete conditional
CIP; heavy arm uses existing explicit source stereo constraints, not fake CIP on
an H-incomplete graph. Neither establishes actual-neighbour/polymer qualification.

## Heavy holdout scoring and evidence isolation

Add a private exact-identity, unchanged-frame heavy holdout metric, no fitting or
symmetry relabelling. Require the reference identities and candidate additions to
equal the complete declared heavy mask exactly once, with matching elements,
angstrom units and unit occupancy. Report atom count, RMS displacement and maximum
displacement for complete raw/rounded candidates EVEN IF geometry rejected; an
incomplete candidate cannot be scored. Scores are descriptive, never acceptance
gates or constructor/ranking inputs. Do not change mask, phase or graph based on
the reference. No native-H RMSD or experimental H quality claim.

Read held-out coordinates only AFTER proposal completion in the driver/scorer;
neither request nor compiler/proposal interface contains them. Retain full saved
request/plan/results and reference evidence separately. Test reference leakage
boundaries and exact mask/identity/unit/duplicate/incomplete-score refusal.
Independently rehydrate final sources/requests/plans/candidates, measure geometry,
H conditional CIP and heavy scores without proposing coordinates. Replay old
controls byte-for-byte; compare known input repetitions as regressions, not new
independent science. Preserve WP2I's earlier PRO interface failures.

Run focused tests, split slow boundaries, Ruff, scoped/source-bound strict typing,
fixture integrity, generated CLI/GUI checks and full fast suite. Keep the ambient
legacy typing failure and unresolved WP2F combined slow timeout distinct. All20
production support, unseen-input/native-quality qualification, partial/measured H,
linked/metal/other chemistry, whole-gap modeling and production/public/default/
dependency/data-license gates remain OPEN regardless of the observed counts.
