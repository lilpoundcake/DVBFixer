# WP2A: pinned CCD authority and prospective geometry roundtrips

Date: 2026-11-19. Status: protocol frozen before component download/execution.
Persistent branch: `feature/chemical-reconstruction`, baseline `49e3568`.
User authorized continuation; no production/default integration or native sampling.

Documentation plan: this protocol/decision, active handoff and agent-task updates,
fixture provenance/checksums and a validation JSON. Pending links: none. Bounded
implementation, source assets and tests accompany these records.

## Source and reuse

Reuse `ccd.load_ccd_component` / `parse_ccd_component` and their exact SHA pins,
not another CIF parser or downloader. Acquisition is an explicit development step
using the existing fixed-endpoint `fetch_ccd_component`, separate from admission,
generation and the new offline authority adapter. Preserve raw downloaded CIF
bytes; pin every file before controls run. No runtime network fallback.

Official source: https://files.rcsb.org/ligands/download/COMPONENT.cif.
Dictionary description: https://www.wwpdb.org/data/ccd.
Usage policy: https://www.wwpdb.org/about/usage-policies (CC0 1.0 data dedication;
checked 2026-11-19). Preserve attribution and source/modified dates. This policy
does not relax PyRosetta/Dunbrack source/version/license or production gates.

Use **only CCD ideal coordinates**, explicit bond orders/charges and all declared
hydrogens as the supplied microspecies. Never import experimental/model coordinates
as geometric targets. Retain every leaving atom; free CCD amino acids are NOT
peptide residues. No automatic removal of OXT/HXT/terminal H, charge edits, pH
selection or linking projection. No claim of polymer-sidechain/H qualification.

Extend the existing heavy-graph converter through a shared implementation with a
separately named all-atom function; preserve old heavy-only behavior. The new
offline adapter verifies the pin before loading and validates full graph chemistry
with RDKit: explicit H/no implicit additions, sanitization, no unsupported radicals,
and 3D-assigned CIP versus each CCD R/S label. Refuse unsupported bond stereo,
unresolved atom stereo, missing coordinates/charges or ambiguous chemistry instead
of silently omitting information absent from the current shared graph DTO.

Only after CIP agreement may the adapter derive an ordered geometric stereo frame
from the source ideal coordinates. Deriving a sign from coordinates alone is not
a check that a label is chemically correct. For heavy-only operations use admitted
heavy-neighbour frames; for complete H controls also require four-ligand stereo
and declared planar sp2/ring constraints where representable. This is an explicit
library/version-bound authority verification, not an independent CIP implementation.

## Plan compiler

Compile deterministic chain placements from graph adjacency and the fixed observed
mask; consult source ideal geometry only to exclude degenerate ordered references.
Never use deleted experimental coordinates, proposal scores or known native RMSD.
Keep the caller's explicit audit bounds and whole retained context. A compiler
refusal does not trigger another engine or change the frozen controls.

The first compiler emits existing chain placements only. Planar tips and new ring
closure solving are NOT added merely to pass this cohort. Ring-containing ideal
roundtrips can pass existing closure checks; that is not a general closure solver.

## Frozen controls (before download)

Components: **ALA, THR, ILE, LYS, PHE, EOH, ACT**. Each runs heavy and all-H arms
under two proper rigid transforms: identity and `(x,y,z) -> (z+7,x-4,y+2)`.
Total: **28 declared controls**. Keep every authority/admission/compiler/geometry
refusal. Do not substitute components, change masks or widen bounds after outcomes.

Heavy masks:

| Component | Missing atoms |
|---|---|
| ALA | CB |
| THR | CB, OG1, CG2 |
| ILE | CB, CG1, CG2, CD1 |
| LYS | CB, CG, CD, CE, NZ |
| PHE | CB, CG, CD1, CD2, CE1, CE2, CZ |
| EOH | Its heavy oxygen (expected insufficient three-heavy-anchor support) |
| ACT | Lexicographically last heavy oxygen (expected unsupported planar branch in the chain-only compiler) |

Heavy observations are all remaining heavy ideal atoms, including free-molecule
leaving groups, transformed by the declared transform. No H is measured in that arm.
The H arm retains every heavy ideal atom and requests every CCD H; no H is deleted
from a real protein. Source identity/masks are exact and case-sensitive.

Bounds declared now: bond tolerance **0.02 Å**, angle tolerance **2°**, severe
minimum separation **0.7 Å**, rigid-anchor RMSD bound **0.02 Å** (unused by chain
construction), stereo normalized minimum volume **0.05**, planar tolerance
**0.02 Å**. These are roundtrip-control bounds, NOT calibrated defaults for users.
Require both raw and generated-only three-decimal rounded audits. Re-audit complete
H-arm candidate CIP with RDKit without calling the constructor. Report transport
error to the same CCD authority explicitly as **authority roundtrip error**, not
held-out/native prediction accuracy. Do not rank/tune candidates using that error.

## Interpretation and next gates

These controls deliberately share the independently sourced CCD ideal authority
with their transformed observations. Success establishes source conversion,
stereo mapping, exact-mask planning and real chemical geometry transport only.
It is NOT evidence of native conformer recovery, general polymer valence/H counts,
protonation, distorted-source closure, complex packing or superiority to templates.
Report all outcomes before any next protocol. Real observed complexes, protein
polymer state authority, mixed generated ownership, broader whole-sidechain/H and
source-preserving publication remain distinct future gates. Public APIs/CLI/GUI,
production defaults and dependency retirement remain disabled.
