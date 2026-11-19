# WP2B: explicit local peptide chemistry and observed-source probes

Date: 2026-11-19. Status: frozen prospective protocol, before projection/control
execution. Baseline `5cd21aa`; persistent `feature/chemical-reconstruction`.

Documentation plan: this protocol, handoff/task-map updates and a validation JSON;
source data and prior protocols remain unchanged. Pending links: none. No new
production route, dependency replacement, optimizer or native sampling cohort.

## Bounded chemistry, not a disguised whole-protein authority

Use pinned existing CCD ALA/THR/ILE/LYS/PHE and the existing Amber19SB RTP source
(`ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419`).
Reuse its parser and the existing canonical force-field atom-name maps. The RTP
supplies only an independent exact atom/bond inventory; partial charges MUST NOT
be interpreted as formal charges or protonation decisions. Retain its upstream
provenance/license notices. No PyRosetta/Dunbrack table access or new data download.

Explicit substitutions join carbonyl C to N with a single bond. Each attachment
declares its exact removed atoms. They must be source leaving atoms connected only
through the named port. For these ordinary CCD amino acids the declared two links
remove left OXT/HXT and right H2; do NOT delete every leaving atom globally.
Preserve all parent graphs, deletions, connections and name mappings as provenance.

A local internal-residue proof uses **two virtual neutral CCD ALA partners**.
The complete explicitly joined three-unit chemical graph must pass RDKit closed-
shell valence, charge/no-implicit-H preservation and 3D CIP comparison. Source
coordinates remain in independent local CCD frames: inter-unit geometry is NOT
validated, no new chemical bond is assigned a fabricated ideal length, and witness
atoms are NEVER observations, generated output or `.dat` atoms. The proof is
conditional local chemistry, not the actual identity/chemistry of neighbouring
source residues or global protein CIP/charge/parameterization qualification.

The derived target has both peptide ports consumed and must match the explicitly
selected RTP inventory via an exact bijective canonical-name map. Do not infer a
variant from pH, sum partial charges into formal charges, or silently swap LYS to
LYN. A different source/state inventory or an unimplemented variant refuses.

**H geometry is not licensed by H inventory.** Retained free-CCD backbone N–H
geometry is not an amide geometry authority. All H-placement operations through
this local projection refuse explicitly for now; target measured H also refuses.
The projected all-H chemistry inventory is useful progress, NOT a qualified H
coordinate generator. Backbone-H geometry, full polymer state and linked assembly
qualification remain separate gates.

## Prospective controls

Chemistry probes: five declared source/template pairs ALA/ALA, THR/THR, ILE/ILE,
LYS/LYS and PHE/PHE, plus LYS/LYN. Use ALA as both chemical witnesses. LYS/LYS
is expected to fail the independent charged-LYS inventory; LYS/LYN is expected to
refuse unsupported variant projection. Retain every result, do not rename a state
or remap hydrogens opportunistically to make it pass. Six declared chemistry probes.

Observed sources already tracked (no source-coordinate inspection for selection):

| Source | Target selection | Heavy mask |
|---|---|---|
| `protein_context/1UBQ.pdb` | exact chain A PHE45 | CB CG CD1 CD2 CE1 CE2 CZ |
| `regressions/1FR2.pdb` | first file-order interior THR in chain A | CB OG1 CG2 |
| `protein_context/1A3N.pdb` | first file-order interior ALA in chain A | CB |

Each has a heavy arm and an H-refusal arm: six declared observed probes. Select
the target from residue identity/file order only, never by geometric success;
selection/missing-reference/occupancy/link/state/geometry refusals remain counted.
Do not crop 1UBQ's partial-occupancy context or strip 1A3N hemes to pass admission.
The H arm must refuse as unqualified, not claim a generated scientific candidate.

Reuse the existing PDB atom parser and declared-link parser. Bind original source
bytes and an explicit withhold mask; only the resulting retained observation DTOs
reach the constructor. No withheld experimental coordinate becomes geometry
authority, no native-RMSD ranking and no held-out accuracy score in this slice.
Retain every other chain, heterogen, water, link and unknown-chemistry diagnostic.
Bind observed predecessor-C/target-N and target-C/successor-N using the existing
`observed_peptide_continuity` owner, without guessing an order for other LINKs.
No phantom witness atom is inserted into the source. No PDB file is materialized.

The observed-heavy policy is prospectively **0.10 Å bond / 12° angle / 0.7 Å
minimum-separation / 0.02 Å anchor bound** (anchor bound unused by chain strategy),
normalized stereo minimum **0.05**, planar tolerance **0.05 Å**. This is a new
explicit observed-probe policy, not a change to the frozen CCD roundtrip policy
or calibrated public defaults. Require both raw and generated-only rounded audits;
refuse rather than move observations, widen bounds or try a second engine.

## Interpretation

Twelve declared probes total. Chemical proof, observed-source admission,
coordinate proposal and geometry qualification are distinct denominators. Unknown
source chemistry remains unassessed; a local witness proof cannot establish
full-source chemistry. All MD/publication/global-source-readiness flags remain
false. No all20, distorted closure, conformer superiority, observed H or general
retained-complex qualification follows. Report every refusal before the next
protocol, and continue heavy/H state and geometry work together.
