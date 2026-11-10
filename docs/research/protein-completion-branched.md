# Typed-kernel VAL/THR paired-branch experiment

Status: **private typed-kernel expansion only**. Production, source/cohort
materialization and private publication target capability remain ALA/SER/free CYS.
The new cases cannot silently enter those downstream consumers.

Continuation: a [separate explicit VAL/THR source adapter](protein-completion-branched-source.md)
now implements global THR-Cβ and rounded readback. Default source/cohort/evidence
entry points remain restricted; the typed-only status and hashes below describe
the earlier checkpoint, not an overwritten execution record.

A [retrospective matched seeded-PDBFixer comparison](protein-completion-template-comparison.md)
has now executed on these and singleton masks. Its measured-chi common gate,
conditional accuracy results and stereo failures are distinct from these earlier
own-only controls; they do not establish general or full-legacy superiority.

## What is implemented

- `protein_completion/branched.py`: explicit paired-branch geometry and an
  independent named-Cβ signed-volume measurement.
- `complete_sidechain`: three staggered χ1 proposals for VAL CG1/CG2 and THR
  OG1/CG2, preserving all observations. CB can be retained or generated from
  observed N/CA/C. No template coordinates, MD optimization or downloads.
- `audit_candidate`: both branch bonds/angles, primary χ1, named-branch
  orientation, Cα chirality and generated-atom clashes. Reflecting the second
  branch while keeping Cα and primary χ1 unchanged fails the new orientation gate.
- Complete observed residues are independently audited on a no-op; wrong fixed
  THR Cβ is rejected, not reflected or normalized.

The admitted masks are exactly `(CG1,CG2)` or `(CB,CG1,CG2)` for VAL and
`(OG1,CG2)` or `(CB,OG1,CG2)` for THR. A single missing distal branch with an
observed partner **refuses** pending an observed-anchor solver. Missing CB with
retained distal atoms still refuses closure. Existing hydrogens, backbone defects,
linked targets and unsupported chemistry keep their prior refusals.

Illustrative geometry: C–C 1.53 Å, C–O 1.43 Å, ideal tetrahedral Cβ branch angles
109.4712°, with second-branch χ offsets +120° (VAL) and -120° (THR). These are
explicit control values, **not an empirical backbone-dependent geometry library**.
The result revision is `canonical-branched-sidechain-staggered-v1`. The optional
SER/CYH authority is explicitly refused for these new residues, never substituted.

## Physical stereochemistry versus atom naming

For vectors from CB, measure `(CA × first_branch) · second_branch`:

- **THR:** positive volume corresponds to natural Cβ **R**; L-threonine is
  **(2S,3R)**. A negative value produces the allo configuration while Cα can
  remain correctly S. Checking only Cα cannot catch this defect.
- **VAL:** CG1 and CG2 are identical methyl groups, so CB is **not a molecular
  stereocentre**. Negative named volume enforces the tested CG1/CG2 orientation,
  not a fictitious R/S configuration. Swapping names is not evidence of a D amino
  acid; exact-name scoring nevertheless distinguishes those labels.

`scripts/branched_stereo_oracle.py` independently assigns RDKit CIP labels from a
local amino-acid graph plus proposal/reference coordinates. Atom-map labels are
removed before assignment so coordinate transport cannot change ligand ranking.
An auxiliary terminal oxygen supplies a graph/conformer placeholder only: it is
never given to reconstruction, scored as generated or published. This oracle is
research-only; core validation has no RDKit dependency.

Independent controls establish:

- Eleven complete 1UBQ reference residues (four VAL, seven THR): CA=S for all,
  CB=R only for THR; VAL CB receives **no CIP label**.
- Every synthetic χ1 proposal has the same expected physical stereochemistry.
- Reflecting THR's second branch changes CB from R to S while CA stays S and
  primary χ remains unchanged. Both RDKit and the measurement-only gate detect it.

This checks target THR Cβ, **not all source/ligand/ILE stereocentres**. The source
adapter is intentionally not expanded to THR/VAL targets until separate source,
rounded-output, context-stereo and publication qualification exists.

## Executed real controls

[Portable evidence](protein-completion-branched-evidence.json), raw records under
`tmp/pdbfixer-research/branched-controls-01/`. Python 3.12.13, RDKit 2025.09.6.

Eight masks were declared before execution on the explicit unmodified 1UBQ ATOM
A1–72 projection. Withheld coordinates enter **only later reporting**, never
generation, candidate selection or the independent stereo oracle for proposals.
The existing fixture control harness is reused with an explicit case-list argument;
its original eight-mask default is unchanged.

| Residue | Missing paired branches: selected exact-name RMSD | Also missing CB |
|---|---:|---:|
| VAL5 | 0.027 Å | 0.117 Å |
| VAL26 | 0.245 Å | 0.215 Å |
| THR7 | 0.472 Å | 0.366 Å |
| THR55 | 0.198 Å | 0.176 Å |

All eight have a selected scoped geometry pass. **12/24 proposals** passed
geometry; all 24 proposals had independently correct physical stereo/named-volume
orientation, including proposals rejected for clashes. Every attempt/rejection is
retained. Original typed observations are immutable; no new PDB writer or
accepted branched structure was produced.

These small-mask results are encouraging **case controls**, not template
superiority, empirical rotamer qualification or broad native-quality recovery.
Exact-name RMSD is not symmetry-minimized for VAL. No PDBFixer/tleap comparator
or held-out protein cohort was executed in this step.

The original ALA/SER eight-mask default was rerun in `own-kernel-controls-05.json`;
historical records/source hashes were not replaced. The six previous private
evidence bundles remain subject to unchanged replay semantics.

Tests: **374 passed in 13.53 s** in the scientific interpreter; twelve-file
package mypy and targeted Ruff passed. Tests include observed-atom preservation,
both branch clashes, no-op wrong-stereo rejection, partial-mask/closure refusal,
wrong empirical authority and explicit source/cohort capability isolation.

```bash
# RDKit oracle checks need the scientific environment.
PYTHONPATH=src python -m pytest tests/test_protein_completion_branched.py tests/test_protein_branched_stereo.py -q
PYTHONPATH=src python scripts/probe_branched_completion.py --output tmp/pdbfixer-research/new-branched-controls --summary tmp/pdbfixer-research/new-branched-summary.json
```

Next gates: broader independent structural cohorts, an explicit source-level THR
Cβ policy, paired-branch rounded readback, and only then deliberate downstream
capability expansion. Single-branch anchored completion, ILE/multi-chi chemistry,
real empirical data and production consumers remain open.
