# Explicit VAL/THR source adapter and global THR Cβ gate

Status: **separate private single-target source experiment**, not a production
backend. Follows the [typed paired-branch kernel](protein-completion-branched.md).
Default `plan_source`, cohort planning and replay-verified evidence publication
still accept only ALA/SER/free-CYS targets; they do not automatically gain VAL/THR.

Continuation: a [separate typed observed-partner singleton kernel](protein-completion-anchored.md)
now exists. This paired-source policy deliberately still refuses singleton masks;
its historical controls/hashes and global-stereo contract are not rewritten.

## Separate capability and source-bound plan

`protein_completion/branched_source.py` supplies:

- `plan_branched_source(source, exact_target)`: explicit VAL/THR admission using
  the existing parser, canonical target graph and peptide/link/layout rules.
- `BranchedSourcePlan`: the complete request, input THR-Cβ checks and policy
  `paired-val-thr-source-global-thr-beta-v1`; all bind the plan digest.
- `materialize_branched_source(...)`: independently rebind/audit the raw request
  and candidate, insert only new atoms, then independently audit rounded output.
- `audit_branched_materialized_source(...)`: require exact source records,
  generated ownership, bounded rounding, target geometry/stereo/clashes, global
  protein Cα and the post-output global THR-Cβ check.
- A returned `MaterializedBranchedSource` can create a fresh canonical `DatRecord`
  owning only generated atoms, preserving chain case and insertion codes.

The default source implementation has an internal explicit capability argument
for reuse; its default remains unchanged. The new entry restricts targets to
VAL/THR and admits only the corresponding audited canonical internal graphs,
alongside the pre-existing ALA/SER/CYS graphs. No new CIF/coordinate reader was
added. Linked targets, nonidentity symmetry and unsupported graphs still refuse.

Both paired distal atoms may be missing, with CB retained or missing. A singleton
missing branch or observed-distal closure remains unsupported. Source scope is
still canonical heavy-atom protein/water, unit occupancy, no altloc/MODEL or
measured H, and the existing explicit layout/ANISOU/serial constraints.

## What the new global gate does

Before generation, every retained THR must have CA/CB/OG1/CG2 with positive named
signed volume **>0.5 Å³**, corresponding to natural CB=R in this canonical graph.
Only the exact admitted target can be `target_pending`; missing THR branch atoms
elsewhere are **undefined**, not a successful stereo check. Chain case and
insertion codes do not share the pending-target exemption.

After rounded readback **no exemption remains**. Every THR must be complete at
that centre and pass the same gate. The unchanged global Cα gate independently
rejects D/nearplanar/undefined protein Cα. No observation is reflected or repaired.
Tests include a wrong or nearplanar context THR CB with a still-correct Cα.

VAL CB is not a stereocentre: its target named-branch orientation is checked by
the existing kernel, not relabelled as physical R/S. ILE secondary centres and
other chemical stereocentres are not globally validated. This gate is not a
global bond/angle/energy or MD-readiness assessment.

PDB readback keeps the existing <=0.000500001 Å per-coordinate rounding bound,
separate 0.15° rounded χ tolerance and unchanged bond/angle/clash thresholds.
Synthetic water-boundary controls pass raw geometry but fail after rounding;
both VAL and THR materializers raise without returning an accepted structure.
ANISOU, bare/sparse-serial TER and observed records remain byte-identical; serial
exhaustion also returns no partial accepted result. Core code does not import RDKit.

## Executed three-structure controls

[Portable evidence](protein-completion-branched-source-evidence.json), complete
retained inputs, plans, attempts, outputs, sidecars and audits under
`tmp/pdbfixer-research/branched-source-01/`.

All **18 cases were declared before execution**. Target selection on the added
structures used residue identities, not proposal accuracy or held-out geometry.

| Explicit ATOM chain-A projection | Target residues | Paired masks | Outcome |
|---|---|---:|---|
| 1UBQ, observed residues through 72 | VAL5/26, THR7/55 | 8 | All raw/rounded/source/sidecar checks passed |
| 1FR2, all observed chain-A ATOM records | VAL19, THR11 | 4 | All passed |
| 1A3N, observed chain-A residues through 59 | VAL1, THR8 | 4 | All passed |
| 1A3N, all observed chain-A ATOM records | VAL1, THR8, whole sidechain withheld | 2 | Admission refused: partial occupancy remains present |

Every projection is explicitly labelled, with unmodified retained ATOM bytes,
declared withheld masks and bare TER/END boundaries. Original deposited inputs
and their hashes are retained separately. No occupancy was rewritten, no partial
residue was silently selected, and these are **not intact-complex admissions**.

The 16 accepted files additionally passed **108 RDKit rounded-coordinate checks**
on each file's THR centres and selected VAL/THR target. These repeat fixed centres
across masks; they are not 108 independent scientific samples. All source THR
centres passed without moving observations. Candidate `.dat` files were read back
through `DatRecord` and exactly matched generated ownership.

Held-out exact-name RMSD was computed only after selection/readback. It ranges
**0.027–1.864 Å**: notably 1A3N VAL1 with the whole sidechain withheld is **1.864 Å**,
versus **0.281 Å** with observed CB retained. A geometry pass is therefore **not a
native-conformation guarantee**. No truth-based reranking, new sampling mode,
force-field relaxation or threshold adjustment was used to fix this poorer case.
No empirical Dunbrack data or PDBFixer/tleap scientific comparator was executed.

Verification: **403 focused tests passed in 15.85 s**; package mypy passed for
thirteen files and targeted Ruff passed. Tests independently assign CIP after
rounding all three synthetic χ proposals, retain prior capability isolation,
and check source/stereo/policy tampering and global undefined/wrong THR centres.
The original eight source controls were also rerun under `source-adapter-04/`;
earlier evidence and source hashes are not overwritten.

```bash
# RDKit checks require the project's scientific environment.
PYTHONPATH=src python -m pytest tests/test_protein_branched_source.py tests/test_protein_branched_stereo.py -q
PYTHONPATH=src python scripts/probe_branched_source.py --output tmp/pdbfixer-research/new-branched-source-controls --summary tmp/pdbfixer-research/new-branched-source-summary.json
```

## Remaining gates

This service returns bytes in memory; only the explicit research harness retains
files. It does **not** expose the new residues through default cohort/evidence
formats, CLI, GUI or any production consumer. Further work needs a deliberately
versioned joint/source-stereo/evidence contract, more demanding masks and broader
scientific comparisons. Singleton-anchor completion, ILE/multi-chi geometry,
licensed empirical statistics and full-chemistry coverage remain open.
