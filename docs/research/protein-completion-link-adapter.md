# Private protein completion: declared links and canonical target graphs

Status: **private executable continuation** of the
[source adapter](protein-completion-source-adapter.md). Production defaults,
consumers, geometry thresholds and PDBFixer dependencies remain unchanged.

Continuation: [bounded multi-target staggered controls](protein-completion-cohort.md)
now share source insertion/readback helpers. The link-control evidence below
keeps its historical source hashes; no empirical or public capability is implied.

## Newly admitted metadata, not new coordinate capabilities

- `protein_completion/links.py::parse_declared_links` resolves LINK/SSBOND atom
  endpoints against the complete observed source. Chain case, residue number,
  insertion code, residue name and atom name must match exactly. Missing atoms,
  altlocs, ambiguous endpoints and malformed/truncated fields refuse.
- Blank or `1555` symmetry is admitted in the deposited frame. Other operators
  refuse; no crystallographic mate is silently substituted. Optional reported
  distances must be finite, positive, and agree with observed distance within
  0.1 Å. Absence of a reported distance never bypasses peptide/disulfide geometry.
- `validation.py::sidechain_bonds` now shares the existing ALA/SER/CYS target
  graph with source admission and candidate auditing. Canonical internal CONECT
  or LINK can be preserved without creating external-link exclusions. Unexpected
  edges refuse; metadata cannot overwrite or remove implicit canonical edges.
  Optional OXT connectivity is admitted, not OXT generation.
- Internal declared graphs for other residue types **refuse** until their graph
  authority exists. This tightens the earlier handling of non-target internal
  CONECT, which was retained without chemical validation.
- Explicit peptide LINK/CONECT must match an independently observed neighbouring
  C–N pair: same case-sensitive chain, file-order adjacency, no intervening TER,
  1.0–1.8 Å. A LINK never invents a peptide partner or bridges a source gap.
- Context disulfides can be declared by SSBOND, LINK or CONECT: canonical CYS SG
  endpoints, observed distance **1.8–2.3 Å**, and **one partner per SG**. Duplicate
  descriptions of the same pair deduplicate in the request while their original
  records remain byte-identical in the PDB. Conflicting partners refuse.
- A disulfide-linked **target** still refuses: single-residue completion is not
  a coupled disulfide solver. The prior conservative nearby-CYS target guard is
  retained for unannotated chemistry.

All CONECT endpoints must exist. Self-links and repeated target serials within
one record refuse; multiplicity/bond-order semantics are not silently interpreted.
Crosslinks other than admitted peptides and context disulfides remain unsupported.

These checks concern heavy-atom connectivity and declared-link geometry. They
are not a general bond-order, valence, protonation, stereochemistry or full-system
geometry audit. In particular, source context-context strain is not repaired.
Observed coordinates and all retained metadata remain immutable.

## Executed controls

[Portable evidence](protein-completion-link-adapter-evidence.json), raw artifacts:
`tmp/pdbfixer-research/source-links-01/`.

A predeclared **explicit 1DQJ projection** retains original ATOM bytes for ALA A9
(CB withheld), CYS A23 and CYS A88 with explicit bare TER boundaries. Metadata
controls use the deposited first SSBOND and an authored canonical A9 N–CA CONECT:

| Control | Result |
|---|---|
| No declared disulfide, control only | ALA9 CB completed and read back |
| Deposited CYS23–CYS88 SSBOND | Completed; SSBOND and observations unchanged |
| SSBOND plus canonical target CONECT | Completed; both records unchanged |
| Nonidentity symmetry in the SSBOND | Refused |
| Disulfide CYS23 selected as target | Refused |

The three accepted controls produce identical generated coordinates. The first
control is **not** a disulfide-detection success: it deliberately omits metadata,
and no general automatic detector is implemented. The full original 1DQJ file
is retained but was not admitted wholesale; no occupancy/altloc edits occurred.
The linked-target refusal happens before the unrelated withheld A9 CB global gate.

All original eight 1UBQ source controls were rerun in
`tmp/pdbfixer-research/source-adapter-02/` with a separate retained summary.
Their earlier evidence remains historical; do not replace old source hashes.

Tests cover independent fixed-column identity resolution, chain `B` versus `b`,
insertion codes, absent endpoints, bad names/residue labels, symmetry, reported
and actual distances, invalid graph edges, CONECT multiplicity, TER barriers,
one-partner conflicts across declarations, linked-target refusal and file/sidecar
readback. The real 1DQJ projection is also an automated regression.

Combined focused suite: **266 passed in 9.84 s**. Nine-file package mypy and
targeted Ruff passed. No full scientific suite or production qualification is
claimed. Sampling remains the unchanged staggered control; no authentic Dunbrack
artifact was obtained or substituted.

Retained source/artifact SHA checks and agent maps passed. The eight 1UBQ output
PDB hashes and generated coordinates match the original source-adapter checkpoint
exactly; the metadata work did not change that control arm.

```bash
PYTHONPATH=src python -m pytest tests/test_protein_completion_links.py tests/test_protein_completion_source.py -q
PYTHONPATH=src python scripts/probe_protein_completion_links.py --output tmp/pdbfixer-research/new-link-controls --summary tmp/pdbfixer-research/new-link-summary.json
```

Next gates remain authentic empirical data/qualification, broader chemical graph
authority, coupled-target composition and a separately audited evidence/publication
contract. No public backend, CLI, registry or multi-target capability is enabled.
