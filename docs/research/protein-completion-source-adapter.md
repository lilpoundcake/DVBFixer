# Private source-bound protein completion adapter

Status: **private executable source/readback slice**, not a production backend.
Follows the [kernel](protein-completion-first-slice.md) and
[optional rotamer adapter](protein-completion-rotamer-adapter.md).
No PDBFixer consumer, dependency or default has changed.

Continuation: [declared-link and canonical-graph admission](protein-completion-link-adapter.md)
now supports a narrowly audited subset of LINK/SSBOND/internal CONECT. The scope
and measurements below describe the original source-adapter checkpoint; its
evidence and source hashes remain historical rather than being rewritten.

## Implemented contract

`src/dvbfixer/protein_completion/source.py` provides:

- `plan_source(bytes, exact_residue)`: parse every retained atom using the existing
  PDB component parser. Derive observations, missing atoms and observed peptide
  links from the bytes, rather than trusting caller-supplied context. Full source
  SHA, chain case, insertion codes and link identities bind the kernel request.
- `materialize_source(bytes, request, candidate, library=...)`: rederive and compare
  the complete request; independently audit the raw candidate; insert only new
  `ATOM` records using the existing fixed-column formatter. Never move, reflect,
  renumber, normalize occupancy or rewrite any original record.
- `audit_materialized_source(...)`: independently parse output, remove only exact
  generated identities and require the remaining physical records to equal the
  source byte-for-byte. Check generated ownership, record types, serials,
  occupancy/B-factor, per-coordinate rounding bounds, geometry and protein Cα.
- `MaterializedCompletion.dat_record()`: create a fresh canonical `DatRecord`
  containing exactly generated atoms, with case-sensitive chain/residue/icode.
  No observed atom receives generated ownership. No upstream sidecar is merged.

Both pre-generation and post-readback gates reject non-L or nearplanar protein
N–CA–C–CB centres (signed volume <=0.5 Å³). An absent target CB is allowed only
before generation; absent centres elsewhere refuse. GLY/water have no such centre.
This does **not** check THR/ILE secondary centres or all chemical stereochemistry.

The adapter returns bytes in memory and writes no files. It has no publication
bundle/CLI registration and no implicit retry, repair or production routing. The
separate research script retains inputs, attempts and outputs explicitly.

## Earlier slice: conservative admitted scope and refusals

- One source frame, one target, ALA/SER/free-CYS sidechain masks only.
- Canonical heavy-atom protein names and water `HOH/O`, with explicit elements,
  occupancy exactly 1 and no altlocs. Measured H, PTMs, AMBER/CHARMM variant names,
  caps, ligands, nucleotides and ions refuse in this first adapter.
- Unique decimal atom/TER serials and unambiguous contiguous residue identities;
  canonical decimal residue numbering within PDB limits. A bare `TER\n` is valid.
- LF-terminated ASCII records; no MODEL/ENDMDL, malformed structural prefixes,
  repeated identities, conflicting residue names or nonempty records after END.
- ANISOU identity and integer fields must match an existing atom; duplicates and
  orphans refuse. Valid records remain unchanged. Headers remain unchanged.
- LINK and SSBOND refuse, rather than being ignored. MASTER refuses because
  retaining its old coordinate counts after insertion would be misleading.
- Every CONECT endpoint must exist. External peptide CONECT must agree with
  file-order, same-chain, TER-aware observed C–N geometry (1.0–1.8 Å). Nonlocal,
  nonpeptide, dangling/self links refuse. Target intra-residue CONECT refuses
  until an independent chemical graph audit exists; non-target internal records
  are retained, not advertised as chemically validated.
- A CYS target with another CYS Cα within 10 Å refuses even without SSBOND/SG,
  conservatively bounding possible disulfide reach. This is admission protection,
  not a general disulfide detector.

Other residues' complete sidechains are not certified; all retained atoms still
enter clash context. Existing context-context clashes are not a full-system
geometry qualification. The shared parser imports the scientific/nonprotein
stack; this new adapter is separate from the dependency-light kernel entry point.
No CIF reader or alternate coordinate parser was introduced.

## PDB rounding is a separate scientific gate

Raw proposals retain the original 1e-5-degree χ metadata gate. Rounded readback
has a separately explicit 0.15-degree χ tolerance, reflecting 0.001 Å fixed-column
precision, **only after** the raw candidate passes and each readback coordinate
differs by <=0.000500001 Å. All bond/angle/chirality/clash thresholds remain
unchanged. No production geometry threshold was changed.

A regression deliberately places a water at the clash boundary: raw CB–O is
2.620063 Å, rounded CB–O is 2.619994 Å. The raw proposal passes, but the written
PDB would fail the unchanged -0.6 Å clearance threshold, so materialization raises
and returns no accepted structure. Rounding is never an excuse to waive a clash.

## Executed evidence

[Portable record](protein-completion-source-adapter-evidence.json), raw artifacts
under `tmp/pdbfixer-research/source-adapter-01/`:

- Eight predeclared masks: ALA28/46 CB and SER20/57/65 OG or CB+OG.
- All eight completed, passed raw and rounded geometry, global protein Cα and
  independent file/sidecar readback. Source records are unchanged; sidecars name
  exactly the 11 newly generated atoms across the eight independent requests.
- The source is an **explicit 1UBQ ATOM A1–72 projection**, plus bare TER/END.
  Original/projection bytes and their hashes, pre-execution plan, per-case inputs,
  full requests, all candidate audits, PDBs and canonical `.dat` are retained.
  The full deposited input is separately refused (MASTER and partial occupancy).
  No occupancies were rewritten and no successful full-input admission is claimed.
- Sampling is the unchanged staggered control, **not authentic Dunbrack evidence**.
  These are source/materialization controls, not eight new native-quality wins.
- Combined narrow/regression suite: **231 passed in 9.08 s**; package mypy passed
  for eight files. Broader backend/full-suite acceptance is not claimed.

The earlier rotamer/control evidence keeps its historical source hashes; changes
to `validation.py` for explicit rounded readback do not rewrite those records.

Reproduce in a new exclusive destination:

```bash
PYTHONPATH=src python -m pytest tests/test_protein_completion_source.py -q
PYTHONPATH=src python scripts/probe_protein_completion_source.py --output tmp/pdbfixer-research/new-source-controls --summary tmp/pdbfixer-research/new-source-summary.json
```

## Remaining gates

Authentic licensed rotamer data and scientific comparison remain open. Next
extensions need independently validated LINK/SSBOND/target graph support,
multiple-target atomic composition and evidence/publication contracts before
production adaptation. This slice does not solve multi-chi construction, coupled
packing, backbone reconstruction, missing OXT/caps, full chemistry or replacement
of all nine PDBFixer rebuilding consumers.
