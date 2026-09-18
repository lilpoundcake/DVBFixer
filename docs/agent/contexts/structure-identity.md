# Structure Identity

Status: partial

Verified on: 2026-09-18

Verified at commit: `425f290eb85760246766f1d1500e51672c640b2d`

## Purpose

This context defines the small, format-independent identity vocabulary in
`dvbfixer.domain.structure_identity` and the policy for allocating portable PDB
chain IDs without reusing reserved identities.

It does not yet provide one identity representation across PDB parsing, `.dat`
sidecars, topology construction, CIF normalization, or toolkit objects.

## Scope

The implemented domain module is
[`structure_identity.py`](../../../src/dvbfixer/domain/structure_identity.py).
It contains immutable `ResidueRef`, `AtomRef`, and `MolecularComponent` values,
the `ComponentKind` vocabulary, and the pure `allocate_chain_ids` function.

Only chain allocation is integrated into production behavior. Parsing,
component discovery, record rewriting, CIF conversion, and topology assembly
remain adapter-local. The broader domain migration is partial, as stated in the
[`domain model`](../../domain-model.md) and [`ARCHITECTURE.md`](../../../ARCHITECTURE.md).

## Capabilities

| Capability | Status | Owner | Evidence |
|---|---|---|---|
| Allocate IDs in caller-supplied alphabet order | implemented | `structure_identity.py::allocate_chain_ids` | `tests/test_scientific_domain.py` |
| Exclude nonblank reserved IDs and fail on exhaustion | implemented | `structure_identity.py::allocate_chain_ids` | `tests/test_scientific_domain.py` |
| Reserve polymer and `REMARK 350` IDs for molecule naming | implemented | `molecule_chains.py::assign_unique_molecule_chains` | `tests/test_split.py` |
| Immutable residue, atom, and component vocabulary | partial | `structure_identity.py` | exported, but unused and untested |
| Preserve compatible CIF chain IDs and map incompatible IDs | implemented | `structure_input.py::_chain_mapping` | `tests/test_structure_input.py` |
| One canonical cross-pipeline identity model | proposed | none | not implemented |
| Connected-component identity for legacy PDB heterogens | proposed | none | current boundary is one residue |

## Entry Points

| Task | Start symbol | Status |
|---|---|---|
| Change chain allocation policy | `structure_identity.py::allocate_chain_ids` | implemented |
| Change unique molecule chain assignment | `molecule_chains.py::assign_unique_molecule_chains` | implemented |
| Change CIF-to-PDB chain mapping | `structure_input.py::_chain_mapping` | implemented |
| Change `.dat` atom identity | `ffutils/dat.py::DatRecord.added_keys` | implemented |
| Change topology parser identity | `top/pipeline.py::_read_pdb` and `top/types.py::PDBResidue` | implemented |

`src/dvbfixer/domain/__init__.py` publicly exports `ResidueRef`, `AtomRef`,
`MolecularComponent`, `ComponentKind`, and `allocate_chain_ids`. Repository-wide
usage shows only `allocate_chain_ids` has a caller; the four type exports are
currently unused domain vocabulary, not an integrated abstraction.

## Contracts

The allocator contract is pure: given a nonnegative `count`, an ordered unique
`alphabet`, and `reserved`, return the first `count` available values as a
tuple. Empty and single-space reservations are ignored. Exhaustion raises
`ValueError`. The production caller supplies the 62 distinct PDB IDs; the
function itself does not validate count or alphabet uniqueness.

Current identity shapes are incompatible and must be converted deliberately:

| Location | Shape | Blank insertion code | Notes |
|---|---|---|---|
| Domain residue | `(chain_id, sequence_number: int, insertion_code, name)` | `" "` | frozen and ordered; name participates in equality |
| Domain atom | `(serial: int, residue: ResidueRef, name)` | inherited | serial is part of equality |
| Molecule-chain adapter | `(chain, resseq: raw str, icode, resname)` | `" "` | resSeq retains four-column padding |
| `.dat` added atom | `(chain, resid: str, icode, atom)` | `""` | resname and serial excluded |
| `.dat` variant | `"chain:resid:icode"` | empty third field | documented shape; loader stores strings without validation |
| Topology residue | fields `chain_id, resname, resseq: int, icode` | `""` | mutable; atoms are name/coordinate tuples |
| Topology molecule count | `(chain, resseq: int, resname)` | omitted | cannot distinguish insertion-code siblings |
| CIF chain map | `{original_chain: pdb_chain}` | not applicable | chain identity only; no residue key |

[`../contracts.toml`](../contracts.toml) records the current
`cif-to-pdb-normalization` and `dat-record-lifecycle` boundaries. No single
contract currently owns every identity conversion in the table above.

## Invariants

- Implemented: allocation preserves alphabet order and excludes every nonblank
  reserved value.
- Implemented for the production alphabet: exhaustion fails rather than
  reusing a reserved chain ID.
- Implemented: chain IDs are case-sensitive; `D` and `d` remain distinct.
- Implemented adapter rule: valid one-character CIF chain IDs are preserved;
  only incompatible IDs are remapped.
- Implemented Split rule: atom serials and CONECT endpoints remain unchanged
  while molecule chain fields are rewritten.
- Implemented: `.dat` added-atom identity is `(chain, resid, icode, atom)`.
- Partial: variant identity is documented as `(chain, resid, icode)`, but legacy
  two-field strings remain accepted and some local maps omit insertion code.
- Missing: a normalized blank-insertion-code and residue-number representation
  shared by every adapter.
- Missing: a MODEL discriminator in the domain identity values.

The binding repository rules are in [`AGENTS.md`](../../../AGENTS.md), notably
CIF normalization at one boundary, case-sensitive chain IDs, insertion-code
preservation, `.dat` access through `DatRecord`, and Split serial preservation.

## Callers

- Direct production caller: `molecule_chains.py` imports only
  `allocate_chain_ids`.
- Indirect caller: `split_chains.py` invokes `assign_unique_molecule_chains`
  for `--unique-molecule-chains` paths.
- `structure_input.py` implements its own preservation-first chain mapping and
  does not call the domain allocator.
- `ffutils/dat.py` and `top/types.py` define independent identity containers and
  do not import the domain types.
- No production or test caller constructs `ResidueRef`, `AtomRef`,
  `MolecularComponent`, or `ComponentKind`.

## Adapters

- `molecule_chains.py` parses fixed-column PDB records, discovers heterogen
  groups, gathers coordinate and assembly-reserved IDs, calls the allocator,
  and rewrites identity-bearing records.
- `structure_input.py` is the sole CIF boundary. Gemmi supplies author chains
  and `label_asym_id`; conversion emits `CIF_CHAIN_MAP` and positional
  `CIF_COMPONENT` remarks for downstream PDB processing.
- `ffutils/dat.py` is the pipeline sidecar adapter. It serializes string residue
  numbers and empty insertion codes and merges downstream values over upstream
  collisions.
- `top/types.py` is a mutable topology-building representation, not an alias or
  adapter for the immutable domain references.

## Side Effects

- `allocate_chain_ids` has no I/O and does not mutate its arguments.
- `assign_unique_molecule_chains` returns rewritten lines, inserts a
  `UNIQUE_MOLECULE_CHAINS` remark, and changes chain fields in applicable
  HETATM, ANISOU, HET, and LINK records; it leaves CONECT text and serials
  unchanged.
- CIF normalization writes a validated temporary PDB and identity remarks,
  translates supported selectors and FASTA headers, and propagates mapping
  remarks to generated PDB outputs.
- `DatRecord.save` writes JSON directly; topology models accumulate mutable
  residue and atom lists.

## Known Divergences

- The exported domain dataclasses are not used by current adapters or pipelines.
- Blank insertion codes are `" "` in the domain/PDB adapter but `""` in `.dat`
  and topology parsing.
- Residue numbers are integers in the domain/topology model, raw padded strings
  in molecule-chain grouping, and unpadded strings in `.dat`.
- Domain residue equality includes residue name; `.dat` added-atom identity does
  not, while molecule-chain grouping and topology counting include it.
- Domain atom equality includes serial; `.dat` identity excludes serial, and
  topology residue atom tuples contain no serial.
- `_count_molecules` in `top/pipeline.py` omits insertion code from its key.
- `DatRecord` does not validate `variant_overrides`; tests still round-trip
  legacy `"A:39"` keys despite the documented `"A:39:"` contract.
- Legacy PDB heterogens are grouped per residue, not by CONECT connected
  component; mmCIF grouping relies on positional `CIF_COMPONENT` remarks.
- The allocator error text assumes the standard 62-ID alphabet even though the
  function accepts any caller-provided alphabet.
- The allocator does not reject duplicate alphabet entries or negative counts;
  its uniqueness contract therefore depends on caller validation.
- None of the identity values carries MODEL identity, alternate location, or
  source-file identity.

## Proposed Work

- Define explicit conversion functions before routing adapter values through
  `ResidueRef`; do not rely on tuple-position compatibility.
- Decide whether residue name and atom serial are identity or metadata for each
  operation, and whether MODEL must be represented.
- Normalize insertion-code blanks only at documented boundaries while
  preserving fixed-column PDB rendering.
- Either integrate and test the four exported domain types or remove them in a
  coordinated API change; do not claim they are active today.
- Validate allocator inputs if it becomes a general public policy rather than a
  helper used with the fixed, unique PDB alphabet.
- Consider sharing allocation policy with CIF mapping without changing CIF's
  preservation-first behavior.
- Add connected-component grouping only with serial/CONECT and mmCIF
  `label_asym_id` regression coverage.

## Focused Verification

```bash
pytest -q tests/test_scientific_domain.py tests/test_split.py \
  tests/test_structure_input.py tests/test_ffutils_dat.py \
  tests/test_top_merge_chains.py
python scripts/check_agent_docs.py
git diff --check -- docs/agent/contexts/structure-identity.md
```
