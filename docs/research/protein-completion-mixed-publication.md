# Mixed-policy evidence: separate schema and complete deterministic replay

Status: **private experimental evidence**, not native-quality or production
qualification. Follows the [mixed joint source contract](protein-completion-mixed.md).
The original cohort, public diffusion and experimental diffusion atom-completion
formats remain unchanged and reject this new kind.

Scientific continuation: a [matched template comparison](protein-completion-template-comparison.md)
now evaluates the earlier single-target controls independently of storage/replay
qualification. It does not change this format or promote its scientific flags.

## Format and source authority

`protein_completion/mixed_evidence.py` supplies `write_mixed_evidence` and
`read_mixed_evidence`. Schema 1 uses kind
`experimental-protein-completion-mixed-evidence` and binds the exact mixed policy.
`public_profile_enabled` and `scientific_qualification` must both remain false.

| Stage | Retained artifacts besides `bundle.json` | Meaning |
|---|---|---|
| Admission refusal | `source.pdb`, `targets.json`, `audit.json` | Bounded original bytes/target identities and replayed refusal; no invented successful plan |
| Raw rejection | Above plus `plan.json`, `attempts.json` | All individual attempts and enumerated joint audits, including failures |
| Rounded rejection | Same as raw rejection | Selected raw cohort exists, but rounded physical output fails; no accepted PDB/.dat |
| Accepted rounded output | Above plus `candidate.pdb`, `candidate.dat` | Independent full-source/geometry/stereo/ownership readback passed |

An empty/oversized source or invalid transport-level target array fails before
snapshotting. Within bounded transport, duplicate/absent targets, partial occupancy
and unsupported chemistry can be retained as admission refusals. Planning is a
pure admission function: its `ValueError` refusals are replayed; I/O/runtime errors
propagate. Only the known rounded mixed-geometry refusal from materialization is
classified as rounded rejection; unexpected failures are not relabelled.

Modes and fixed partners are **rederived from source bytes**, never supplied as
scientific authority. The writer optionally accepts an existing mixed plan/result
pair, but canonical replay must reproduce both completely before staging. Missing
attempts, alternative selection, edited observations/partners, wrong mode or
changed global THR-Cβ evidence refuse. Canonical JSON distinguishes booleans from
integers even where Python equality does not.

## Reading is more than hashing

The reader requires the committed index, exact stage-specific flat inventory,
bounded byte counts, SHA-256 and single-linked regular artifacts. It then repeats
source admission, mode/partner derivation, all generation/selection and joint
materialization. Every stored machine-owned payload/index must equal canonical
replay, including reasons for refused admission and full rounded stereo audits.

Accepted PDB bytes are independently re-audited, including new–new clashes,
unchanged observations, exact ownership, bounded rounding, global Cα and THR-Cβ.
The canonical `.dat` is decoded through `DatRecord` and contains only new atoms,
never fixed partners. Editing a payload and updating its index hash is insufficient.

Hashes are **integrity, not authentication**. A fully replaced, internally valid
record is not a signed provenance claim. Consumers can additionally pin the expected
source with `expected_source_sha256`. No hostile-host or multi-tenant security
qualification is claimed. The source here is the supplied masked projection;
original deposited files and mask/projection provenance remain in the retained
research record, not an invented guarantee of full-complex admission.

## Commit and failure lifecycle

Preflight precedes staging/reservation. The writer reuses bounded private-file
reading, durable exclusive writes and the existing manifest-last commit helper.
It independently reads the staged snapshot before committing its readiness index
last. Existing/crash-left/competing destinations are never adopted or overwritten.
Caller inputs remain untouched. Failed staging retains written diagnostics; generic
commit rollback can remove only its own transferred entries. No automatic resume
or cleanup of user-owned reservations exists.

Limits: 16 MiB source, 1–4 transported targets, 64 MiB per artifact, 128 MiB total,
64 kB index. Missing indices, extra files, traversal references, malformed sizes,
symlinks and hardlinks refuse. Old and new readers reject each other's kinds,
including when all chosen residue types could otherwise fit both capabilities.

## Retained-control execution

[Portable record](protein-completion-mixed-publication-evidence.json); committed
private records under `tmp/pdbfixer-research/mixed-publication-01/`.

All eight previous `mixed-cohort-02` cases were checksum-verified, replayed and
committed/read back: **five accepted**, **one admission refusal**, **one joint
rejection**, and **one rounded refusal**. Saved plans/attempts and accepted PDB/
canonical sidecars match retained controls. Refusals contain no accepted structures.
These are unchanged deterministic controls, **not new masks or scientific seeds**.
The harder 1A3N VAL1 1.864 Å whole-sidechain result remains visible in the original
research evidence; storage/replay does not improve its native fit.

Verification: **544 focused tests passed in 22.08 s**, including generic commit
tests and both old/new evidence suites. Seventeen-file package mypy and targeted
Ruff passed. New tests cover mode/partner/stereo/attempt tampering after rehashing,
fixed-partner ownership, admission reasons, source pins, bool/int identity/audit
forgery, individual/joint/rounded failures, competing/crash-left paths, staging
I/O errors, budget limits and cross-reader incompatibility.

```bash
# Cross-reader tests require the scientific environment.
PYTHONPATH=src python -m pytest tests/test_protein_mixed_evidence.py tests/test_diffusion_bundle_commit.py -q
PYTHONPATH=src python scripts/probe_mixed_evidence.py --retained tmp/pdbfixer-research/mixed-cohort-02 --output tmp/pdbfixer-research/new-mixed-evidence --summary tmp/pdbfixer-research/new-mixed-evidence-summary.json
```

No public API/CLI/GUI/backend was registered. Next gates are broader scientific
comparisons and chemistry/mask coverage, not relabelling these private records as
qualified. Licensed empirical data, ILE/multi-chi support and all nine production
PDBFixer rebuilding adapters remain open; existing chirality guards stay intact.
