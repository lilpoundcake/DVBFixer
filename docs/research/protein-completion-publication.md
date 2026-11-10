# Private cohort evidence: manifest-last commit and independent replay

Status: executable **private evidence format**, not production qualification or a
public preparation/diffusion backend. Follows the
[bounded cohort controls](protein-completion-cohort.md).

Continuation: [typed VAL/THR paired-branch controls](protein-completion-branched.md)
are now implemented, but are **not enrolled** in this source/cohort/publication
capability. Existing ALA/SER/free-CYS replay semantics remain unchanged.

## Format and admission

`protein_completion/evidence.py` implements `write_cohort_evidence` and
`read_cohort_evidence`. Its index has schema 1 and kind
`experimental-protein-completion-cohort-evidence`; both `public_profile_enabled`
and `scientific_qualification` must be **false**. Existing public diffusion and
experimental diffusion atom-completion readers reject this kind.

Every committed record contains:

| Artifact | Meaning |
|---|---|
| `source.pdb` | Exact bytes supplied to the cohort, never rewritten |
| `targets.json` | Exact, ordered target identities |
| `plan.json` | Complete immutable source-bound requests/context/links/masks |
| `attempts.json` | All individual and enumerated joint attempts, including failures |
| `audit.json` | Raw/rounded stage, refusal or accepted generated atoms, private-only flags |
| `candidate.pdb`, `candidate.dat` | Present **only** for accepted rounded readback |
| `bundle.json` | Index: exact flat inventory, byte counts, SHA-256, source/plan/policy binding |

`accepted` means this private geometry/source/rounding contract passed, **not**
native quality, MD readiness, empirical rotamer qualification or public readiness.
Raw joint rejection and known rounded-geometry rejection retain source, plan and
all attempts with status `rejected`, and contain no accepted PDB or `.dat`.
Unsupported pre-admission inputs, encoding failures, programming errors and I/O
errors propagate; they are not relabelled as scientific geometry rejections.

The format binds the supplied source, not undocumented preprocessing. In these
controls that source is an explicitly labelled projection; deposited originals
and projection/mask provenance remain in the separate retained research record.
No frozen request or input occupancy is silently changed.

## Why hashes alone are insufficient

The reader first checks fixed filenames, sizes, single-linked regular files,
digests and complete inventory. It then rederives the plan from the source,
replays **all** deterministic attempts/selection, rematerializes and independently
audits the stored output. Every machine-owned JSON payload must equal its
canonical replay, including retained rejection findings. Unexpected fields,
omitted attempts, altered status or selected index cannot be hidden by updating
artifact hashes.

`.dat` decoding stays in `DatRecord`; exact canonical bytes additionally protect
derived totals/fields that its compatibility decoder may normalize. Observed
atoms cannot acquire generated ownership. The reader also runs the independent
combined PDB readback audit on stored candidate bytes.

Canonical-byte comparison matters: Python considers `True == 1`. The writer
rejects a forged numeric audit flag even when dataclass equality would pass.
Exact replay intentionally fails closed if a future implementation changes old
semantics. New behavior needs an explicit policy/version decision, not rewriting
historical records to force acceptance.

SHA-256 is **integrity, not authentication**. A fully substituted, internally
valid record is not a signed claim of authorship. Consumers can bind an external
expected source digest through `expected_source_sha256`. No remote hostile-host,
multi-tenant or general publication security qualification is claimed.

## Filesystem behavior

Preflight replay and artifact limits occur before staging/reservation. Payloads
are bounded at 64 MiB each and 128 MiB total; the index is bounded at 64 kB.
The implementation reuses the existing private-artifact reader and
`model/diffusion/bundle_commit.py::commit_manifest_last`:

1. Create an owned private staging directory; exclusively write/fsync payloads.
2. Write the staged index and independently read/replay the entire staged record.
3. Exclusively reserve the destination, transfer durable payloads, then commit
   `bundle.json` last and sync the parent.

Destination directory existence alone is never readiness. Existing destinations,
crash-left reservations and competing reservations are not adopted or overwritten.
Caller files are never edited/deleted. Staging I/O failures retain already-written
diagnostics; commit rollback can remove only its own transferred entries, as
specified by the shared helper. Failed staging may be incomplete and is not
automatically resumed. Readback rejects symlinks, hardlinks, traversal filenames,
extra/missing artifacts and missing indices.

## Executed retained-control publication

[Portable record](protein-completion-publication-evidence.json), committed private
records under `tmp/pdbfixer-research/cohort-evidence-01/`:

- Reused all six `source-cohort-01` inputs and checksum-verified saved requests,
  attempts and accepted outputs. No new masks or stochastic/scientific seeds.
- Four accepted: three real 1UBQ projection cohorts and the synthetic two-SER
  feasible-combination control. Stored PDB and canonical sidecar content match
  the earlier retained outputs exactly.
- Two rejected: synthetic no-joint-solution and rounded-boundary refusal. Both
  preserve full attempts without accepted candidate files.
- Every private record passed readback with an external source SHA pin. Prior
  records and their source hashes remain historical and unchanged.

Combined focused tests: **341 passed in 13.82 s**, including shared manifest-last
commit tests. Eleven-file package mypy and targeted Ruff passed. Tampering tests
cover rehashed semantic edits, source pinning, scientific/public promotion,
observed-atom `.dat` ownership, missing attempts, malformed references, symlinks,
hardlinks, staging failures, competing reservations and no-overwrite behavior.

The public-reader compatibility test uses the scientific environment: the system
Python lacks OpenMM. It passed in the project scientific interpreter; no broad
skip or environment repair was used. No full-suite/backend acceptance is claimed.

```bash
# Run tests in the project's scientific environment.
PYTHONPATH=src python -m pytest tests/test_protein_completion_evidence.py tests/test_diffusion_bundle_commit.py -q
PYTHONPATH=src python scripts/probe_protein_completion_evidence.py --retained tmp/pdbfixer-research/source-cohort-01 --output tmp/pdbfixer-research/new-cohort-evidence --summary tmp/pdbfixer-research/new-cohort-evidence-summary.json
```

The helper refuses existing output destinations. Replaying unchanged deterministic
controls is evidence validation, not independent scientific replication. Actual
licensed rotamer data, broader chemistry, native-quality comparison and production
consumer adapters remain open; this format does not close those gates.
