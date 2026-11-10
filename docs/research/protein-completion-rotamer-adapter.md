# Single-chi rotamer authority adapter: implemented, real-data gate open

Status: private executable adapter; **no authentic Dunbrack statistics have been
loaded or bundled**. This follows the [first completion slice](protein-completion-first-slice.md)
and [replacement requirements](pdbfixer-usage-and-replacement-scope.md). Production
backends, dependencies and scientific gates remain unchanged.

## What is now implemented

- `protein_completion/rotamers.py`: bounded, checksum-checked **offline** reader
  for the traditional single-chi table format; immutable parsed views; explicit
  authority metadata; periodic phi/psi queries, probabilities, circular means and
  variances. Raw data, parsed view, metadata and lookup policy are all digest-bound.
- `protein_completion/backbone.py`: phi/psi from **observed** N/CA/C and explicit
  previous-C / next-N peptide links. There is no sequence-number adjacency guess,
  default terminal angle, source-gap bridge or hypothetical neighbour placement.
- `complete_sidechain(request, library=...)`: use interpolated mean chi values
  instead of the three-angle control; retain all attempts; reject bad geometry
  first, then select by empirical prior with capped clearance as tie break.
- `audit_candidate(..., library=...)`: independently recompute the backbone
  angles and library choices, rejecting wrong authority, altered probability,
  sigma/mean/phi, wrong category, false chi and unlabelled control candidates.

No supplied authority means the **original explicitly labelled staggered control**.
Supplying an authority but lacking data/support produces `unsupported`, **not a
fallback to the control**. The optional empirical path uses result revision
`canonical-sidechain-rotamer-adapter-v2`; frozen public diffusion schemas are not
modified. Candidate records now have optional immutable rotamer evidence, and
results expose the separate model digest; neither field is written to a public
structure/bundle by this prototype.

The underlying coordinate capability remains ALA/SER/free CYS. Single-chi
statistics are admitted for **SER** or explicitly free-cysteine **CYH** only.
Combined/disulfide CYS statistics are not silently labelled free CYS. Other
residues, nonrotameric chi densities, multi-chi construction, linked sidechains,
backbone defects, cap/OXT creation and observed-distal closure remain unsupported.
ALA has no sidechain chi and needs no rotamer library.

## File and interpolation contract

The parser reads plain ASCII local files up to 16 MiB. It requires the documented
traditional 17 fields:

```text
type phi psi count r1 r2 r3 r4 probability chi1 chi2 chi3 chi4 sigma1 sigma2 sigma3 sigma4
```

Only three single-chi labels with zero unused chi fields are admitted. It checks
grid spacing/range, finite parameters, nonnegative normalized cell probabilities,
positive sigma, duplicate rows and conflicting periodic boundary aliases.
External-local data require a full **36×36** periodic ten-degree grid; sparse
synthetic controls stay labelled synthetic. Missing required interpolation cells
refuse rather than substituting a nearby bin. Headers are not executed; there
is no pickle, archive extraction, dynamic code, hidden download or auto-registration.

For phi/psi between four nodes, let `w` be the periodic bilinear weights and `p`
the corresponding probabilities for one rotamer label:

```text
P = sum(w*p)
a = w*p / P
mu = atan2(sum(a*sin(mu_node)), sum(a*cos(mu_node)))
variance = sum(a * (sigma_node² + wrap(mu_node-mu)²))
```

Zero-probability labels are omitted, not assigned a fabricated epsilon. Undefined
antipodal circular means refuse. The returned sigma includes node variance and
between-node variation; it is retained as evidence but **not yet used to sample
additional subrotamers**. This is our versioned circular-mixture interpolation,
not a claim to reproduce the authors' analytic KDE, bicubic derivatives, or all
SCWRL/Rosetta lookup details.

Within the independently feasible candidate set, probability ranks modes first;
clearance resolves ties. This is a prior-driven prototype, **not a calibrated
posterior, H-bond/solvation energy or coupled rotamer packing implementation**.
Even a high-prior candidate is never selected when its geometry audit fails.

## Data acquisition and licensing: what is actually known

The intended first real authority is **Dunbrack 2010, 5% smoothing, traditional
single-chi SER and free-CYS categories**, following the authors' recommendation
of 5% as a useful default. It is not yet a selected byte-identified artifact.

- The [publisher landing page](https://dunbrack.fccc.edu/lab/bbdep2010) as exposed
  by indexed search text states **CC BY 4.0 from 2019-07-25** and requires credit.
- Direct retrieval of the landing page, tutorial and linked terms currently
  times out on this host. Cached/indexed descriptions are useful research
  context, **not exact bytes of a statistical artifact or its license evidence**.
- A publicly indexed Gray Lab copy has a
  [legacy academic license](https://graylab.jhu.edu/Sergey/for.Rocco/nobu-3/ExtendedOpt1-5/AcademicLicenseBBDep2010.txt)
  with noncommercial/nontransfer/redistribution restrictions. The exact terms
  bytes were retained, SHA-256
  `0c4c76ad93c504c3a23e12c6d5364763a5f2314bd8d3ac6cd02f29c62998e61a`.
- This **does not establish that the modern publisher library remains under
  those old restrictions**, nor that all old copies have identical data or
  current redistribution authority. No old mirror statistics were downloaded,
  accepted on the user's behalf, relicensed or substituted as the chosen dataset.
- MEDFORD and unrelated libraries are not Dunbrack substitutes; no rows from
  those repositories were loaded. No 935+ MiB third-party archive was fetched.

Before including or qualifying a real authority, retain exact source/archive
identity, chosen variant/subclass, data SHA-256, applicable license/attribution
evidence, and the relation of any derived single-residue extraction to the
original bytes. Do not fill manifest fields with invented URLs/hashes or treat
a file's presence as permission to redistribute it.

Reader metadata is explicitly **caller-attested provenance**, not a licensing
decision engine. `license_evidence_sha256` identifies supplied evidence; the
numeric parser does not independently authenticate the licensing claim. No
real-data manifest is currently approved. Since no database rows are distributed,
this slice does not add an assumed blanket license notice for an unselected file.
Planned scientific credit when used: Shapovalov & Dunbrack, *Structure* 2011,
19:844–858, DOI [10.1016/j.str.2011.03.019](https://doi.org/10.1016/j.str.2011.03.019).

## Executed verification, clearly separated from real-data qualification

The [portable evidence](protein-completion-rotamer-adapter-evidence.json) contains
the actual authority-page probe and an **explicitly synthetic** end-to-end
transport record:

- Two authored tables, SER/CYH, each 1,296 cells with arbitrary smooth test
  probabilities/means; no experimental statistics or training data.
- Four controls (SER/CYS, terminal-only/whole-sidechain), 12 generated candidates;
  all pass the scoped audit and independent evidence readback. Observations and
  exact masks remain immutable. This is software transport/geometry evidence,
  **not an empirical Dunbrack result or scientific completion benchmark**.
- Boundary wrapping, probability normalization, circular interpolation, zero
  probabilities, conflicting aliases, forged metadata, wrong categories,
  missing neighbours, incomplete cells, checksum failures and geometry-before-
  prior ranking are tested.
- Independent Biopython checks match phi/psi/chi signs on 1UBQ SER20/57/65 and
  generated chi geometry. These use real coordinates to check **conventions**,
  not to claim library sampling or native-conformation accuracy.
- The original eight-mask staggered control was rerun. Generated coordinates,
  audits, statuses and selected indices match the original `controls-03` record
  exactly. New provenance fields do not turn that rerun into a new empirical arm.
- Focused combined suite: **185 passed in 8.18 s**. Package mypy passed for seven
  source files (`--follow-imports=silent`), targeted Ruff passed. Earlier scientific
  controls and their expected values were not rewritten to accommodate this slice.

Raw records are under `tmp/pdbfixer-research/rotamer-authority-01/`,
`rotamer-transport-01/`, `rotamer-transport-02/`, and
`own-kernel-controls-04.json`. The portable record digest-binds the executed
source files; earlier first-slice source hashes are **historical**, not expected
to match the extended implementation. Full-suite, actual artifact compatibility,
empirical ranking accuracy and production acceptance remain unclaimed.

## Reproduction and next gate

These commands use new exclusive destinations inside the worktree:

```bash
PYTHONPATH=src python -m pytest tests/test_protein_rotamers.py tests/test_protein_completion.py tests/test_pdbfixer_usage_inventory.py -q
python scripts/probe_rotamer_authority.py --output tmp/pdbfixer-research/new-authority-probe
PYTHONPATH=src python scripts/probe_rotamer_transport.py --output tmp/pdbfixer-research/new-synthetic-transport --authority-probe tmp/pdbfixer-research/new-authority-probe/evidence.json --summary tmp/pdbfixer-research/new-synthetic-summary.json
```

The authority probe retrieves only pages/terms; it does not accept a license or
download statistics. Synthetic tests and generation work without network access.
The scientific-interpreter test suite also includes
`tests/test_protein_rotamer_conventions.py` for independent Biopython checks.

Next: obtain a legally usable, exact publisher artifact (or explicitly reviewed
local copy); pin it and validate the reader against its real headers/rows and
independent expected values. Only then run a declared real-library comparison
against the unchanged staggered control and PDBFixer/`tleap` where applicable.
Do not turn the synthetic table into a fallback, and do not promote the private
adapter or remove production chirality guards to bypass this open gate.
