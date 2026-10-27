# Protein-complex gaps: geometry before observed-anchor relaxation

Date: 2026-10-27. Branch: `research/9cqd-anchor-relaxation`.
Purpose: detect bad raw and generated-only geometry **before** allowing observed
anchors to move. This is a geometry stress test, not reconstruction-accuracy or
neural-profile qualification.

## Frozen protocol and cohort

Nine experimental X-ray structures, 15 selected internal gaps, seed **7** only.
Missing sequences were checked against deposited REMARK 465 and matched to the
SEQRES-based placement and original author-numbered observed flanks.

- Biological assembly **1**, rendered from REMARK 350 without renumbering.
  1BRS uses A+D, 1I2M uses A+B; 1KTZ includes both BIOMT operators.
- Residue-coherent occupancy selection of altlocs. No coordinate repair of the
  experimental source before sampling.
- **Protein-only** assembly copies: HETATM/water and associated CONECT/LINK
  records removed explicitly. Results do not assess removed glycans, ligands,
  Mg/GNP or their interactions. This especially limits 1N8Z, 1G9M, 1GC1 and 1HE8.
- Authentic parent SEQRES targets, terminal gaps cropped with `no_terminal`.
  Selected gaps use the shipped independent per-gap views. Other missing regions
  are neither predicted conditioning nor silently deleted from parent targets.
  These are isolated selected-gap challenges, not fully reconstructed complexes;
  interactions between independently generated, unmerged loops are not assessed.
- Protpardelle cc89, native Apple MPS, 500 steps, step scale 1.2, s_churn 200,
  neutral chain phase, source-physical positions, no per-step reinjection and no
  CPU sampling fallback. Generated-only refinement uses OpenMM CPU, no solvent,
  production iteration/restart settings and a fixed **240 s** process timeout.
- Observed-anchor relaxation is **disabled throughout**, equivalent to
  `--no-relax-anchors`. No replacement seeds or threshold tuning after outcomes.
- Lengths outside 5–10 residues are extrapolation beyond the established cc89
  length evidence. Old public structures may overlap pretraining; this is not a
  blind independent accuracy holdout. Missing loops have no coordinate reference.

Frozen cc89 checkpoint SHA-256:
`dfc9895b399ec4497bf6d646502168725f01dcbd55bc0b541fc3e3cc1f2f0483`;
config SHA-256:
`e9999ace79bf3044351cc982a624fc3459add3a4f91cbf97ff5b437b8942eb9d`.

| PDB / selected gap | Length | Missing sequence | Observed flanks |
|---|---:|---|---|
| 1BRS D:64–65 | 2 | `EN` | THR63 / GLY66 |
| 1N8Z C:303–305 | 3 | `EDG` | ALA302 / THR306 |
| 1N8Z C:361–364 | 4 | `GDPA` | ASP360 / SER365 |
| 1I2M A:32–36 | 5 | `TGEFE` | LEU31 / LYS37 |
| 1I2M B:233–238 | 6 | `SRGSRG` | LYS232 / HIS239 |
| 1JPS T:84–90 | 7 | `ESTGSAG` | VAL83 / GLU91 |
| 1AHW C:83–90 | 8 | `VESTGSAG` | ASN82 / GLU91 |
| 1N8Z C:102–110 | 9 | `NNTTPVTGA` | LEU101 / SER111 |
| 1N8Z C:581–590 | 10 | `GVKPDLSYMP` | SER580 / ILE591 |
| 1G9M G:398–409 | 12 | `STWSTEGSNNTE` | ASN397 / GLY410 |
| 1GC1 G:397–409 | 13 | `NSTWSTKGSNNTE` | PHE396 / GLY410 |
| 1KTZ A:55–72 | 18 | `DTTHSTVLGLYNTLNPEA` | ALA54 / SER73 |
| 1HE8 A:436–456 | 21 | `GKAPALSSKASAESPSSESKG` | CYS435 / LYS457 |
| 1HE8 A:322–352 | 31 | `EWPLADDCTGVTGYHEQLTIHGKDHESVFTV` | GLU321 / SER353 |
| 1HE8 A:528–619 | 92 | `ALPKHQPTPDPEGDRVRAEMPNQLRKQLEAIIATDPLNPLTAEDKELLWHFRYESLKHPKAYPKLFSSVKWGQQEIVAKTYQLLARREVWDQ` | ILE527 / SER620 |

## Which stage passed?

**13/14 sampled raw candidates fail geometry.** One passes without minimization,
10 pass after generated-only minimization, one retains bad geometry, and two
generated-only processes time out. One source is rejected before sampling.
No observed-anchor minimization was attempted.

J = bad junction C–N; C = severe clash; R = gross Ramachandran outlier;
A = backbone-angle finding; ω = non-planar amide. Counts are detailed findings,
not duplicate aggregate gates. Local scope includes the selected generated
region, flanks/outer neighbors and their contacts; PASS is not whole-complex MD
readiness. Junction C–N must be 1.20–1.45 Å; coarse internal break detection uses
1.80 Å and severe-clash overlap uses 0.90 Å. Other shipped geometry gates remain.

| PDB / gap | Passed before any minimization? | Passed after generated-only? | Passed after anchor relaxation? | Raw findings |
|---|---|---|---|---|
| 1BRS D:64–65 | Input rejected | Not run | Not run | Incomplete GLN D:58 context |
| 1N8Z C:303–305 | No | Yes | Not needed | J1 |
| 1N8Z C:361–364 | No | Yes | Not needed | J2 |
| 1I2M A:32–36 | No | Yes | Not needed | C1 |
| 1I2M B:233–238 | No | Yes | Not needed | C1 |
| 1JPS T:84–90 | No | Yes | Not needed | J2 |
| 1AHW C:83–90 | **Yes** | Not needed | Not needed | None |
| 1N8Z C:102–110 | No | Yes | Not needed | R1, C1 |
| 1N8Z C:581–590 | No | Yes | Not needed | C3, J1 |
| 1G9M G:398–409 | No | Yes | Not needed | C1 |
| 1GC1 G:397–409 | No | **No: G:398 Ramachandran** | **Not run** | R1, C2, J2 |
| 1KTZ A:55–72 | No | Yes | Not needed | A1, ω1, C6 |
| 1HE8 A:436–456 | No | Yes | Not needed | C1 |
| 1HE8 A:322–352 | No | **Unknown: timeout 240 s** | Not run | C18 |
| 1HE8 A:528–619 | No | **Unknown: timeout 240 s** | Not run | C85, J1 |

All observed atom records remain exact in the 14 raw PDBs and 11 saved stage-1
outputs. Detectable D-Cα count is zero in all 25. Timeout is not a geometric
stage-1 failure, and “not run” is not PASS or FAIL. The anchor stage was disabled;
the cohort provides no measurement of repair efficacy with mobile anchors.

## 1BRS: input-policy refusal, not a failed two-residue reconstruction

Observed **GLN D:58 lacks OE1 and NE2**. The automatically chosen single-gap
profile requires complete heavy-atom context and refuses this input before
neural generation. This is not evidence that the gap cannot be reconstructed.

The versioned independent/sparse-context path can represent incomplete observed
sidechains when backbone N/CA/C/O is complete, without inventing conditioning
coordinates. That alternative was **not executed for 1BRS in this frozen run**.
A separate explicitly declared sparse-context test, or an explicitly repaired
source with separate provenance, is possible follow-up. Do not retroactively
classify the unattempted gap as a model failure or silently complete the source.

## 1GC1: minimization helped, but did not fix the conformation

| Measurement | Raw | Generated-only output |
|---|---:|---:|
| Left junction C–N | 1.094 Å | 1.319 Å |
| Right junction C–N | 1.459 Å | 1.347 Å |
| Severe clash findings | 2 | 0 |
| Generated SER G:398 φ / ψ | 171.0° / −15.2° | 159.7° / −80.1° |
| Gross Ramachandran outlier | Present | **Still present** |

Bond distances and clashes alone would falsely accept this output. Local
minimization does not guarantee escape from an incorrect loop conformation.
This single-seed outcome does not establish that the gap is unrepairable.

The residual finding names **generated G:398**, not observed flanks G:396/G:410.
The shipped `problematic_anchors` returns an **empty set** for the saved stage-1
assessment. It therefore provides no current justification for automatically
releasing observed anchors. A separate resampling or fixed-anchor conformational
search is a reasonable follow-up, not a claim that moving anchors will solve it.

## Diagnostics and the raw/final boundary bug

All 14 neural jobs completed 500 steps and wrote raw PDB/trace files, then exited
with code 1 in post-sampling diagnostics: `validate_runner_result` incorrectly
required evidence of a **future** refinement stage. These are neither 14 failed
neural generations nor 14 geometry failures. Original logs/outcomes are retained.

The correction adds explicit `phase="raw"` diagnostics to the Protpardelle,
Protenix and Boltz checkpoint adapters. Raw mode keeps strict observed-coordinate
checks and rejects post-refinement traces; it cannot grant a source-motion
exception. `build_validated_result` still uses final mode and requires the staged
request's evidence. Regression tests cover good/bad raw geometry, observed drift,
missing final evidence and rejection of refined traces in raw mode.

**No neural resampling:** all raw PDB/trace hashes were verified and the same
coordinates were evaluated with the explicitly archived replay utility and
`--no-relax-anchors`. Missing original runner resource/probe result objects are
not fabricated. This recovery establishes coordinate geometry, not complete
public sampling-profile/hardware qualification. See the tracked
[replay utility](../../scripts/replay_staged_refinement.py).

Two additional analysis safeguards:

- Incomplete source PDBs can acquire artificial peptide bonds directly across an
  annotated missing loop in OpenMM's topology. Related source bond/angle/ω
  findings are labelled expected gap artifacts, not experimental defects.
  No corresponding exception is applied to the filled raw candidate. After this
  classification, these 14 local source neighborhoods have no other findings;
  that is not a claim that the entire experimental sources are flawless.
- After a failed generated-only attempt, the runner may return raw for inspection.
  Audit the **actual saved stage-1 PDB**, not that returned raw. A synthetic
  partial-repair failure self-test confirms this distinction.

## Sources and evidence custody

Original downloads, not generated coordinates:

| Source | SHA-256 |
|---|---|
| [1BRS.pdb](https://files.rcsb.org/download/1BRS.pdb) | 0076aae54252706a8c8c02f7b1e71ce2394b0e928b413590b4ac12aade9e8f6c |
| [1N8Z.pdb](https://files.rcsb.org/download/1N8Z.pdb) | ed934d36503335c6764b08880f9e50933d99eb13b9bf145339cfcaee150eb21e |
| [1I2M.pdb](https://files.rcsb.org/download/1I2M.pdb) | d11d60e277912fc62f2379ae4d2a64fb5d7970dbc9742197a11ea6f4cf4a096c |
| [1JPS.pdb](https://files.rcsb.org/download/1JPS.pdb) | 22cb05ba2059722d2a571337de9b4f0ef6adc8946b13d16e5997d972a9210e14 |
| [1AHW.pdb](https://files.rcsb.org/download/1AHW.pdb) | bc73f266b34041673df285e1e0c66a507d9452642ff64451424167d32e117cdf |
| [1G9M.pdb](https://files.rcsb.org/download/1G9M.pdb) | b1d4b3cfb26ca9c5bb2ac5c83da0a3255e047ba7812e1f0a1ca285a58396541d |
| [1GC1.pdb](https://files.rcsb.org/download/1GC1.pdb) | 73fb7c47f6b4c5ce0fe51888cf5e32bd5ee3de96bcc57e2fbaf10692119bff60 |
| [1KTZ.pdb](https://files.rcsb.org/download/1KTZ.pdb) | 54adccd77517b42c9d6019afe4eee8c0a6aa9e7fd6966a0c1121eda7b3e9426e |
| [1HE8.pdb](https://files.rcsb.org/download/1HE8.pdb) | 1814c7e5e4a5104301564d41d01d82482958bf20b46ef9eab16021feebc4b489 |

Optional local evidence remains under `.artifacts/protein-complex-gaps/`, not
required by CI or distributed as tracked regression fixtures. It includes the
Russian README, source/assembly/altloc provenance, frozen requests, code snapshots,
all original failures, 25 review PDBs and independent reports. The scripts
`benchmark.py`, `recover_raw.py`, `summarize.py` and `make_report.py` live there.
Raw artifacts and the coordinator environment are not available in a fresh clone;
the hashes below identify this execution, not a claim of self-contained replay.

Paths below are relative to that directory's `runs/pre-anchor-v2/`:

| Evidence | SHA-256 |
|---|---|
| `protocol.json` | 22ce48f2514f4d75fb555487938a26195f56d228efd4d5ae70fa5e3c4f2a563c |
| `raw-recovery-manifest.json` | 49beb212b82342f13767d5f354645faedba5fc647b7dbfbfef82df80c92cb1c9 |
| `independent-summary.json` | 53912171c6fb678ef1079c9abdc15930b6ae1c24d470f02582543a3c6e6272a5 |
| `results.csv` | 79e036e53e9b35e758057e4a285dba64cc94c24ef8fded9d4629a0f4c71e887f |
| `review/manifest.json` | 85e3f2b202cb6b89acc181eb3d10e2465372d9d67587c3911c3ce6e4ee847990 |

Key review files: `review/raw-bad/` (13 raw failures),
`review/stage1-bad/1GC1-G-397-409-stage1.pdb` (the real residual failure),
and `residual-anchor-triage.json` (empty observed-anchor eligibility).
`verification.json` records unchanged original downloads, recovered raw/trace
bindings, zero observed movement and disabled anchor policies.

Verification after the raw/final fix: **1590 passed, 2 optional skips, 31 slow
tests deselected**; Ruff passes. Historical executed Protenix/Boltz driver hashes
remain unchanged; separate current-adapter hash fields track the code correction.
No new Boltz/Protenix neural runs or anchor-release experiments are claimed.
The commit-time focused regression selection additionally passes **95 tests**.
