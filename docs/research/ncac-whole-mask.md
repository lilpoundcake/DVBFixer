# N–CA–C fitting on all eleven whole-sidechain masks

Status: **retrospective research, not a production repair or broad qualification**.
Extends the [four-case anchor controls](template-anchor-controls.md) to every
whole-sidechain VAL/THR mask in the existing 2CI2/1R69/1CTF cohort. Four cases had
already executed under NCAC; seven are new only to this variant. Original native
and baseline outcomes were already known: this is not a blinded independent cohort.

## Frozen protocol and audit order

`scripts/research_ncac_expansion.py` froze all eleven masks in original enrollment
order before reconstruction. Both stock and NCAC use normal relaxation and the
unchanged helper's seeds 1–5. NCAC uses the previous proper N–CA–C fit: observed O
is excluded from fitting, not moved or removed from the model. No reflection,
new score selector, threshold adjustment or installed-library patch was added.

The [geometry-only ledger](ncac-whole-mask-evidence-geometry.json) is written before
the auditor opens native reference PDBs. The subsequent
[accuracy evidence](ncac-whole-mask-evidence.json) binds that ledger by digest.
This audit-order separation does not undo earlier knowledge of cohort outcomes.
Raw and rounded accuracy are reporting only; they cannot choose a proposal/seed.
Own-scored comparisons replay the earlier fixed-heavy score choice, not a refit.

Independent readback verifies source identities/coordinates at every stage,
reported iterates, fit transport, proper used-anchor optimum, template/hybrid
geometry, and every helper-returned seed under the unchanged common geometry
gate. All 23 stock seed states reproduce the original untraced baseline exactly;
all four seen NCAC complete traces and payloads repeat exactly. There are 22 method
arms, 34 returned states, 170 stage snapshots, 67 target minimizer iterates,
22 checked NCAC fits and zero dynamics calls. All NCAC arms return on seed 1.

## Geometry coverage is not native accuracy

| Group | Masks | Stock passes | NCAC passes | Own-scored passes | Native passes |
|---|---:|---:|---:|---:|---:|
| All | 11 | 4 | **11** | 11 | 9 |
| Previously tested NCAC | 4 | 0 | 4 | 4 | 4 |
| New only to NCAC | 7 | 4 | 7 | 7 | 5 |

All NCAC constructions and returns have CA S; THR endpoints also have CB R.
Observed atoms stay fixed. These eleven passes are a bounded geometry result,
not proof that NCAC always prevents inversion or handles other missing-atom classes.

Native failures remain explicit: 2CI2 VAL32 CG1–C is 2.758456 Å against the
unchanged 2.8 Å cutoff; 1CTF THR76 CA–CB–CG2 is 125.118487°, outside the illustrative
tetrahedral ±12° gate. Neither native case is dropped or used to loosen a threshold.

Rounded generated-atom exact-name RMSD medians, Å:

| Matched accepted comparison | Count | NCAC | Other method |
|---|---:|---:|---:|
| NCAC vs own-scored, all | 11 | 0.322945 | 0.293937 |
| NCAC vs own-scored, seen four | 4 | 1.113282 | 0.157531 |
| NCAC vs own-scored, new seven | 7 | 0.322945 | 0.383067 |
| NCAC vs stock, mutual common passes | 4 | 0.218371 | 0.215798 |

The all-raw stock median 1.252770 Å includes invalid proposals and is **diagnostic**,
not an accepted-output accuracy comparison. There are no stock/NCAC mutual passes
in the seen-four stratum; its matched median is null, not zero. Generated-only VAL
symmetry medians equal the exact-name medians here; raw and rounded values remain
separate in evidence. These comparisons do not establish general superiority.

## Per-mask rounded exact-name RMSD

Stock numbers include rejected outputs; consult the geometry ledger for admission.

| Mask | Stock | NCAC | Own-scored |
|---|---:|---:|---:|
| 2CI2 THR22 | 1.771746 | 1.945968 | 0.146591 |
| 2CI2 VAL28 | 1.241830 | 0.280596 | 0.168472 |
| 2CI2 VAL32 | 0.343435 | 0.322945 | 0.568045 |
| 2CI2 THR55 | 1.273464 | 0.116132 | 0.207788 |
| 1R69 VAL6 | 0.113927 | 0.122943 | 1.910512 |
| 1R69 VAL24 | 0.193651 | 0.210458 | 0.293937 |
| 1R69 THR26 | 2.577022 | 1.848568 | 0.383067 |
| 1R69 THR27 | 2.490213 | 1.970131 | 0.137056 |
| 1CTF VAL56 | 1.252770 | 0.369833 | 0.114035 |
| 1CTF VAL66 | 0.237945 | 0.226284 | 0.352396 |
| 1CTF THR76 | 1.836512 | 2.061262 | 0.480019 |

Correct chirality does not settle distal conformation. Several NCAC THR proposals
remain roughly 2 Å from native despite passing geometry. Conversely own-scored
VAL6 is much worse than NCAC. Neither observation licenses native-based selection.
The next separate experiment should freeze relative rotamer alternatives and
score all geometry-feasible candidates without native information; it must not
alter this cohort's original proposals or claim independent enrollment.

## Reproduction and limits

Frozen artifacts: `tmp/pdbfixer-research/ncac-expansion-plan-01.json` and
`ncac-expansion-run-01/`. Original fixtures and previous evidence are unchanged.
No public PDB/`.dat`/bundle, GUI capability, dependency or production consumer changes.
Only whole-sidechain VAL/THR masks with exactly observed N/CA/C/O are covered, not
retained distal atoms, ILE/multi-χ chemistry, loops, PTMs or full legacy preparation.

```bash
PYTHONPATH=src:. python scripts/audit_ncac_expansion.py \
  --plan tmp/pdbfixer-research/ncac-expansion-plan-01.json \
  --run tmp/pdbfixer-research/ncac-expansion-run-01 \
  --output <new-directory>/ncac-whole-mask-evidence.json
```

Use fresh output paths; the companion geometry filename is derived automatically.
Both publications must replay byte-identically with the same basenames. Focused
suite: **678 passed in 27.11 s**; this is not full-repository/platform acceptance.
