# NCAC plus fixed relative rotamers: bounded accuracy improvement

Status: **retrospective eleven-mask conformation experiment**, not a public backend
or independent-family qualification. Continues the [whole-mask study](ncac-whole-mask.md):
NCAC repaired geometry on these cases but did not choose a good distal conformation.

## Protocol frozen before generating alternatives

`scripts/research_ncac_rotamers.py freeze` bound all eleven saved NCAC proposals,
source/implementation digests and exactly three relative rotations: 0/+120/−120°.
Both generated distal branches rotate rigidly together around CA–CB. Generated CB
and all observed atoms remain unchanged. No minimization, packing, hydrogen
optimization, native-derived angle tuning or score-weight fitting is introduced.
This is a three-state relative-χ experiment, **not an empirical rotamer library**.

Every candidate receives the unchanged raw/rounded/common stereo and steric audit.
All 33 proposals and three refusals are retained. Only the complete set of 30
geometry-feasible candidates enters the unchanged fixed-heavy PyRosetta runner:
`ref2015`, seed 7, revision `1a56185c2592611dec4c9c75ddc9468cd2227c1f`.
Fresh import rebuilds scorer-owned hydrogens without repacking or H optimization.
The runner-facing plan contains no native PDB or native-control cell. Missing
feasible scores block selection; stable candidate index breaks total-score ties.

The [independent selection ledger](ncac-relative-rotamers-evidence-selection.json)
is committed for all cases before the auditor reads native reference PDBs. The
[accuracy evidence](ncac-relative-rotamers-evidence.json) binds it by digest.
Source records, axial/radial transport, actual signed rotation angles, rigid
branch-pair distances, common audits and all 30 physical/topology import snapshots
are independently checked. The auditor does not regenerate or rescore candidates.
Scorer-owned terminal extras, virtual coordinates and context hydrogens are
identical across feasible alternatives in every case. Input heavy coordinates
remain unchanged, including generated CB. CA stays S and THR CB stays R.

## Results

All eleven pools yield a selected candidate. Selection changes four masks, all
THR; all four changes improve native RMSD and none worsen it on this dataset.

| Group | Masks | Changed / improved / worsened | Zero-angle median | Selected median | Prior own-scored median |
|---|---:|---:|---:|---:|---:|
| All | 11 | 4 / 4 / 0 | 0.322945 | **0.226284** | 0.293937 |
| Four previously seen NCAC mechanism cases | 4 | 2 / 2 / 0 | 1.113282 | 0.198433 | 0.157531 |
| Seven new only to NCAC | 7 | 2 / 2 / 0 | 0.322945 | 0.322945 | 0.383067 |

Numbers are matched rounded generated-atom exact-name RMSD medians in Å. Every
selected candidate and previous own-scored comparator passes the common gate;
the zero-angle outputs exactly repeat the original eleven-mask accuracy. The
seen-four group still has a worse median than prior own-scored. The unchanged
seven-case median does not imply that no individual case improved.

| Changed mask | Relative angle | Zero-angle RMSD | Selected RMSD |
|---|---:|---:|---:|
| 2CI2 THR22 | +120° | 1.945968 | 0.176830 |
| 1R69 THR26 | −120° | 1.848568 | 0.337043 |
| 1R69 THR27 | +120° | 1.970131 | 0.220035 |
| 1CTF THR76 | +120° | 2.061262 | 0.331298 |

Three rejected proposals remain in the geometry ledger and are never scored:
2CI2 VAL32 +120/−120° clash with I/33 OE2; 2CI2 THR55 −120° clashes with
I/69 CE2. Consult exact component identities in private geometry
evidence. Native VAL32/THR76 failures from the parent study remain failures;
this experiment neither retunes their thresholds nor drops those masks.

## Interpretation and limits

On this retained cohort, distal conformation selection explains substantial NCAC
THR errors without altering CB placement or observed coordinates. This is useful
mechanism evidence, not proof of broad accuracy, calibrated probabilities or
general superiority. Native outcomes were already known when this follow-up was
designed, despite frozen candidate generation and native-free score selection.
The next accuracy gate must use separately enrolled inputs and a frozen algorithm,
not more angle tuning on these same errors. Multi-χ/ILE support, empirical priors,
calibrated context chemistry and a full legacy comparator remain open.

Production dependencies, installed libraries, scientific thresholds, source
fixtures and consumers are unchanged. No public PDB/dat/bundle or GUI capability.
Authentic Dunbrack asset lineage remains unresolved: direct official-site HTTPS
connection and page fetches timed out; no bundled PyRosetta tables were extracted.

Retained artifacts under `tmp/pdbfixer-research/`:
`ncac-rotamers-protocol-01.json`, `ncac-rotamers-score-plan-01/`,
`ncac-rotamers-score-run-01/`. Earlier audit versions remain diagnostic history.

```bash
PYTHONPATH=src:. python scripts/audit_ncac_rotamers.py \
  --plan tmp/pdbfixer-research/ncac-rotamers-score-plan-01 \
  --run tmp/pdbfixer-research/ncac-rotamers-score-run-01 \
  --output <existing-new-directory>/ncac-relative-rotamers-evidence.json
```

Use fresh destination files. Focused suite: **689 passed in 26.46 s**; not a full
repository or cross-platform qualification run.
