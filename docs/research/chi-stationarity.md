# Finite local χ probes: residual descent at saved stopping-control endpoints

Status: **native-free finite sensitivity diagnosis**, not a new optimizer,
proposal selector, analytic derivative validation or convergence certificate.
Follows the [stopping-policy controls](chi-stopping-controls.md); the
[research-only licensing and production gate](../plans/pdbfixer-reconstruction-handoff.md#current-authorization-installed-pyrosetta-for-research-only-2026-11-03)
is unchanged.

## Frozen scope and numerical policy

`research_chi_stationarity.py plan` froze **all 42 mode-C endpoints** from all
15 original executable pools, including both parent geometry rejections. The
seven unsupported original masks remain inherited nonexecutions. Mode C is used
as the final predeclared stopping control, **not** as a profile chosen for native
accuracy. No native coordinates, RMSD reporting or candidate selection.

For every target χ (including proton χ), score independent endpoint clones at
**±0.01° and ±0.005°**. No sequential stepping or optimization. The objective,
weights, residue types and topology remain unchanged. The run keeps every probe,
not just downhill ones, and never adopts a perturbed state as an output.

- Initial complete pose snapshots and numerical score metadata must reproduce
  exactly, with JSON container normalization only.
- Materialize the saved endpoint by setting **only its mobile atoms** from its
  saved trace. Protected atoms are never moved and restored. Require the entire
  resulting atom trace to match the saved endpoint exactly, with unchanged graph.
- Recompute its energy within a predeclared **1e−8 absolute score-unit bound**.
  This is saved-coordinate materialization, **not** bit-identical final-score
  replay or another minimizer execution.
- Each signed perturbation must preserve protected coordinates within 1e−9 Å
  and preserve graph identity; actual protected movement is zero.
- Report a resolved finite descent only when its score change is **below −1e−7**
  units. This numerical reporting floor is not a scientific acceptance threshold.

For each step size, report central slope `(Eplus−Eminus)/(2h)` and curvature
`((Eplus−E0)+(Eminus−E0))/h²`, with `h` in radians. Both step sizes remain visible;
their difference is reported, not used to adapt a step-size grid.

## Results

[Independent saved-evidence readback](chi-stationarity-evidence.json):

| Measure | Result |
|---|---:|
| Enrolled / materialized / audited endpoints | 42 / 42 / 42 |
| Execution or protected-motion blockers | 0 |
| Target χ coordinates / signed perturbations | 62 / 248 |
| Endpoints with resolved finite descent | **19/42** |
| Endpoints with descent along a heavy-moving χ | **14/42** |
| Protected maximum movement | **0.0 Å** |
| Largest absolute central slope | 0.425411 score units/radian |
| Most negative signed energy change | −0.0000725595 score units |
| Maximum saved/materialized endpoint score difference | 5.68e−14 score units |

The heavy-moving subset is identified from retained atom movement, not inferred
solely from χ numbering: at least one generated heavy atom moves more than
1e−9 Å under that probe. The 19 and 14 counts overlap; they are not additive.
The remaining 23 endpoints are **not certified stationary**: absence of resolved
descent in these four finite tests per χ is not a convergence proof.

These results establish that nearby lower-score states remain at 19 saved mode-C
endpoints under this finite diagnostic. They do not identify whether iteration
limits, failed line searches, analytic derivatives, coordinate parameterization
or another numerical mechanism caused termination. In particular, a central
finite slope is not automatically the exact gradient quantity tested inside
Rosetta. No direct analytic-gradient comparison or iteration trace was executed.

## Ownership, geometry and limits

The same private ownership policy applies: generated distal heavy atoms and
target sidechain H may move; generated CB, all observations, backbone H, context
H/extras/virtuals stay protected. Exact saved endpoint coordinates and every full
probe atom trace are retained privately. Source observations are never restored
after movement. Independent readback re-derives allowed mobility and verifies
all 248 atom inventories and protected movement.

These probes are energy diagnostics, **not geometry-qualified proposals**. Their
raw/rounded common geometry was not individually qualified, and no new claim of
probe chirality or geometry acceptance is made. Parent mode-C geometry results
(40/42) remain separate and unchanged. Including its two rejected endpoints is
intentional and is not evidence that they became acceptable. There is no new
candidate ranking, output structure or best-mode selection.

All observations here are repeated numerical tests on previously evaluated
inputs, not independent scientific samples or accuracy qualification. No native
RMSD, geometry cutoffs, weights, starting angles or tolerance settings were used
to select these probes or altered afterward.

## Reproduction and verification

Private plan/run: `tmp/pdbfixer-research/chi-stationarity-plan-01.json` and
`chi-stationarity-run-01/`. `audit_chi_stationarity.py` does not invoke PyRosetta,
minimize, generate coordinates or evaluate scores. It binds the complete endpoint
and signed-probe inventory, rechecks full atom traces and independently computes
finite differences. No library means/sigmas/raw weights/tables are published.

```bash
PYTHONPATH=src:. python scripts/audit_chi_stationarity.py \
  --plan tmp/pdbfixer-research/chi-stationarity-plan-01.json \
  --run tmp/pdbfixer-research/chi-stationarity-run-01 \
  --output <fresh-report.json>
```

**751 focused tests passed in 26.38 s**; targeted Ruff and 17-file private-package
mypy pass. This is not the full repository/platform suite. Production consumers,
defaults, dependencies, fixtures and public capability remain unchanged; no commit
or push. Original tables/provenance/license rights and explicit user approval
remain mandatory before production, separately from scientific qualification.

Subsequent [logging-only exact replay](chi-termination.md) identifies explicit
iteration-cap warnings for two of these endpoints. The other 17 finite-descent
endpoints remain unexplained; no new numerical profile or accuracy claim follows.
