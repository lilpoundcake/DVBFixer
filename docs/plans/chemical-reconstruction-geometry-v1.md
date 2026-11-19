# WP1: authority-bound geometry proposals and independent audits

Status: bounded private WP1 implemented and checked on the persistent
`feature/chemical-reconstruction`; general scientific qualification remains open.
User authorized continuation after WP0; no production integration or new native
scientific cohort is part of this slice.

Documentation plan: this decision, active handoff continuation, task-map discovery
and a separate WP1 validation JSON; narrowly ignore retained test reports. These
four bounded records accompany implementation/tests. Pending links: none.

## Geometry authority decision

Use the **explicit selected state's supplied ideal graph** as a geometry authority
for this private slice. Bind its existing provenance, graph digest and declared
units through the unchanged shared request digest. This is not a new calibrated
protein geometry library, a residue-name lookup or permission to copy held-out
experimental coordinates into an authority. Synthetic tests label their supplied
graphs as synthetic. Real scientific qualification requires independently sourced,
documented geometry data; PyRosetta production/data-license gates stay unchanged.

Internal-coordinate plans explicitly name the generated atom and its three ordered
references. The existing `protein_completion.geometry.extend` supplies chain
construction; `two_angle_tip` supplies two-anchor cone intersection with an explicit
ordered stereo-frame constraint. Lengths come from graph bonds, angles and chain
torsion from ideal geometry. No `l_cb` default constants, empirical priors, native
reference scoring, optimization or stochastic fallback is added. Sequential plans
can describe arbitrary admitted missing subsets, not just one/two atoms. The first
strategy supports protein and isolated-component categories; that is a versioned
strategy envelope, not the shared architecture's chemistry limit.

A second adapter reuses the existing Class A `_fit_coordinates` rigid fit, without
calling authority resolution, microstate selection or its legacy geometry judge.
It remains the existing connected CHNO in-memory adapter subset. All absent graph
atoms must be owned by the request; it must not build unowned H and discard them.
Retain the complete shared request/context around the target-only fit and audit
every returned observation against its unchanged input. NumPy is imported lazily
only for this explicit strategy; WP0 imports/admission stay toolkit-free.
The existing Class A DTO requires pH metadata: this adapter requires it explicitly
rather than inventing a value. It is retained transport metadata only; neither
legacy authority/microstate resolution nor pH adjustment is called.

## Execution and ownership

Plans bind request digest, explicit audit tolerances and placement order. The first
executors refuse temporary/mobile ownership and generated atoms on external links:
those need a later linked/temporary construction policy, not strip-and-splice or
an unreported context deletion. Other context and declared links remain retained.
No new state is chosen and no existing coordinate is returned as an addition.

Retain generated-only raw proposals even when validation rejects them. Degenerate
construction reports its reason and any partial additions instead of reflecting
an observation or trying another engine. There is no conformer ranking in WP1.

Independent audits do not call the constructor: validate exact generated mask,
identity/element/unit/finite coordinates and request binding; measure supplied graph
bonds/angles and every explicit constraint; reject severe generated-atom overlaps
against the **whole retained context**. Topological 1–2/1–3 pairs are excluded from
the simple minimum-separation check, including declared covalent links. The caller's
positive minimum separation is only a severe-overlap guard, NOT calibrated VDW,
1–4 chemistry, coordination energetics or a global steric score. Unknown retained
chemistry stays explicitly unassessed.
Tolerances have no inferred/calibrated defaults. Minimum separation must be at
least `1e-6` Å to keep the bounded spatial-grid arithmetic numerically usable;
that engineering lower bound is not a chemical contact cutoff.

Repeat the audit with only generated coordinates rounded to PDB three-decimal
precision; preserve observations exactly, in their original units/order. Reject
unrepresentable generated coordinates and require both raw and rounded passes.
This is a numerical rounding control, not a source PDB writer or proof that headers,
links and every observation survive a serialized file. No bytes are published.

Successful execution means **declared geometry passed this private audit**, not
chemical valence/CIP correctness, global canonical L-chirality assessment, native
conformation, general H-placement qualification or MD readiness. Keep geometry-
approval/MD/source-coverage flags false. Any fully observed source defect detected
by the declared checks is a refusal, not permission to repair that source.

## Verification and gates

Synthetic protein cone, multi-step missing subset, Class A rigid-transform and
explicit-H controls; broken authority, masks/dependencies, wrong/inverted proposals,
context clashes, degenerate anchors and rounding-induced failure controls. Preserve
WP0 pure-import guards and old protein/Class A evidence regressions. No real
structure sampling or new claims for all20/rings/backbone/glycans/metals follow.
Continue the parallel heavy/H workstream in WP2 only after this bounded slice is
checked. The full engine/qualification roadmap remains in the active handoff.

## Implementation owners and private usage

| Owner | Boundary |
|---|---|
| [`proposals.py`](../../src/dvbfixer/chemical_reconstruction/proposals.py) | Immutable plans, explicit tolerances, generated-only candidates and retained outcomes |
| [`construction.py`](../../src/dvbfixer/chemical_reconstruction/construction.py) | Explicit internal/cone or Class A rigid strategy; no optimizer/fallback |
| [`authority_geometry.py`](../../src/dvbfixer/chemical_reconstruction/authority_geometry.py) | Nonmutating angstrom views of the supplied graph/unit authority |
| [`candidate_audit.py`](../../src/dvbfixer/chemical_reconstruction/candidate_audit.py) | Independent mask/binding/geometry/context audit |
| [`geometry controls`](../../tests/test_chemical_reconstruction_geometry.py) | Positive and injected-negative synthetic controls, not native qualification |

Given a typed shared `request` and explicitly ordered `placements` whose authority
has already been resolved by the caller:

```python
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.proposals import AuditPolicy, GeometryPlan, Strategy

# Example bounds for a synthetic control only, not recommended chemical defaults.
policy = AuditPolicy(0.002, 0.2, 0.2, 0.01)
plan = GeometryPlan(request.digest, Strategy.INTERNAL, policy, placements)
result = propose(request, plan)
assert not result.geometry_approved
assert not result.md_ready
```

**Do not** inspect only `coordinates_produced`: failed/partial proposals deliberately
retain coordinates. **Do** inspect status and both independent audits and preserve
raw/rejected results. `declared-geometry-passed` is not a publication, parameterization
or scientific qualification decision. Rounding rejection never triggers widening
of bounds or editing of observations.

The [WP1 validation record](chemical-reconstruction-wp1-validation.json) binds the
frozen source/test hashes and actual checks. [Active handoff](chemical-reconstruction-handoff.md)
owns the next work; WP0's earlier validation remains a historical source snapshot.

Verification: 34 new synthetic controls; 465 focused regressions and 3400 full-fast
tests passed, three existing optional Torch/Biotite skips and 33 slow deselected.
Strict mypy passes for all nine private modules; Ruff and agent maps pass. This
does not qualify the supplied ideal graphs as a real chemical geometry library.
