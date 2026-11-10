"""Fixed convenience-fixture enrollment for a graph-defined methyl control.

No coordinates parsed, source repaired, H placed or outcomes evaluated here.
Known occupancy refusal control stays enrolled; no replacement inputs permitted.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write

SOURCES = (
    ('1CRN', 'tests/fixtures/protein_context/scoring_followup/1CRN.pdb'),
    ('1CTF', 'tests/fixtures/protein_context/independent_followup/1CTF.pdb'),
    ('1TIM', 'tests/fixtures/protein_context/1TIM.pdb'),
    ('1UBQ', 'tests/fixtures/protein_context/1UBQ.pdb'),
    ('1EMV', 'tests/fixtures/regressions/1EMV.pdb'),
    ('1FR2', 'tests/fixtures/regressions/1FR2.pdb'),
    ('2VLN', 'tests/fixtures/regressions/2VLN.pdb'),
    ('2VLQ', 'tests/fixtures/regressions/2VLQ.pdb'),
)


def enroll(output):
    records = []
    inputs = {}
    for label, name in SOURCES:
        path = Path(name)
        digest = sha(path.read_bytes())
        inputs[name] = digest
        records.append(dict(id=label, source=name, source_sha256=digest))
    for name in ('scripts/enroll_methyl_diversity.py', 'scripts/prospective_independent_cohort.py',
                 'src/dvbfixer/protein_completion/source.py', 'src/dvbfixer/protein_completion/branched_source.py'):
        inputs[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.methyl-diversity-enrollment.v1', sources=records, inputs=inputs,
        policy=dict(
            timing='convenience input set fixed before graph-initializer implementation/source admission/H-control outcomes; sources previously used in other research, not blinded or temporal holdout',
            projection='unchanged project_first: first coordinate model and first ATOM chain, retain exact ATOM H/altloc/occupancy/numbering and all touching LINK/SSBOND/CONECT; original bytes retained separately',
            admission='reuse canonical heavy-only source layout/global CA/global THR-beta checks and declared-link endpoint validation; retain every refusal, no stripping H/altloc/occupancy or crop/replacement to gain entry',
            preparation='one legacy prepare per admitted projection, standard naming, unchanged PROPKA/Reduce/pH; PDBFixer and H placement explicitly Reference, H seed1; require resource readiness first',
            source_capture='two predeclared H seeds1,2, forced protein-H rebuild and Reference on prepared/dat; observe native H Context and first whole-System start without performing full native refinement',
            groups='all graph methyls: existing carbon with exactly three newly mobile H and one heavy bonded neighbor, canonical protein context; second-shell heavy frame chosen within residue first then full atom identity; no guessed frame',
            angular='canonical tetrahedral H angles, each original C-H radius; all nonmethyl initial coordinates and existing heavy atoms exact',
            arms=['unchanged', 'graph_tetrahedral_zero', 'graph_tetrahedral_energy_selected'],
            selection='all methyls phase0 first; energy-selected arm makes exactly one pass over parent identities in lexicographic case-sensitive order; each group evaluates 0..345 by15 degrees in original System energy, retains best, no revisit or final/native scoring',
            ties='1e-8 kJ/mol, shortest signed canonical phase then signed angle',
            physics='original Reference H-only tolerance1.0/max50; no force/exclusion/mass changes, H regeneration during arms, additional seed/profile/iteration search or integration',
            ownership='all three touched H must be absent from pre-H observations and mobile in original H System; no measured/kept H rewrite',
            audit='all admissions, groups, initial energy grids, callback paths, original endpoint equality and adverse geometry; broader whole-System continuation requires its own declared transport/control plan',
            limits='research-only new protein contexts, not original OpenCL THR2 replay, production default/repair, full-complex or broad chemistry qualification; 1UBQ remains an explicit known occupancy-scope control',
            replacement_sources_allowed=False, maximum_sources=8, H_seeds=[1, 2]),
        production_changed=False, native_coordinates_read=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    enroll(args.output)
