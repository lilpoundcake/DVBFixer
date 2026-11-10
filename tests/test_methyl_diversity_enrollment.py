"""The fixed H-initializer cohort is enrolled without parsing/evaluating input."""
import hashlib
import json
from pathlib import Path

from scripts.enroll_methyl_diversity import SOURCES, enroll


def test_convenience_inputs_and_known_refusal_control_are_fixed(tmp_path):
    output = tmp_path/'enrollment.json'
    enroll(output)
    record = json.loads(output.read_bytes())
    assert [r['id'] for r in record['sources']] == ['1CRN', '1CTF', '1TIM', '1UBQ', '1EMV', '1FR2', '2VLN', '2VLQ']
    assert len(SOURCES) == 8
    for source in record['sources']:
        assert source['source_sha256'] == hashlib.sha256(Path(source['source']).read_bytes()).hexdigest()
    assert not record['policy']['replacement_sources_allowed']
    assert record['policy']['H_seeds'] == [1, 2]
    assert record['policy']['arms'] == ['unchanged', 'graph_tetrahedral_zero', 'graph_tetrahedral_energy_selected']
    assert not record['production_changed'] and not record['native_coordinates_read']
