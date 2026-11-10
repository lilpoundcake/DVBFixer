"""Native-reference files are excluded until complete all-arm selection commits."""
import ast
import json
from pathlib import Path

from scripts.audit_prospective_guard_cohort import ARMS, freeze_manifest


def test_manifest_binds_masked_sources_but_defers_native_coordinate_files(tmp_path):
    root=tmp_path/'masks'
    root.mkdir()
    for name in ('x-original.pdb','x-reference.pdb','masked.pdb','plan.json'):
        (root/name).write_text('synthetic transport control\n')
    output=tmp_path/'manifest.json'
    freeze_manifest(dict(masks=str(root)),output)
    record=json.loads(output.read_bytes())
    assert str(root/'masked.pdb') in record['inputs']
    assert str(root/'plan.json') in record['inputs']
    assert str(root/'x-original.pdb') not in record['inputs']
    assert str(root/'x-reference.pdb') not in record['inputs']


def test_no_native_accuracy_call_occurs_in_selection():
    tree=ast.parse(Path('scripts/audit_prospective_guard_cohort.py').read_text())
    selection=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='select')
    assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('accuracy','_parse_atom_lines') for n in ast.walk(selection))
    reporting=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='report')
    assert isinstance(reporting.body[0],ast.Assign)
    assert reporting.body[0].value.func.id=='select'


def test_all_predeclared_arms_are_reported_not_an_accuracy_chosen_subset():
    assert ARMS==('stock','own','mean_prior','mean_fixed','mean_c','mean_guard')
