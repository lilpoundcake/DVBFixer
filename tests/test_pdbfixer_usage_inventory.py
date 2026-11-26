"""Stdlib-only audit regressions; never import scientific command parsers."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def scanner():
    spec = importlib.util.spec_from_file_location(
        'pdbfixer_usage_inventory', ROOT/'scripts/inventory_pdbfixer_usage.py',
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_aliases_count_but_comments_and_strings_do_not(scanner):
    source = '''
from pdbfixer import PDBFixer as Constructor
from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry as rebuild
def work():
    fixer = Constructor(filename='source.pdb')
    fixer.findMissingAtoms()
    fixer.missingResidues = {}
    rebuild(fixer)
    # fixer.addMissingHydrogens()
    message = 'fixer.addSolvent()'
'''
    records = scanner.scan_source(source, 'example.py')['uses']
    assert [(r['kind'], r['name']) for r in records] == [
        ('constructor', 'PDBFixer'), ('method', 'findMissingAtoms'),
        ('state_write', 'missingResidues'),
        ('seeded_rebuild', 'rebuild_missing_atoms_with_retry'),
    ]
    assert all(r['function'] == 'work' for r in records)


def test_current_ten_rebuild_sites_are_explicit(scanner):
    result = scanner.inventory(ROOT)
    assert result['totals']['constructor'] == 12
    sites = [(f['path'], u['line']) for f in result['files'] for u in f['uses']
             if u['kind'] == 'seeded_rebuild']
    assert len(sites) == 10
    assert {path for path, _ in sites} == {
        'src/dvbfixer/prepare/pipeline.py', 'src/dvbfixer/minimize/pipeline.py',
        'src/dvbfixer/protonate.py', 'src/dvbfixer/top/pipeline.py',
        'src/dvbfixer/model/diffusion/observed_completion.py',
        'src/dvbfixer/model/diffusion/boundary_refinement.py',
        'src/dvbfixer/model/diffusion/rfdiffusion_v1.py',
        'src/dvbfixer/model/diffusion/atom_completion_benchmark.py',
        'src/dvbfixer/chemical_reconstruction/template_comparison.py',
    }


def test_raw_rebuild_is_centralized_and_pdbfixer_hydrogens_unused(scanner):
    methods = [(f['path'], u['name']) for f in scanner.inventory(ROOT)['files']
               for u in f['uses'] if u['kind'] == 'method']
    raw = [path for path, name in methods if name == 'addMissingAtoms']
    assert raw == ['src/dvbfixer/ffutils/geometry.py']*2
    assert not {'addMissingHydrogens', 'addSolvent', 'addMembrane',
                'downloadTemplate', 'registerTemplate', 'removeChains'} & {name for _, name in methods}


def test_diagnostics_only_discovers_without_building(scanner):
    record = next(f for f in scanner.inventory(ROOT)['files']
                  if f['path'] == 'src/dvbfixer/diagnose/pipeline.py')
    assert not any(u['kind'] == 'seeded_rebuild' for u in record['uses'])
    assert {u['name'] for u in record['uses'] if u['kind'] == 'method'} == {
        'findMissingResidues', 'findMissingAtoms', 'findMissingTerminals',
    }


def test_retained_source_inventory_is_current(scanner):
    # New pure metadata adapter factoring is a separately retained amendment;
    # no old production/scientific source inventory is overwritten.
    retained = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory-metadata-chirality.json').read_text())
    historical = ROOT/retained['historical_inventory']['path']
    assert hashlib.sha256(historical.read_bytes()).hexdigest() == retained['historical_inventory']['sha256']
    assert scanner.inventory(ROOT)['files'] == retained['files']


def prior_to_comparison_totals(scanner):
    """Historical tests compare their historical scope, not the new private arm."""
    current = scanner.inventory(ROOT)
    files = [f for f in current['files'] if f['path'] != 'src/dvbfixer/chemical_reconstruction/template_comparison.py']
    previous = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory-glycam-fallback.json').read_bytes())
    assert files == previous['files']
    return {kind: sum(u['kind'] == kind for f in files for u in f['uses']) for kind in current['totals']}


def test_comparison_amendment_adds_only_private_seeded_baseline(scanner):
    path = ROOT/'docs/research/pdbfixer-usage-inventory-comparison.json'
    current = json.loads(path.read_bytes())
    previous_path = ROOT/current['historical_inventory']['path']
    previous = json.loads(previous_path.read_bytes())
    assert current['historical_inventory']['sha256'] == hashlib.sha256(previous_path.read_bytes()).hexdigest()
    actual = scanner.inventory(ROOT)
    assert [f for f in actual['files'] if f['path'] != 'src/dvbfixer/chemical_reconstruction/template_comparison.py'] == [
        f for f in current['files'] if f['path'] != 'src/dvbfixer/chemical_reconstruction/template_comparison.py']
    before = next(f for f in current['files'] if f['path'] == 'src/dvbfixer/chemical_reconstruction/template_comparison.py')
    after = next(f for f in actual['files'] if f['path'] == before['path'])
    for key in ('imports', 'uses'):
        assert [{k: v for k, v in row.items() if k != 'line'} for row in before[key]] == [
            {k: v for k, v in row.items() if k != 'line'} for row in after[key]]
    assert current['totals'] == scanner.inventory(ROOT)['totals']
    assert prior_to_comparison_totals(scanner) == previous['totals']
    assert current['added_private_sources'] == ['src/dvbfixer/chemical_reconstruction/template_comparison.py']
    assert current['added_totals'] == {'constructor': 1, 'seeded_rebuild': 1, 'method': 1}
    assert not current['historical_template_probe_rerun'] and not current['historical_scientific_results_changed']
    assert not current['public_backend_or_dependency_replacement']
    assert hashlib.sha256((ROOT/'scripts/amend_pdbfixer_usage_comparison.py').read_bytes()).hexdigest() == current['generator_source_sha256']


def test_metadata_chirality_amendment_is_exact_current_inventory_with_no_use_delta(scanner):
    path = ROOT/'docs/research/pdbfixer-usage-inventory-metadata-chirality.json'
    current = json.loads(path.read_bytes())
    parent = ROOT/current['historical_inventory']['path']
    assert hashlib.sha256(parent.read_bytes()).hexdigest() == current['historical_inventory']['sha256']
    previous = json.loads(parent.read_bytes())
    assert current['files'] == scanner.inventory(ROOT)['files']
    assert current['totals'] == previous['totals'] == scanner.inventory(ROOT)['totals']
    assert [r['path'] for r in current['changed_source_bindings']] == ['src/dvbfixer/chemical_reconstruction/template_comparison.py']
    for old, new in zip(previous['files'], current['files'], strict=True):
        assert old['path'] == new['path']
        for key in ('imports', 'uses'):
            assert [{k: v for k, v in r.items() if k != 'line'} for r in old[key]] == [
                {k: v for k, v in r.items() if k != 'line'} for r in new[key]]
    assert not current['historical_template_probe_rerun'] and not current['historical_scientific_results_changed']
    assert not current['public_backend_or_dependency_replacement']
    generator = ROOT/'scripts/amend_pdbfixer_usage_metadata_chirality.py'
    assert hashlib.sha256(generator.read_bytes()).hexdigest() == current['generator_source_sha256']


def test_glycam_amendment_keeps_oh_inventory_and_usage_semantics(scanner):
    path = ROOT/'docs/research/pdbfixer-usage-inventory-hydrogen-postflight.json'
    old = json.loads(path.read_bytes())
    new = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory-glycam-fallback.json').read_bytes())
    assert new['historical_inventory']['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert new['totals'] == old['totals'] == prior_to_comparison_totals(scanner)
    assert [c['path'] for c in new['changed_source_bindings']] == ['src/dvbfixer/prepare/pipeline.py']
    for a, b in zip(old['files'], new['files'], strict=True):
        assert a['path'] == b['path']
        for key in ('imports', 'uses'):
            assert [{k:v for k, v in r.items() if k != 'line'} for r in a[key]] == [
                {k:v for k, v in r.items() if k != 'line'} for r in b[key]]
    assert not new['historical_template_probe_rerun']
    assert not new['historical_scientific_results_changed']
    assert hashlib.sha256((ROOT/'scripts/amend_pdbfixer_usage_glycam_fallback.py').read_bytes()).hexdigest() == new['generator_source_sha256']


def test_hydrogen_amendment_preserves_usage_semantics_and_historical_evidence(scanner):
    historical_path = ROOT/'docs/research/pdbfixer-usage-inventory-stereo-diagnostics.json'
    historical = json.loads(historical_path.read_text())
    current = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory-hydrogen-postflight.json').read_text())
    assert current['historical_inventory']['sha256'] == hashlib.sha256(historical_path.read_bytes()).hexdigest()
    assert current['totals'] == historical['totals'] == prior_to_comparison_totals(scanner)
    assert {c['path'] for c in current['changed_source_bindings']} <= {
        'src/dvbfixer/ffutils/geometry.py', 'src/dvbfixer/prepare/pipeline.py', 'src/dvbfixer/diagnose/structural.py',
    }
    for before, after in zip(historical['files'], current['files'], strict=True):
        assert before['path'] == after['path']
        for key in ('imports', 'uses'):
            before_semantic = [{k:v for k, v in row.items() if k != 'line'} for row in before[key]]
            after_semantic = [{k:v for k, v in row.items() if k != 'line'} for row in after[key]]
            assert before_semantic == after_semantic
    assert not current['historical_template_probe_rerun']
    assert not current['historical_scientific_results_changed']
    assert current['production_hydrogen_postflight_changed']
    generator = ROOT/'scripts/amend_pdbfixer_usage_hydrogen.py'
    assert hashlib.sha256(generator.read_bytes()).hexdigest() == current['generator_source_sha256']


def test_inventory_amendment_changes_only_diagnostic_source_hash(scanner):
    old = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory.json').read_text())
    current = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory-stereo-diagnostics.json').read_text())
    assert current['totals'] == old['totals'] == prior_to_comparison_totals(scanner)
    assert [c['path'] for c in current['changed_source_bindings']] == ['src/dvbfixer/diagnose/pipeline.py']
    assert len(old['files']) == len(current['files'])
    for before, after in zip(old['files'], current['files'], strict=True):
        assert before['path'] == after['path']
        assert before['imports'] == after['imports'] and before['uses'] == after['uses']
    assert not current['historical_template_probe_rerun'] and not current['scientific_results_changed']
    generator = ROOT/'scripts/amend_pdbfixer_usage_inventory.py'
    assert hashlib.sha256(generator.read_bytes()).hexdigest() == current['generator_source_sha256']


def test_stock_keep_sets_do_not_imply_reconstruction_templates():
    probe = json.loads((ROOT/'docs/research/pdbfixer-usage-inventory.json').read_text())['template_probe']
    assert len(probe['stock_templates']) == 30
    assert set(probe['protein_keep_classification']) <= set(probe['stock_templates'])
    assert 'I' in probe['rna_keep_classification'] and not probe['stock_template_lookup']['I']
    assert 'DI' in probe['dna_keep_classification'] and not probe['stock_template_lookup']['DI']
    assert not any(probe['stock_template_lookup'].values())
    assert len(probe['deletion_detection_controls']) == 30
    assert all(r['withheld_atom_detected'] for r in probe['deletion_detection_controls'].values())
    assert not probe['network_allowed'] and not probe['coordinates_generated']
