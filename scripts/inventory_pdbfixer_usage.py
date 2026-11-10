"""Inventory direct PDBFixer dependencies without importing scientific pipelines.

Static AST inventory follows explicit imports and the current fixer/_fixer
receiver convention, not arbitrary Python dataflow. Comments/docstrings do not
count as calls. Optional stock-template probing imports PDBFixer only and forbids
network access; it demonstrates discovery scope, not successful reconstruction.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
from pathlib import Path
from unittest.mock import patch


def scan_source(source: str, path: str) -> dict:
    tree = ast.parse(source)
    constructors, rebuilders = set(), set()
    imports = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.module and (node.module == 'pdbfixer' or node.module.startswith('pdbfixer.')):
            imports.extend({'module': node.module, 'name': name.name, 'line': node.lineno} for name in node.names)
            constructors.update(name.asname or name.name for name in node.names if name.name == 'PDBFixer')
        if node.module == 'dvbfixer.ffutils.geometry':
            rebuilders.update(name.asname or name.name for name in node.names if name.name == 'rebuild_missing_atoms_with_retry')
    records = []

    class Visitor(ast.NodeVisitor):
        def __init__(self):
            self.scope = []

        def visit_FunctionDef(self, node):
            self.scope.append(node.name)
            self.generic_visit(node)
            self.scope.pop()

        visit_AsyncFunctionDef = visit_FunctionDef

        def emit(self, node, kind, name):
            records.append({'kind': kind, 'name': name, 'line': node.lineno,
                            'function': '.'.join(self.scope) or '<module>'})

        def visit_Call(self, node):
            fn = node.func
            if isinstance(fn, ast.Name):
                if fn.id in constructors:
                    self.emit(node, 'constructor', 'PDBFixer')
                if fn.id in rebuilders:
                    self.emit(node, 'seeded_rebuild', 'rebuild_missing_atoms_with_retry')
            if isinstance(fn, ast.Attribute) and isinstance(fn.value, ast.Name) and fn.value.id in {'fixer', '_fixer'}:
                self.emit(node, 'method', fn.attr)
            self.generic_visit(node)

        def visit_Attribute(self, node):
            if isinstance(node.value, ast.Name) and node.value.id in {'fixer', '_fixer'}:
                if node.attr in {'topology', 'positions', 'missingAtoms', 'missingTerminals',
                                 'missingResidues', 'nonstandardResidues', 'platform'}:
                    self.emit(node, 'state_write' if isinstance(node.ctx, ast.Store) else 'state_read', node.attr)
            self.generic_visit(node)

    Visitor().visit(tree)
    return {'path': path, 'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'imports': imports, 'uses': sorted(records, key=lambda record: (record['line'], record['kind']))}


def inventory(root: Path) -> dict:
    files = []
    for path in sorted((root/'src/dvbfixer').rglob('*.py')):
        result = scan_source(path.read_text(), str(path.relative_to(root)))
        if result['imports'] or result['uses']:
            files.append(result)
    return {'schema': 'dvbfixer.pdbfixer-usage.v1',
            'limits': 'Explicit imports and fixer/_fixer receivers; not a general call-graph or runtime coverage proof.',
            'files': files,
            'totals': {kind: sum(use['kind'] == kind for file in files for use in file['uses'])
                       for kind in ('constructor', 'seeded_rebuild', 'method')}}


def probe_templates() -> dict:
    import inspect
    from importlib.metadata import version

    from openmm import app
    from pdbfixer import PDBFixer
    from pdbfixer.pdbfixer import dnaResidues, proteinResidues, rnaResidues

    def no_network(*args, **kwargs):
        raise AssertionError('Template inventory must not download CCD or other data')

    pdb = ('ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n'
           'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C\n'
           'ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C\nEND\n')
    with patch('pdbfixer.pdbfixer.urlopen', no_network):
        fixer = PDBFixer(pdbfile=io.StringIO(pdb))
        stock = sorted(fixer.templates)
        trials = {}
        for name in stock:
            template = fixer.templates[name]
            text = io.StringIO()
            app.PDBFile.writeFile(template.topology, template.positions, text, keepIds=True)
            atom_lines = [line for line in text.getvalue().splitlines(True) if line.startswith(('ATOM  ', 'HETATM'))]
            # One actual template heavy atom is withheld; this is inventory only.
            removed = atom_lines[-1][12:16].strip()
            candidate = PDBFixer(pdbfile=io.StringIO(''.join(atom_lines[:-1])+'TER\nEND\n'))
            candidate.findMissingResidues()
            candidate.missingResidues = {}
            candidate.findMissingAtoms()
            detected = sorted(a.name for atoms in candidate.missingAtoms.values() for a in atoms)
            terminals = sorted(a for atoms in candidate.missingTerminals.values() for a in atoms)
            trials[name] = {'removed_atom': removed, 'detected_missing_atoms': detected,
                            'detected_terminal_atoms': terminals,
                            'withheld_atom_detected': removed in detected or removed in terminals}
        lookup = {name: fixer._getTemplate(name) is not None for name in (
            'I', 'DI', 'MSE', 'NHE', 'HOH', 'NAG', 'DAN', 'HEM', 'ZN', 'NLN', 'OLS', 'OLT', 'UNK')}
    return {'pdbfixer_version': version('pdbfixer'),
            'source_sha256': hashlib.sha256(Path(inspect.getfile(PDBFixer)).read_bytes()).hexdigest(),
            'stock_templates': stock, 'protein_keep_classification': proteinResidues,
            'rna_keep_classification': rnaResidues, 'dna_keep_classification': dnaResidues,
            'stock_template_lookup': lookup, 'deletion_detection_controls': trials,
            'network_allowed': False, 'coordinates_generated': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--probe-templates', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    result = inventory(root)
    if args.probe_templates:
        result['template_probe'] = probe_templates()
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
