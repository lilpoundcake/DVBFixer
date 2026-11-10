"""Private stage-traced PDBFixer/OpenMM experiments; never a production backend.

Run with PYTHONPATH=src and an intact scientific interpreter. Output reservation
is exclusive. All monkey patches are process-local and restored. Raw rebuilding
is invoked only by the production seeded wrapper; rejected attempts are retained.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import inspect
import io
import json
import platform
import random
import sys
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import numpy as np
import openmm as mm
import pdbfixer.pdbfixer as pf
from openmm import app, unit

from dvbfixer.ffutils.geometry import (
    build_ca_chirality_force,
    detect_misplaced_hydrogens,
    find_d_residues,
    rebuild_missing_atoms_with_retry,
)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def key(atom):
    r = atom.residue
    return [r.chain.id, r.id, r.insertionCode, atom.name]


def volumes(topology, positions):
    xyz = np.asarray(positions.value_in_unit(unit.nanometer))
    result = []
    for r in topology.residues():
        a = {a.name: a.index for a in r.atoms()}
        if r.name != 'GLY' and all(n in a for n in ('N', 'CA', 'C', 'CB')):
            n, ca, c, cb = (xyz[a[n]] for n in ('N', 'CA', 'C', 'CB'))
            result.append([r.chain.id, r.id, r.insertionCode, r.name,
                           float(np.dot(np.cross(n-ca, c-ca), cb-ca))])
    return result


def proper_overlay(points1, points2):
    """Independent proper rigid fit using row-vector covariance, returned column form."""
    p, q = np.asarray(points1), np.asarray(points2)
    if len(p) < 3:
        raise ValueError('research proper fit requires >=3 anchors')
    u, _, vt = np.linalg.svd((q-q.mean(0)).T @ (p-p.mean(0)))
    correction = np.diag([1., 1., np.linalg.det(u @ vt)])
    return -q.mean(0), (u @ correction @ vt).T, p.mean(0)


class Trace:
    def __init__(self, directory):
        self.directory = directory
        self.events = []
        self.topology = None

    def snapshot(self, label, positions, context=None):
        xyz = np.asarray(positions.value_in_unit(unit.nanometer))
        number = len(self.events)
        np.save(self.directory / f'{number:03d}-{label}.npy', xyz)
        atoms = list(self.topology.atoms())
        event = {'stage': label, 'coordinates': f'{number:03d}-{label}.npy',
                 'volumes_nm3': volumes(self.topology, positions),
                 'd_residues': find_d_residues(self.topology, positions)}
        (self.directory / f'{number:03d}-atoms.json').write_text(json.dumps([key(a) for a in atoms]))
        if context is not None:
            state = context.getState(getEnergy=True, getForces=True)
            event['energy_kj_mol'] = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            event['max_force_kj_mol_nm'] = float(np.linalg.norm(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer), axis=1).max())
        self.events.append(event)


def heavy_case(text, directory, proper=False, guard=False, backbone_only=False):
    directory.mkdir()
    (directory/'input.pdb').write_text(text)
    trace = Trace(directory)
    overlays = []
    iterations = []
    original_overlay = pf._overlayPoints
    original_minimize = mm.LocalEnergyMinimizer.minimize
    original_integrator = mm.LangevinIntegrator
    original_context = mm.Context

    def overlay(p, q):
        if backbone_only:
            # Only used for the declared N/CA/C/O-only sidechain cohort.
            if len(p) != 4:
                raise ValueError('backbone-only control requires exactly N/CA/C/O')
            p, q = p[:3], q[:3]
        out = proper_overlay(p, q) if proper and len(p) >= 3 else original_overlay(p, q)
        overlays.append({'anchors': len(p), 'determinant': float(np.linalg.det(out[1]))})
        return out

    class Fixer(pf.PDBFixer):
        def _addAtomsToTopology(self, *args):
            result = super()._addAtomsToTopology(*args)
            trace.topology = result[0]
            trace.snapshot('construction', result[1])
            return result

    class Integrator(original_integrator):
        def step(self, steps):
            super().step(steps)
            trace.snapshot('dynamics', self._context.getState(getPositions=True).getPositions(), self._context)

    def minimize(context, *args):
        class Reporter(mm.MinimizationReporter):
            def report(self, iteration, x, grad, values):
                xyz = np.asarray(x).reshape(-1, 3)
                iterations.append({'iteration': iteration,
                                   'volumes_nm3': volumes(trace.topology, xyz*unit.nanometer)})
                return False

        context.getIntegrator()._context = context
        trace.snapshot('before-minimize', context.getState(getPositions=True).getPositions(), context)
        original_minimize(context, *args, reporter=Reporter())
        trace.snapshot('after-minimize', context.getState(getPositions=True).getPositions(), context)

    def make_context(system, *args):
        if guard:
            force, count = build_ca_chirality_force(trace.topology, fixer._research_initial)
            system.addForce(force)
            fixer._research_guard_count = count
        (directory/'system.xml').write_text(mm.XmlSerializer.serialize(system))
        return original_context(system, *args)

    fixer = Fixer(pdbfile=io.StringIO(text), platform=mm.Platform.getPlatformByName('Reference'))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    # Explicit gaps are supplied by the experiment, not a second sequence mapper.
    if 'REMARK 999 RESEARCH GAP10-13' in text:
        fixer.missingResidues = {(0, 9): ['GLY', 'LYS', 'THR', 'ILE']}
        fixer.findMissingAtoms()
    with ExitStack() as stack:
        stack.enter_context(patch.object(pf, '_overlayPoints', overlay))
        # A construction-only preview supplies guard reference coordinates.
        preview = fixer._addAtomsToTopology(True, True)
        fixer._research_initial = preview[1]
        stack.enter_context(patch.object(mm, 'LangevinIntegrator', Integrator))
        stack.enter_context(patch.object(mm, 'Context', make_context))
        stack.enter_context(patch.object(mm.LocalEnergyMinimizer, 'minimize', minimize))
        try:
            rebuild_missing_atoms_with_retry(fixer, max_attempts=1)
            outcome = 'returned'
        except Exception as exc:
            outcome = f'{type(exc).__name__}: {exc}'
    with (directory/'returned.pdb').open('w') as handle:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    result = {'input_sha256': digest(text.encode()), 'proper_fit': proper, 'guard': guard,
              'backbone_only': backbone_only, 'minimizer_iterations': iterations,
              'outcome': outcome, 'overlays': overlays, 'events': trace.events,
              'final_d': find_d_residues(fixer.topology, fixer.positions)}
    (directory/'result.json').write_text(json.dumps(result, indent=2))
    return result


def hydrogen_case(text, directory, seed, iterations=50, strip=False, forcefield=True):
    directory.mkdir()
    if strip:
        text = ''.join(l for l in text.splitlines(True) if not (l.startswith('ATOM') and l[76:78].strip() == 'H'))
    (directory/'input.pdb').write_text(text)
    pdb = app.PDBFile(io.StringIO(text))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    original = mm.LocalEnergyMinimizer.minimize
    states = []

    def minimize(context, tolerance, maximum):
        states.append(np.asarray(context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer)))
        original(context, tolerance, iterations)
        states.append(np.asarray(context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer)))
        (directory/'system.xml').write_text(mm.XmlSerializer.serialize(context.getSystem()))

    random.seed(seed)
    with patch.object(mm.LocalEnergyMinimizer, 'minimize', minimize):
        modeller.addHydrogens(app.ForceField('amber14-all.xml') if forcefield else None,
                              platform=mm.Platform.getPlatformByName('Reference'))
    for i, xyz in enumerate(states):
        np.save(directory/f'stage-{i}.npy', xyz)
    atoms = list(modeller.topology.atoms())
    (directory/'atoms.json').write_text(json.dumps([key(a) for a in atoms]))
    findings = detect_misplaced_hydrogens(modeller.topology, modeller.positions)
    xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
    ser = {a.name: a.index for a in atoms if a.residue.id == '126'}
    original_xyz = {tuple(key(a)): np.asarray(pdb.positions[a.index].value_in_unit(unit.nanometer)) for a in pdb.topology.atoms()}
    max_fixed = max(float(np.linalg.norm(xyz[a.index]-original_xyz[tuple(key(a))])) for a in atoms if tuple(key(a)) in original_xyz)
    result = {'seed': seed, 'iterations': iterations, 'strip': strip, 'forcefield': forcefield,
              'input_sha256': digest(text.encode()), 'max_existing_movement_nm': max_fixed,
              'hg_og_angstrom': float(np.linalg.norm(xyz[ser['HG']]-xyz[ser['OG']])*10),
              'hg_oxt_angstrom': float(np.linalg.norm(xyz[ser['HG']]-xyz[ser['OXT']])*10),
              'findings': [{'atom': key(f.h_atom), 'reason': f.reason, 'distance_nm': f.parent_distance_nm} for f in findings]}
    with (directory/'output.pdb').open('w') as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle, keepIds=True)
    (directory/'result.json').write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    resources = Path(pf.__file__).parent
    if not (resources/'soft.xml').is_file():
        raise RuntimeError('PDBFixer soft.xml missing: do not repair prefix piecemeal')
    evidence = {'python': sys.version, 'machine': platform.platform(), 'openmm': mm.__version__,
                'platform': 'Reference', 'script_sha256': digest(Path(__file__).read_bytes()),
                'sources': {}, 'heavy': {}, 'hydrogens': {}}
    for name, obj in [('pdbfixer', pf), ('modeller', app.modeller)]:
        source = Path(obj.__file__)
        evidence['sources'][name] = {'sha256': digest(source.read_bytes()), 'symbols': {}}
        (args.output/f'{name}-source.py').write_bytes(source.read_bytes())
    for name in ('_overlayPoints',):
        obj = getattr(pf, name)
        lines, start = inspect.getsourcelines(obj)
        evidence['sources']['pdbfixer']['symbols'][name] = [start, start+len(lines)-1]
    for name in ('findMissingResidues', 'findMissingAtoms', '_addAtomsToTopology', '_addMissingResiduesToChain', '_createForceField', 'addMissingAtoms'):
        lines, start = inspect.getsourcelines(getattr(pf.PDBFixer, name))
        evidence['sources']['pdbfixer']['symbols'][name] = [start, start+len(lines)-1]
    lines, start = inspect.getsourcelines(app.Modeller.addHydrogens)
    evidence['sources']['modeller']['symbols']['addHydrogens'] = [start, start+len(lines)-1]
    evidence['resources'] = {str(p.relative_to(resources)): digest(p.read_bytes()) for p in resources.rglob('*') if p.is_file() and p.suffix in ('.pdb', '.xml')}
    original = (root/'tests/fixtures/protein_context/1UBQ.pdb').read_text()
    evidence['fixture_sha256'] = digest(original.encode())
    subset = [l for l in original.splitlines(True) if l.startswith('ATOM') and l[21] == 'A' and int(l[22:26]) <= 72 and l[76:78].strip() != 'H']
    sidechain = ''.join(l for l in subset if not (int(l[22:26]) == 27 and l[12:16].strip() in ('CG', 'CD', 'CE', 'NZ')))+'TER\nEND\n'
    gap = 'REMARK 999 RESEARCH GAP10-13\n'+''.join(l for l in subset if not 10 <= int(l[22:26]) <= 13)+'TER\nEND\n'
    for case, text in [('lys27', sidechain), ('gap10-13', gap)]:
        for proper, guard in [(False, False), (True, False), (True, True)]:
            name = f'{case}-proper{int(proper)}-guard{int(guard)}'
            evidence['heavy'][name] = heavy_case(text, args.output/name, proper, guard)
    tree = ast.parse((root/'tests/test_prepare_broken_geom.py').read_text())
    text = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '_SER_INPUT_HG_MISSING_OXT_PRESENT' for t in n.targets))
    for seed in (1, 7, 19, 42, 101):
        for iterations, strip, ff in [(50, False, True), (500, False, True), (50, True, True), (50, False, False)]:
            name = f'h-seed{seed}-iter{iterations}-strip{int(strip)}-ff{int(ff)}'
            evidence['hydrogens'][name] = hydrogen_case(text, args.output/name, seed, iterations, strip, ff)
    (args.output/'evidence.json').write_text(json.dumps(evidence, indent=2, allow_nan=False))
    print(json.dumps({'heavy': {k: v['final_d'] for k,v in evidence['heavy'].items()}, 'hydrogens': evidence['hydrogens']}, indent=2))


if __name__ == '__main__':
    main()
