"""Causal pair-force ablation on a retained research_pdbfixer hydrogen run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    results = {}
    for seed in (1, 7, 19, 42, 101):
        source = args.run/f'h-seed{seed}-iter50-strip0-ff1'
        atoms = json.loads((source/'atoms.json').read_text())
        index = {a[3]: i for i, a in enumerate(atoms) if a[1] == '126'}
        initial = np.load(source/'stage-0.npy')
        system_text = (source/'system.xml').read_text()
        pair = [index['HG'], index['OXT']]
        for intervention in ('baseline', 'remove-pair-coulomb', 'repulsive-core'):
            system = mm.XmlSerializer.deserialize(system_text)
            nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
            qh, sh, eh = nb.getParticleParameters(pair[0])
            qo, so, eo = nb.getParticleParameters(pair[1])
            if any({int(nb.getExceptionParameters(i)[0]), int(nb.getExceptionParameters(i)[1])} == set(pair) for i in range(nb.getNumExceptions())):
                raise RuntimeError('Pair already excepted; ablation assumption invalid')
            if intervention == 'remove-pair-coulomb':
                nb.addException(*pair, 0.0, (sh+so)/2, (eh*eo)**0.5)
            elif intervention == 'repulsive-core':
                force = mm.CustomBondForce('k*(r0/r)^12')
                force.addGlobalParameter('k', 1.0)
                force.addGlobalParameter('r0', 0.15)
                force.addBond(*pair, [])
                system.addForce(force)
            integrator = mm.VerletIntegrator(0.0)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            context.setPositions(initial*unit.nanometer)
            mm.LocalEnergyMinimizer.minimize(context, 1.0, 500)
            state = context.getState(getPositions=True, getEnergy=True)
            xyz = np.asarray(state.getPositions().value_in_unit(unit.nanometer))
            name = f'seed{seed}-{intervention}'
            np.save(args.output/f'{name}.npy', xyz)
            results[name] = {'seed': seed, 'intervention': intervention,
                'hg_charge_e': float(qh.value_in_unit(unit.elementary_charge)),
                'oxt_charge_e': float(qo.value_in_unit(unit.elementary_charge)),
                'hg_epsilon_kj_mol': float(eh.value_in_unit(unit.kilojoule_per_mole)),
                'pair_epsilon_kj_mol': float(((eh*eo)**0.5).value_in_unit(unit.kilojoule_per_mole)),
                'og_oxt_angstrom': float(np.linalg.norm(xyz[index['OG']]-xyz[index['OXT']])*10),
                'hg_og_angstrom': float(np.linalg.norm(xyz[index['HG']]-xyz[index['OG']])*10),
                'hg_oxt_angstrom': float(np.linalg.norm(xyz[pair[0]]-xyz[pair[1]])*10),
                'energy_kj_mol': float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))}
            del context, integrator
    (args.output/'evidence.json').write_text(json.dumps(results, indent=2, allow_nan=False))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
