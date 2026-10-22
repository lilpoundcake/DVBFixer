import { expect, it } from 'vitest'
import { parsePDB } from 'molstar/lib/mol-io/reader/pdb/parser'
import { trajectoryFromPDB } from 'molstar/lib/mol-model-formats/structure/pdb'
import { Structure, StructureProperties, StructureSelection, QueryContext } from 'molstar/lib/mol-model/structure'
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler'
import { altlocSerialExpression, altlocStateOpacity, displayedAltloc } from './altloc-preview'
import type { AltlocRow } from './altloc-types'

it('Selected only uses the tapped state, never the checked output state', () => {
  const row = { selected: 'A', states: [{ label: 'A' }, { label: 'B' }] } as AltlocRow
  expect(displayedAltloc(row, 'B')).toBe('B')
  expect(displayedAltloc(row, '')).toBe('A')
  expect(row.selected).toBe('A')
})

it('Compare makes only the tapped alternative opaque, regardless of the output selection', () => {
  const row = { selected: 'A', states: [{ label: 'A' }, { label: 'B' }, { label: 'C' }] } as AltlocRow
  const opacities = (label: string) => row.states.map(state => altlocStateOpacity('compare', state.label, displayedAltloc(row, label)))
  expect(opacities('B')).toEqual([0.25, 1, 0.25])
  expect(opacities('C')).toEqual([0.25, 0.25, 1])
  row.selected = 'C'
  expect(opacities('B')).toEqual([0.25, 1, 0.25])
  expect(altlocStateOpacity('all', 'A', 'B')).toBe(1)
})

it('Mol* retains A/B and uses exact Python serial masks rather than author residue guesses', async () => {
  const atom = (serial: number, label: string, x: number) =>
    `ATOM  ${String(serial).padStart(5)}  CA ${label}ALA C  82A   ${x.toFixed(3).padStart(8)}   2.000   3.000  0.50 20.00           C  \n`
  const parsed = await parsePDB(atom(11, 'A', 1) + atom(22, 'B', 2) + 'END\n').run()
  if (parsed.isError) throw new Error(parsed.message)
  const trajectory = await trajectoryFromPDB(parsed.result).run()
  const structure = Structure.ofModel(trajectory.representative)
  expect(structure.elementCount).toBe(2)
  const selection = compile(altlocSerialExpression([22]))(new QueryContext(structure))
  const selected = StructureSelection.unionStructure(selection)
  expect(selected.elementCount).toBe(1)
  const location = StructureElementLocation(selected)
  expect(StructureProperties.atom.label_alt_id(location)).toBe('B')
  expect(StructureProperties.atom.id(location)).toBe(22)
})

// Read a selected element without involving WebGL or a browser.
import { StructureElement } from 'molstar/lib/mol-model/structure'
function StructureElementLocation(structure: Structure) {
  const unit = structure.units[0]
  return StructureElement.Location.create(structure, unit, unit.elements[0])
}
