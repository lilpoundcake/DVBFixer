import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder'
import type { AltlocRow } from './altloc-types'

/** Compare highlights the tapped state, independently of output checkboxes. */
export function altlocStateOpacity(mode: 'all' | 'selected' | 'compare', label: string, tapped: string): number {
  return mode === 'compare' && label !== tapped ? 0.25 : 1
}

/** Viewing a tapped state is independent of the checkbox used for saving. */
export function displayedAltloc(row: AltlocRow, tappedLabel: string): string {
  return row.states.some(state => state.label === tappedLabel) ? tappedLabel : row.states[0]?.label || ''
}

/** Serial masks come from Python's normalized, model-specific inventory. */
export function altlocSerialExpression(serials: readonly number[]) {
  return MS.struct.generator.atomGroups({
    'atom-test': MS.core.set.has([MS.set(...serials), MS.struct.atomProperty.macromolecular.id()]),
  })
}
