import type { AltlocRow } from './altloc-types'

/** UI state only. Policies, occupancy and atom selection are evaluated in Python. */
export function chooseAltloc(rows: readonly AltlocRow[], choices: Record<string, string>, key: string, label: string): Record<string, string> {
  if (!rows.find(row => row.key === key)?.states.some(state => state.label === label)) throw new Error('State is absent from this analysis')
  return { ...choices, [key]: label }
}

export function altlocSelectionComplete(rows: readonly AltlocRow[], choices: Record<string, string>): boolean {
  return Object.keys(choices).length === rows.length && rows.every(row => row.states.some(state => state.label === choices[row.key]))
}
