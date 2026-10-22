import { describe, expect, it } from 'vitest'
import { chooseAltloc, altlocSelectionComplete } from './altloc-selection'
import type { AltlocRow } from './altloc-types'

const row: AltlocRow = { key: '[1,"C","82","A",1]', model: 1, chain: 'C', residueNumber: '82', insertionCode: 'A', occurrence: 1,
  selected: null, method: 'manual', sharedSerials: [1], warnings: [], states: ['A', 'B'].map((label, index) => ({ label, residueName: 'SER', occupancy: .5, serials: [index + 2], atomNames: ['CB'], atomCount: 1 })) }

describe('altloc checkboxes', () => {
  it('keeps exactly one checked state and replaces rather than appends', () => {
    const first = chooseAltloc([row], {}, row.key, 'A')
    const second = chooseAltloc([row], first, row.key, 'B')
    expect(first[row.key]).toBe('A')
    expect(second).toEqual({ [row.key]: 'B' })
    expect(altlocSelectionComplete([row], second)).toBe(true)
  })
  it('blocks incomplete/stale selections without guessing or computing occupancy', () => {
    expect(altlocSelectionComplete([row], {})).toBe(false)
    expect(altlocSelectionComplete([row], { [row.key]: 'Z' })).toBe(false)
    expect(altlocSelectionComplete([row], { [row.key]: 'A', stale: 'B' })).toBe(false)
    expect(() => chooseAltloc([row], {}, row.key, 'Z')).toThrow('absent')
  })
})
