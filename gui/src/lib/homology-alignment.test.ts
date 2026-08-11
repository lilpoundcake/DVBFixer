import { describe, expect, it } from 'vitest'
import { alignmentColumnsForResidues, alignmentColumnsToSpans, comparisonToReference, consensusFor, updateColumnSelection } from './homology-alignment'

describe('homology alignment consensus', () => {
  it('uses Clustal identity, strong, weak, gap, and mismatch marks', () => {
    expect(consensusFor([
      { sequence: 'ANCSG' },
      { sequence: 'AQAG-' },
    ])).toBe('*:.. ')
  })

  it('marks differences independently against the target row', () => {
    expect(comparisonToReference('ANCSG', 'AQWG-')).toBe('*:×. ')
    expect(comparisonToReference('ANCSG', 'ANCWG')).toBe('***×*')
  })

  it('extends from the newest additive anchor without clearing old ranges', () => {
    const firstRange = updateColumnSelection([1], 'ABCDEFG', 3, 1, 'extend')
    const withSecondAnchor = updateColumnSelection(firstRange, 'ABCDEFG', 6, 6, 'toggle')
    expect(updateColumnSelection(withSecondAnchor, 'ABCDEFG', 5, 6, 'extend')).toEqual([1, 2, 3, 5, 6])
  })

  it('converts only explicit transient columns into modeling spans', () => {
    expect(alignmentColumnsToSpans([5, 1, 2, 5, 7])).toEqual([
      { start: 1, end: 3 },
      { start: 5, end: 6 },
      { start: 7, end: 8 },
    ])
    expect(alignmentColumnsToSpans([])).toEqual([])
  })

  it('maps 3D residue ordinals to gapped alignment columns', () => {
    expect(alignmentColumnsForResidues('A-CD-E', new Set([2, 4]))).toEqual([2, 5])
    expect(alignmentColumnsForResidues('A-CD-E', new Set())).toEqual([])
  })
})
