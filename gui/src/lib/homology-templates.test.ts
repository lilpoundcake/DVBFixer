import { describe, expect, it } from 'vitest'
import { bestMatchingChainId, chainIdentityLabel, pairwiseIdentityPercent, targetSequences } from './homology-templates'

describe('homology template helpers', () => {
  it('parses target chains and removes alignment/stop markers', () => {
    expect(targetSequences('>A target\nAC-D*\n>B\nGG?\n')).toEqual({ A: 'ACD', B: 'GG' })
  })

  it('computes global pairwise identity percentages', () => {
    expect(pairwiseIdentityPercent('ACDE', 'ACDE')).toBe(100)
    expect(pairwiseIdentityPercent('ACDE', 'ACDF')).toBe(75)
    expect(pairwiseIdentityPercent('ACDE', 'ACD')).toBe(75)
    expect(pairwiseIdentityPercent('', 'ACD')).toBeNull()
  })

  it('formats available and unavailable chain identities', () => {
    expect(chainIdentityLabel('A', 4, 'ACDE', 'ACDF')).toBe('A · 4 aa · 75.0% identity')
    expect(chainIdentityLabel('B', 3, '', 'ACD')).toBe('B · 3 aa · — identity')
  })

  it('selects the template chain with the highest target identity', () => {
    const chains = [{ id: 'A', sequence: 'AAAA' }, { id: 'B', sequence: 'ACDE' }]
    expect(bestMatchingChainId(chains, 'ACDF')).toBe('B')
    expect(bestMatchingChainId(chains, '', 'B')).toBe('B')
  })
})
