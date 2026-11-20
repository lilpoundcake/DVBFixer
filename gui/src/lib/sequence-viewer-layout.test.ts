import { describe, expect, it } from 'vitest'
import { SEQUENCE_LAYOUT, sequenceResiduesPerLine } from './sequence-viewer-layout'

const { cellWidth, blockSize, blockGap, gutterWidth } = SEQUENCE_LAYOUT
const lineWidth = (count: number) => gutterWidth + count * cellWidth + Math.floor((count - 1) / blockSize) * blockGap

describe('adaptive sequence wrapping', () => {
  it('wraps at complete ten-residue blocks including gutter and inter-block gaps', () => {
    for (const count of [10, 20, 40, 60, 100, 160]) {
      expect(sequenceResiduesPerLine(lineWidth(count))).toBe(count)
      expect(sequenceResiduesPerLine(lineWidth(count) - 0.25)).toBe(count === 10 ? 9 : count - 10)
      expect(sequenceResiduesPerLine(lineWidth(count) + 0.25)).toBe(count)
    }
  })

  it('uses the actual dock width rather than retaining a 60-residue limit', () => {
    expect(sequenceResiduesPerLine(500)).toBe(30)
    expect(sequenceResiduesPerLine(350)).toBe(20)
    expect(sequenceResiduesPerLine(1000)).toBe(70)
    expect(sequenceResiduesPerLine(500)).toBe(30)
  })

  it('falls back to whole readable cells in very narrow or hidden docks', () => {
    expect(sequenceResiduesPerLine(100)).toBe(5)
    expect(sequenceResiduesPerLine(0)).toBe(1)
    for (let width = gutterWidth + cellWidth; width <= 2000; width++) {
      const count = sequenceResiduesPerLine(width)
      expect(lineWidth(count)).toBeLessThanOrEqual(width)
      if (count >= blockSize) expect(count % blockSize).toBe(0)
    }
  })
})
