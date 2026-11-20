/** Fixed readable cells; only the number of residues per line changes on resize. */
export const SEQUENCE_LAYOUT = {
  blockSize: 10,
  cellWidth: 12,
  blockGap: 8,
  gutterWidth: 36,
} as const

/** Width is the sequence content box, excluding its scroll padding. */
export function sequenceResiduesPerLine(width: number): number {
  const { blockSize, cellWidth, blockGap, gutterWidth } = SEQUENCE_LAYOUT
  const available = Math.max(0, width - gutterWidth)
  const blocks = Math.floor((available + blockGap) / (blockSize * cellWidth + blockGap))
  // Very narrow docks still show whole cells, without shrinking text.
  return blocks > 0 ? blocks * blockSize : Math.max(1, Math.floor(available / cellWidth))
}
