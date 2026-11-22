import { renderToStaticMarkup } from 'react-dom/server'
import { Box, ThemeProvider } from '@mui/material'
import { describe, expect, it } from 'vitest'
import { theme } from '../theme'
import { ALIGNMENT_LABEL_WIDTH, alignmentCellSx, alignmentRowActionsSx, alignmentSelectionSx } from './alignment'

describe('MSA residue selection presentation', () => {
  it('keeps three standard row actions inside the sticky label gutter', () => {
    expect(alignmentRowActionsSx).toMatchObject({ flexWrap: 'wrap', maxWidth: '100%', gap: 0.5 })
    // Three ordinary 80 px actions, two 4 px gaps and the label's 8 px side insets.
    expect(3 * 80 + 2 * 4 + 2 * 8).toBeLessThanOrEqual(ALIGNMENT_LABEL_WIDTH)
  })
  it('uses flat fills to distinguish transient selection from modeling spans', () => {
    expect(alignmentSelectionSx(false, true)).toEqual({ boxShadow: 'none', bgcolor: '#dce8fc' })
    expect(alignmentSelectionSx(true, true)).toEqual({ boxShadow: 'none', bgcolor: '#294f90' })
    for (const modelSelected of [false, true]) {
      expect(alignmentSelectionSx(modelSelected, false)).toEqual({ boxShadow: 'none' })
    }
  })

  it('preserves the fixed scientific grid and reserves outlines for keyboard focus', () => {
    expect(alignmentCellSx).toMatchObject({ width: 16, height: 22, flexShrink: 0, border: 0, outline: 'none' })
    expect(alignmentCellSx['&:focus-visible']).toEqual({ outline: '2px solid #263f78', outlineOffset: -2 })
    const html = renderToStaticMarkup(<ThemeProvider theme={theme}>
      <Box sx={{ display: 'flex' }}>{['A', 'C', 'D'].map(letter => <Box key={letter} component="button"
        sx={{ ...alignmentCellSx, ...alignmentSelectionSx(true, true) }}>{letter}</Box>)}</Box>
    </ThemeProvider>)
    expect(html).toContain('background-color:#294f90')
    expect(html).toContain('box-shadow:none')
    expect(html).toContain(':focus-visible')
    expect(html).not.toContain('inset 0 0 0 2px')
  })
})
