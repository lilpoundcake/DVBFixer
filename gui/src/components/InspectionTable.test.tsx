import { renderToStaticMarkup } from 'react-dom/server'
import { TableBody, TableCell, TableHead, TableRow, ThemeProvider } from '@mui/material'
import { describe, expect, it } from 'vitest'
import { theme } from '../theme'
import { InspectionBadge, InspectionRow, InspectionTable, ResidueLabel } from './InspectionTable'

describe('shared inspection presentation', () => {
  it('renders common density, stripes, selected and keyboard-focus states', () => {
    const html = renderToStaticMarkup(<ThemeProvider theme={theme}>
      <InspectionTable aria-label="Inspector">
        <TableHead><TableRow><TableCell>Residue</TableCell></TableRow></TableHead>
        <TableBody><InspectionRow selected onClick={() => {}}><TableCell>
          <ResidueLabel chain="d" name="SER" number="82" insertionCode="A" atom="CB" />
        </TableCell></InspectionRow></TableBody>
      </InspectionTable>
    </ThemeProvider>)
    expect(html).toContain('height:32px')
    expect(html).toContain('nth-of-type(even)')
    expect(html).toContain(':focus-visible')
    expect(html).toContain('Mui-selected')
    expect(html).toContain('tabindex="0"')
    expect(html).toContain('d:SER82A · CB')
  })

  it('uses the same readable badge geometry for severity and interaction types', () => {
    for (const label of ['severe', 'Hydrogen Bond']) {
      const html = renderToStaticMarkup(<InspectionBadge label={label} color="#c62828" />)
      expect(html).toContain('height:18px')
      expect(html).toContain('font-size:10px')
      expect(html).not.toContain('text-transform:uppercase')
      expect(html).toContain(label)
    }
  })
})
