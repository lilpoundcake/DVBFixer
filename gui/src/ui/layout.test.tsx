import { renderToStaticMarkup } from 'react-dom/server'
import { ThemeProvider } from '@mui/material/styles'
import Button from '@mui/material/Button'
import IconButton from '@mui/material/IconButton'
import TextField from '@mui/material/TextField'
import Box from '@mui/material/Box'
import AddIcon from '@mui/icons-material/Add'
import { describe, expect, it } from 'vitest'
import { theme } from '../theme'
import { UI, actionGroupSx, commandActionGroupSx, commandInputSx, compactIconButtonSx, formActionRowSx, inspectionTableSx, panelToolbarSx, primaryActionSx, workflowContentSx, workflowControlSx, workflowStageSx, workflowToolbarSx, workspaceRunColors, workspaceRowColors } from './layout'

describe('shared GUI control geometry', () => {
  it('compacts whitespace without shrinking existing text or controls', () => {
    expect(theme.typography.body2.fontSize).toBe('0.75rem')
    expect(theme.typography.caption.fontSize).toBe('0.7rem')
    expect(theme.components?.MuiButton?.styleOverrides?.root).toMatchObject({ fontSize: '0.75rem', height: 32 })
    expect(theme.components?.MuiInputBase?.styleOverrides?.root).toMatchObject({ fontSize: '0.75rem' })
    expect(UI.toolbarMinHeight).toBe(40)
  })
  it('separates run headers from nested file stripes without borrowing selection colors', () => {
    expect(workspaceRunColors[0].header).not.toBe(workspaceRunColors[1].header)
    for (const band of workspaceRunColors) {
      expect(band.files[0]).not.toBe(band.files[1])
      expect(band.files).not.toContain(band.header)
      expect(band.header).not.toBe(workspaceRowColors.selected)
      expect(band.files).not.toContain(workspaceRowColors.selected)
    }
  })
  it('shares Interactions row density, alternating backgrounds and distinct selection styling', () => {
    expect(inspectionTableSx['& > thead > tr > th, & > tbody > tr > td']).toMatchObject({ px: 1, py: 0.5, fontSize: '0.75rem', lineHeight: '18px' })
    const row = inspectionTableSx['& > tbody > tr']
    expect(row.height).toBe(28)
    expect(row.bgcolor).toBe('#ffffff')
    expect(row['&:nth-of-type(even)'].bgcolor).toBe('#f8fafc')
    expect(inspectionTableSx['& > thead > tr > th, & > tbody > tr > td']).not.toHaveProperty('bgcolor')
    expect(inspectionTableSx['& > thead > tr > th:not(:last-child), & > tbody > tr > td:not(:last-child)']).toEqual({
      borderRight: '1px solid', borderRightColor: 'rgba(29,50,97,0.07)',
    })
    expect(row['&.MuiTableRow-hover:hover'].bgcolor).toBe('action.hover')
    expect(row['&.Mui-selected'].bgcolor).toBe('rgba(29,50,97,0.24)')
    expect(row['&.Mui-selected:hover'].bgcolor).toBe('rgba(29,50,97,0.31)')
    expect(inspectionTableSx['& > tbody > tr > td.MuiTableCell-paddingCheckbox'].py).toBe(0)
  })
  it('uses one height for buttons, icon buttons, toggles and input controls', () => {
    expect(theme.components?.MuiButton?.styleOverrides?.root).toMatchObject({ height: UI.controlHeight, minHeight: UI.controlHeight })
    expect(theme.components?.MuiIconButton?.styleOverrides?.root).toMatchObject({ width: UI.controlHeight, height: UI.controlHeight })
    expect(theme.components?.MuiToggleButton?.styleOverrides?.root).toMatchObject({ height: UI.controlHeight })
    expect(theme.components?.MuiOutlinedInput?.styleOverrides?.root).toMatchObject({ '&:not(.MuiInputBase-multiline)': { height: UI.controlHeight } })
    expect(theme.components?.MuiInputLabel?.styleOverrides?.sizeSmall).toMatchObject({
      transform: 'translate(14px, 6px) scale(1)',
      '&.MuiInputLabel-shrink': { transform: 'translate(14px, -8px) scale(0.75)' },
    })
  })

  it('separates outlined buttons and toggle states from the toolbar surface', () => {
    expect(theme.palette.background.paper).not.toBe('#ffffff')
    expect(theme.components?.MuiButton?.styleOverrides?.outlined).toMatchObject({
      backgroundColor: '#ffffff',
      '&:hover': { backgroundColor: '#f7f9fc' },
      '&.Mui-disabled': { backgroundColor: '#f7f8fa' },
    })
    expect(theme.components?.MuiToggleButton?.styleOverrides?.root).toMatchObject({
      backgroundColor: '#ffffff',
      '&.Mui-selected': { backgroundColor: '#dce5f4', color: theme.palette.primary.main },
      '&.Mui-selected:hover': { backgroundColor: '#cedcf0' },
      '&.Mui-disabled': { backgroundColor: '#f7f8fa' },
    })
    expect(theme.components?.MuiOutlinedInput?.styleOverrides?.root).toMatchObject({
      backgroundColor: '#ffffff',
      '&.Mui-disabled': { backgroundColor: '#f7f8fa' },
    })
    expect(workflowToolbarSx.bgcolor).toBe('#ffffff')
    expect(panelToolbarSx.bgcolor).toBe('background.paper')
    expect(panelToolbarSx.py * 8).toBeGreaterThanOrEqual(4)
  })

  it('keeps workflow actions 160×32 and dense row actions explicitly 24×24', () => {
    expect(primaryActionSx).toMatchObject({ width: 160, height: 32 })
    expect(workflowToolbarSx).not.toHaveProperty('& .MuiButton-root')
    expect(compactIconButtonSx).toMatchObject({ width: 24, height: 24, minWidth: 24 })
  })

  it('shares compact workflow insets without stacking tab, stage and toolbar whitespace', () => {
    expect(workflowContentSx).toMatchObject({ px: 1, py: 0.5, minHeight: 0, overflow: 'auto' })
    expect(workflowStageSx).toEqual({ display: 'flex', flexDirection: 'column', gap: 0.5 })
    expect(workflowToolbarSx).toMatchObject({ px: 0, py: 0.5, minHeight: 40 })
    // 4 px content inset + 4 px toolbar inset above; 4 px toolbar + 4 px gap below.
    expect((workflowContentSx.py + workflowToolbarSx.py) * 8).toBe(8)
    expect((workflowToolbarSx.py + workflowStageSx.gap) * 8).toBe(8)
  })

  it('lets command inputs yield width before wrapping the right-aligned action group', () => {
    expect(commandInputSx).toMatchObject({ flex: '1 1 160px', minWidth: 'min(160px, 100%)', maxWidth: '100%' })
    expect(commandActionGroupSx).toMatchObject({ ml: 'auto', flexShrink: 0, maxWidth: '100%', justifyContent: 'flex-end', flexWrap: 'wrap' })
    expect(primaryActionSx.width).toBe(160)
  })

  it('lets workflow fields yield width while reserving 160 px only for primary actions', () => {
    expect(workflowControlSx).toMatchObject({ flex: '1 1 120px', minWidth: 'min(120px, 100%)', maxWidth: 160 })
    expect(primaryActionSx).toMatchObject({ width: 160, height: 32 })
    expect(theme.components?.MuiButton?.styleOverrides?.root).toMatchObject({ minWidth: 80, height: 32 })
    // Three short fields, one main action and the Model row's three 8 px gaps.
    expect(3 * 120 + primaryActionSx.width + 3 * 8).toBeLessThanOrEqual(560)
  })

  it('wraps toolbars and action groups rather than shifting controls with transforms', () => {
    expect(panelToolbarSx).toMatchObject({ alignItems: 'center', flexWrap: 'wrap', gap: 0.5, px: 1, minHeight: 40, py: 0.5 })
    expect(panelToolbarSx['&:has(.MuiInputLabel-root)']).toEqual({ py: 1, rowGap: 1, minHeight: 48 })
    expect(workflowToolbarSx['&:has(.MuiInputLabel-root)']).toBe(panelToolbarSx['&:has(.MuiInputLabel-root)'])
    expect(actionGroupSx).toMatchObject({ ml: 'auto', alignItems: 'center', flexWrap: 'wrap' })
    expect(formActionRowSx).toMatchObject({ alignItems: 'flex-start', gap: 1 })
    expect(formActionRowSx).not.toHaveProperty('transform')
  })

  it('renders labelled and disabled controls with shared sizing and visible keyboard focus', () => {
    const markup = renderToStaticMarkup(<ThemeProvider theme={theme}>
      <Box sx={panelToolbarSx}>
        <TextField label="Input structure" value="" onChange={() => {}} />
        <Button variant="contained" startIcon={<AddIcon />}>Import files</Button>
        <IconButton aria-label="Add workspace" disabled><AddIcon /></IconButton>
      </Box>
    </ThemeProvider>)
    expect(markup).toContain('height:32px')
    expect(markup).toContain('aria-label="Add workspace"')
    expect(markup).toContain('disabled=""')
    expect(markup).toContain('Mui-focusVisible')
    expect(markup).toContain('outline-offset:2px')
    expect(markup).toContain('Input structure')
  })
})
