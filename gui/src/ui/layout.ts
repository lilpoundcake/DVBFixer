/** Shared desktop control geometry. Keep docs/gui-design-system.md in sync. */
export const UI = {
  controlHeight: 32,
  compactHeight: 24,
  iconSize: 18,
  compactIconSize: 14,
  actionWidth: 160,
  fieldWidth: 160,
  wideFieldWidth: 320,
  toolbarMinHeight: 48,
} as const

export const toolbarSx = {
  display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: 1,
  px: 1.5, py: 1, minHeight: UI.toolbarMinHeight, boxSizing: 'border-box',
  flexShrink: 0, bgcolor: 'background.paper',
} as const

export const panelToolbarSx = {
  ...toolbarSx, borderBottom: 1, borderColor: 'divider',
} as const

/** Workspace and scientific inspectors share the exact same row palette. */
export const workspaceRowColors = {
  base: 'transparent', stripe: 'rgba(29,50,97,0.045)', hover: 'action.hover',
  selected: 'rgba(29,50,97,0.24)', selectedHover: 'rgba(29,50,97,0.31)',
} as const

/** Shared inspection-table appearance, based on the Interactions density. */
export const inspectionTableSx = {
  '& > thead > tr > th, & > tbody > tr > td': {
    px: 1.5, py: 0.75, fontSize: '0.75rem', lineHeight: '18px',
  },
  '& > thead > tr > th': { fontWeight: 700, bgcolor: 'background.paper' },
  '& > tbody > tr > td[align="right"]': { fontVariantNumeric: 'tabular-nums' },
  '& > tbody > tr': {
    height: 32,
    bgcolor: workspaceRowColors.base,
    '&:nth-of-type(even)': { bgcolor: workspaceRowColors.stripe },
    '&.MuiTableRow-hover:hover': { bgcolor: workspaceRowColors.hover },
    '&.Mui-selected': { bgcolor: workspaceRowColors.selected },
    '&.Mui-selected:hover': { bgcolor: workspaceRowColors.selectedHover },
    '&:focus-visible': { outline: '2px solid', outlineColor: 'primary.main', outlineOffset: -2 },
  },
  '& > tbody > tr > td.MuiTableCell-paddingCheckbox': { px: 0.5, py: 0 },
} as const

export const actionGroupSx = {
  display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: 1, ml: 'auto',
} as const

/** Use when fields have helper text below them; actions align to the input, not its caption. */
export const formActionRowSx = {
  display: 'flex', alignItems: 'flex-start', flexWrap: 'wrap', gap: 1,
} as const

export const controlSx = { width: UI.fieldWidth, minWidth: UI.fieldWidth, maxWidth: '100%' } as const
export const wideControlSx = { width: UI.wideFieldWidth, minWidth: 0, maxWidth: '100%' } as const
export const primaryActionSx = {
  width: UI.actionWidth, minWidth: UI.actionWidth, height: UI.controlHeight, flexShrink: 0,
} as const
export const controlWrapperSx = { display: 'inline-flex', alignItems: 'center' } as const
export const actionTooltipWrapperSx = { ...controlWrapperSx, ...primaryActionSx } as const

export const workflowToolbarSx = {
  ...toolbarSx,
  '& .MuiButton-root': primaryActionSx,
} as const

export const compactIconButtonSx = {
  width: UI.compactHeight, minWidth: UI.compactHeight, height: UI.compactHeight, minHeight: UI.compactHeight, p: 0.5,
  '& .MuiSvgIcon-root': { fontSize: UI.compactIconSize },
} as const

export const inlineActionSx = {
  height: UI.compactHeight, minHeight: UI.compactHeight, minWidth: 0, px: 1,
} as const
