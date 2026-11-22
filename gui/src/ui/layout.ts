/** Shared desktop control geometry. Keep docs/gui-design-system.md in sync. */
export const UI = {
  controlHeight: 32,
  compactHeight: 24,
  iconSize: 18,
  compactIconSize: 14,
  actionWidth: 160,
  fieldWidth: 160,
  wideFieldWidth: 320,
  toolbarMinHeight: 40,
  contentGap: 4,
  floatingLabelOutset: 8,
} as const

export const toolbarSx = {
  display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: 0.5,
  px: 1, py: 0.5, minHeight: UI.toolbarMinHeight, boxSizing: 'border-box',
  flexShrink: 0, bgcolor: 'background.paper',
  // Floating labels extend above the input border. Reserve room inside the
  // colored band, including when controls wrap onto another row.
  '&:has(.MuiInputLabel-root)': {
    '--ui-label-clearance': `${UI.floatingLabelOutset}px`,
    py: UI.floatingLabelOutset / 8,
    rowGap: (UI.contentGap + UI.floatingLabelOutset) / 8,
    minHeight: UI.controlHeight + 2 * UI.floatingLabelOutset,
  },
} as const

export const panelToolbarSx = {
  ...toolbarSx, borderBottom: 1, borderColor: 'divider',
} as const

/** Workspace and scientific inspectors share hover/selection semantics. */
export const workspaceRowColors = {
  base: 'transparent', stripe: 'rgba(29,50,97,0.045)', hover: 'action.hover',
  selected: 'rgba(29,50,97,0.24)', selectedHover: 'rgba(29,50,97,0.31)',
} as const

/** Inspector data stays lighter than the tinted header, even on a tinted parent. */
export const inspectionRowColors = {
  ...workspaceRowColors, base: '#ffffff', stripe: '#f8fafc',
} as const

/** Two-level Workspace bands: run identity is independent of child row parity. */
export const workspaceRunColors = [
  { header: 'rgba(43,97,153,0.09)', files: ['background.paper', 'rgba(43,97,153,0.045)'], edge: 'rgba(43,97,153,0.55)' },
  { header: 'rgba(113,83,148,0.09)', files: ['background.paper', 'rgba(113,83,148,0.045)'], edge: 'rgba(113,83,148,0.55)' },
] as const

/** Shared compact inspection-table appearance; toolbar controls stay 32 px. */
export const inspectionTableSx = {
  '& > thead > tr > th, & > tbody > tr > td': {
    px: 1, py: 0.5, fontSize: '0.75rem', lineHeight: '18px',
  },
  '& > thead > tr > th': { fontWeight: 700, bgcolor: 'background.paper' },
  '& > thead > tr > th:not(:last-child), & > tbody > tr > td:not(:last-child)': {
    borderRight: '1px solid', borderRightColor: 'rgba(29,50,97,0.07)',
  },
  '& > tbody > tr > td[align="right"]': { fontVariantNumeric: 'tabular-nums' },
  '& > tbody > tr': {
    height: 28,
    bgcolor: inspectionRowColors.base,
    '&:nth-of-type(even)': { bgcolor: inspectionRowColors.stripe },
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

/** Floating labels extend above their fields; repeated rows need 12 px clearance. */
export const labelledFieldRowsSx = {
  display: 'flex', flexDirection: 'column', gap: (UI.contentGap + UI.floatingLabelOutset) / 8,
} as const

/** Reserve label space for direct fields/field rows, not arbitrary nested panels.
 * Toolbars declare the clearance already provided by their own top padding. */
export const formStackSx = {
  display: 'flex', flexDirection: 'column', gap: UI.contentGap / 8,
  '& > .MuiFormControl-root:has(> .MuiInputLabel-root), & > :has(> .MuiFormControl-root > .MuiInputLabel-root)': {
    marginTop: `calc(${UI.floatingLabelOutset}px - var(--ui-label-clearance, 0px))`,
  },
} as const

export const sequenceRecordRowSx = {
  display: 'grid', gridTemplateColumns: `${UI.controlHeight}px minmax(0, 1fr) max-content`,
  gap: 1, alignItems: 'center',
} as const

export const controlSx = { width: UI.fieldWidth, minWidth: UI.fieldWidth, maxWidth: '100%' } as const
export const wideControlSx = { width: UI.wideFieldWidth, minWidth: 0, maxWidth: '100%' } as const
/** Command inputs yield space to actions before the row needs to wrap. */
export const commandInputSx = {
  ...wideControlSx, flex: `1 1 ${UI.fieldWidth}px`, minWidth: `min(${UI.fieldWidth}px, 100%)`,
} as const
export const commandActionGroupSx = {
  ...actionGroupSx, flexShrink: 0, maxWidth: '100%', justifyContent: 'flex-end',
} as const
/** Short workflow fields can yield modest width without shrinking their text. */
export const workflowControlSx = {
  ...controlSx, flex: '1 1 120px', minWidth: 'min(120px, 100%)', maxWidth: UI.fieldWidth,
} as const
export const primaryActionSx = {
  width: UI.actionWidth, minWidth: UI.actionWidth, height: UI.controlHeight, flexShrink: 0,
} as const
export const controlWrapperSx = { display: 'inline-flex', alignItems: 'center' } as const
export const actionTooltipWrapperSx = { ...controlWrapperSx, ...primaryActionSx } as const

/** One shared content inset and stage gap; do not stack extra tab/toolbar gutters. */
export const workflowContentSx = {
  flex: 1, minHeight: 0, overflow: 'auto', px: 1, py: 0.5,
} as const

export const workflowStageSx = formStackSx

export const workflowToolbarSx = {
  ...toolbarSx,
  // The containing workflowContentSx already supplies the horizontal inset.
  px: 0,
  // Workflow controls belong to the white content canvas, not a tinted panel header.
  bgcolor: '#ffffff',
} as const

export const compactIconButtonSx = {
  width: UI.compactHeight, minWidth: UI.compactHeight, height: UI.compactHeight, minHeight: UI.compactHeight, p: 0.5,
  '& .MuiSvgIcon-root': { fontSize: UI.compactIconSize },
} as const

export const inlineActionSx = {
  height: UI.compactHeight, minHeight: UI.compactHeight, minWidth: 0, px: 1,
} as const
