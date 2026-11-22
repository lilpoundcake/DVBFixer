export const ALIGNMENT_LABEL_WIDTH = 270

/** Keep whole row actions inside the sticky label, never over comparison cells. */
export const alignmentRowActionsSx = {
  display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: 0.5, maxWidth: '100%',
} as const

/** Scientific residue cells keep their grid geometry, not toolbar-button sizing. */
export const alignmentCellSx = {
  width: 16, height: 22, p: 0, border: 0, boxSizing: 'border-box', flexShrink: 0,
  fontFamily: 'inherit', fontSize: 12, outline: 'none',
  '&:focus-visible': {
    outline: '2px solid #263f78', outlineOffset: -2,
  },
} as const

/** A continuous fill, never a separate inset frame around every selected letter. */
export function alignmentSelectionSx(modelSelected: boolean, transientSelected: boolean) {
  return {
    boxShadow: 'none',
    ...(transientSelected ? { bgcolor: modelSelected ? '#294f90' : '#dce8fc' } : {}),
  } as const
}
