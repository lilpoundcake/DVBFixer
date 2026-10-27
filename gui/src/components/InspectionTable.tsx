import { Box, Table, TableContainer, TableRow } from '@mui/material'
import type { TableProps, TableRowProps } from '@mui/material'
import { inspectionTableSx } from '../ui/layout'

/** Shared scrolling, density and interaction states for scientific inspectors. */
export function InspectionTable({ children, ...props }: Omit<TableProps, 'sx' | 'size' | 'stickyHeader'>) {
  return <TableContainer sx={{ flex: 1, minHeight: 0, overflow: 'auto' }}>
    <Table {...props} size="small" stickyHeader sx={inspectionTableSx}>{children}</Table>
  </TableContainer>
}

export function InspectionRow({ onClick, children, ...props }: TableRowProps) {
  return <TableRow {...props} hover onClick={onClick} tabIndex={onClick ? 0 : undefined}
    onKeyDown={event => {
      if (onClick && event.target === event.currentTarget && (event.key === 'Enter' || event.key === ' ')) {
        event.preventDefault()
        event.currentTarget.click()
      }
    }} sx={{ cursor: onClick ? 'pointer' : undefined }}>{children}</TableRow>
}

export function ResidueLabel({ chain, name, number, insertionCode = '', atom }: {
  chain: string; name: string; number: string | number; insertionCode?: string; atom?: string
}) {
  return <Box component="span" sx={{ fontFamily: 'ui-monospace, monospace', whiteSpace: 'nowrap' }}>
    {chain.trim() || '_'}:{name}{number}{insertionCode}{atom ? ` · ${atom}` : ''}
  </Box>
}

export function InspectionBadge({ label, color }: { label: string; color: string }) {
  return <Box component="span" sx={{ display: 'inline-flex', alignItems: 'center', height: 18,
    px: 0.75, borderRadius: 1, fontSize: '10px', lineHeight: '18px', fontWeight: 600,
    whiteSpace: 'nowrap', color, bgcolor: `${color}20` }}>{label}</Box>
}
