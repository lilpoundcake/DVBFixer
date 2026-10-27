import { Box, Button, Checkbox, FormControl, InputLabel, MenuItem, Select, TableBody, TableCell, TableHead, TableRow, Tooltip, Typography } from '@mui/material'
import { useWorkspaceStore } from '../stores/workspaceStore'
import { useAltlocStore, EMPTY_ALTLOC_SESSION } from '../stores/altlocStore'
import { compactIconButtonSx, controlSx, panelToolbarSx, primaryActionSx } from '../ui/layout'
import { InspectionRow, InspectionTable, ResidueLabel } from './InspectionTable'
import { openAltlocReview } from '../lib/altloc-operations'

/** Manual review controls and table; routine CLI configuration stays in DVBfixer. */
export function AltlocTable() {
  const id = useWorkspaceStore(state => state.active?.id || '')
  const session = useAltlocStore(state => state.sessions[id] || EMPTY_ALTLOC_SESSION)
  const choose = useAltlocStore(state => state.choose)
  const focus = useAltlocStore(state => state.focus)
  return <Box sx={{ display: 'flex', flexDirection: 'column', height: '100%', minHeight: 0 }}>
    <Box sx={panelToolbarSx}>
      <FormControl size="small" sx={controlSx} disabled={!session.report}><InputLabel>3D display</InputLabel>
        <Select label="3D display" value={session.mode} onChange={event => {
          useAltlocStore.getState().patch(id, { mode: event.target.value as typeof session.mode })
          openAltlocReview(id)
        }}>
          <MenuItem value="all">All alternatives</MenuItem><MenuItem value="selected">Selected only</MenuItem><MenuItem value="compare">Compare</MenuItem>
        </Select>
      </FormControl>
      <Button variant="outlined" sx={{ ...primaryActionSx, ml: 'auto' }} disabled={!session.report} onClick={() => openAltlocReview(id)}>Review states</Button>
    </Box>
    <InspectionTable aria-label="Alternate-state selection">
      <TableHead><TableRow>
        {['Keep', 'Model', 'Chain', 'Residue', 'State', 'Occupancy', 'Atoms', 'Notes'].map(label => <TableCell key={label} align={label === 'Occupancy' || label === 'Atoms' ? 'right' : 'left'}>{label}</TableCell>)}
      </TableRow></TableHead>
      <TableBody>
        {!session.report && <TableRow><TableCell colSpan={8}>Choose a structure and Analyze in DVBfixer → Altlocs.</TableCell></TableRow>}
        {session.report && !session.report.residues.length && <TableRow><TableCell colSpan={8}>No alternate locations.</TableCell></TableRow>}
        {session.report?.residues.flatMap(row => row.states.map(state => <InspectionRow key={`${row.key}:${state.label}`}
          selected={session.activeKey === row.key && session.focusLabel === state.label}
          onClick={() => { focus(id, row.key, state.label); window.dispatchEvent(new Event('dvbfixer:focus-altloc-viewer')) }}
          >
          <TableCell padding="checkbox"><Checkbox sx={compactIconButtonSx} checked={session.choices[row.key] === state.label} disabled={session.busy}
            inputProps={{ 'aria-label': `Keep model ${row.model} ${row.chain.trim() || '_'}:${row.residueNumber}${row.insertionCode} occurrence ${row.occurrence} state ${state.label}` }}
            onClick={event => event.stopPropagation()} onChange={() => {
              choose(id, row.key, state.label)
              window.dispatchEvent(new Event('dvbfixer:focus-altloc-viewer'))
            }} /></TableCell>
          <TableCell>{row.model}</TableCell><TableCell>{row.chain.trim() || '_'}</TableCell>
          <TableCell><Tooltip title={`Occurrence ${row.occurrence}`}><Box component="span"><ResidueLabel chain={row.chain} name={state.residueName} number={row.residueNumber} insertionCode={row.insertionCode} /></Box></Tooltip></TableCell>
          <TableCell>{state.label}</TableCell><TableCell align="right">{state.occupancy === null ? '—' : state.occupancy.toFixed(2)}</TableCell>
          <TableCell align="right">{state.atomCount}</TableCell>
          <TableCell sx={{ color: row.warnings.length ? 'warning.main' : 'text.secondary', maxWidth: 240 }}><Tooltip title={row.warnings.join(' ') || 'Shared atoms retained'}><Box sx={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{row.warnings.join(' ') || 'Shared atoms retained'}</Box></Tooltip></TableCell>
        </InspectionRow>))}
      </TableBody>
    </InspectionTable>
    {session.report && <Box sx={{ px: 1.5, py: 0.5, borderTop: 1, borderColor: 'divider' }}><Typography variant="caption" color="text.secondary">
      {session.report.residues.length} alternate residues · {session.report.residues.filter(row => !session.choices[row.key]).length} unresolved · {session.report.warnings.length} warnings
    </Typography></Box>}
  </Box>
}
