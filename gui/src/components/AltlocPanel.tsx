import { useEffect } from 'react'
import { Alert, Box, Button, FormControl, InputLabel, MenuItem, Select, TextField, Typography } from '@mui/material'
import { useWorkspaceStore } from '../stores/workspaceStore'
import { useAltlocStore, EMPTY_ALTLOC_SESSION } from '../stores/altlocStore'
import type { AltlocPolicy, AltlocRequest } from '../lib/altloc-types'
import { actionGroupSx, formActionRowSx, primaryActionSx, controlSx, wideControlSx } from '../ui/layout'
import { altlocSelectionComplete } from '../lib/altloc-selection'
import { runAltlocOperation } from '../lib/altloc-operations'
import { useStructureStore } from '../stores/structureStore'
import { openedAltlocArtifact } from '../lib/altloc-input'

/** Routine CLI launcher only. The dockable Altlocs table owns manual editing. */
export function AltlocPanel() {
  const workspace = useWorkspaceStore(state => state.active)
  const id = workspace?.id || ''
  const session = useAltlocStore(state => state.sessions[id] || EMPTY_ALTLOC_SESSION)
  const configure = useAltlocStore(state => state.configure)
  const fileName = useStructureStore(state => state.fileName)
  const secondaryFileName = useStructureStore(state => state.secondaryFileName)
  const loadTargetSlot = useStructureStore(state => state.loadTargetSlot)
  const { input, inputMode, policy, prefer, fallback, report, choices, busy, error, result, reportFile } = session
  const artifacts = workspace?.artifacts.filter(item => item.kind === 'structure' && !item.hidden && /\.(pdb|ent|cif|mmcif)$/i.test(item.file)) || []
  const validInput = artifacts.some(item => item.id === input)
  const openedFile = loadTargetSlot === 'secondary' ? secondaryFileName : fileName
  useEffect(() => {
    if (!id || busy || (validInput && inputMode === 'manual') || !openedFile) return
    const opened = openedAltlocArtifact(workspace?.artifacts || [], openedFile, input, inputMode)
    if (opened) configure(id, { input: opened, inputMode: 'auto' })
  }, [id, busy, validInput, inputMode, input, openedFile, workspace, configure])
  const resolved = !!report && altlocSelectionComplete(report.residues, choices)
  return <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1.5 }}>
    <Typography variant="caption" color="text.secondary">Select one alternate state per residue without changing coordinates or occupancy. Manual editing opens the separate Altlocs table; visualization uses 3D Structure.</Typography>
    <Box sx={formActionRowSx}>
      <FormControl size="small" sx={wideControlSx} disabled={busy}>
        <InputLabel>Input structure</InputLabel>
        <Select label="Input structure" value={validInput ? input : ''} onChange={event => configure(id, { input: event.target.value, inputMode: 'manual' })}>
          {artifacts.map(item => <MenuItem key={item.id} value={item.id}>{item.name}</MenuItem>)}
        </Select>
      </FormControl>
      <Box sx={actionGroupSx}>
        <Button variant="outlined" sx={primaryActionSx} disabled={busy || !validInput || (policy === 'prefer' && !prefer.trim())} onClick={() => void runAltlocOperation(id, 'analyze')}>Analyze</Button>
        <Button variant="contained" sx={primaryActionSx} disabled={busy || !validInput || (policy === 'manual' && !resolved) || (policy === 'prefer' && !prefer.trim())}
          onClick={() => void runAltlocOperation(id, policy === 'manual' || resolved ? 'apply' : 'run')}>
          {busy ? 'Working…' : policy === 'manual' || resolved ? 'Apply selection' : 'Run Altlocs'}
        </Button>
      </Box>
    </Box>
    <Typography variant="subtitle2">Selection</Typography>
    <Box sx={{ ...formActionRowSx, flexWrap: 'wrap', gap: 1 }}>
      <FormControl size="small" sx={controlSx} disabled={busy}><InputLabel>Policy</InputLabel>
        <Select label="Policy" value={policy} onChange={event => configure(id, { policy: event.target.value as AltlocPolicy })}>
          <MenuItem value="manual">Manual</MenuItem><MenuItem value="occupancy">Occupancy</MenuItem>
          <MenuItem value="first">First encountered</MenuItem><MenuItem value="prefer">Prefer label</MenuItem>
        </Select>
      </FormControl>
      {policy === 'prefer' && <>
        <TextField size="small" label="Preferred label" sx={controlSx} disabled={busy} value={prefer} slotProps={{ htmlInput: { maxLength: 1 } }} onChange={event => configure(id, { prefer: event.target.value })} />
        <FormControl size="small" sx={controlSx} disabled={busy}><InputLabel>Missing label</InputLabel><Select label="Missing label" value={fallback} onChange={event => configure(id, { fallback: event.target.value as AltlocRequest['fallback'] })}>
          <MenuItem value="error">Manual choice</MenuItem><MenuItem value="occupancy">Use occupancy</MenuItem><MenuItem value="first">Use first</MenuItem>
        </Select></FormControl>
      </>}
    </Box>
    {error && <Alert severity="error">{error}</Alert>}
    {result && <Alert severity="success">{result}</Alert>}
    {report && <>
      {Object.entries(report.chainMap).some(([from, to]) => from !== to) && <Alert severity="info">CIF chain mapping: {JSON.stringify(report.chainMap)}</Alert>}
      {report.warnings.filter(warning => !report.residues.some(item => warning.startsWith(item.key))).map(warning => <Alert key={warning} severity="warning">{warning}</Alert>)}
    </>}
    {reportFile && <Button sx={{ alignSelf: 'flex-start' }} onClick={() => {
      useWorkspaceStore.getState().setTextPreview({ workspaceId: id, file: reportFile, name: 'Altloc selection report' })
      window.dispatchEvent(new Event('dvbfixer:open-text-viewer'))
    }}>Open selection report</Button>}
  </Box>
}
