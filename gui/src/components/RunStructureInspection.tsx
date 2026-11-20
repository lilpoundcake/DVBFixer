import { useEffect, useRef, useState } from 'react'
import Alert from '@mui/material/Alert'
import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import Dialog from '@mui/material/Dialog'
import DialogTitle from '@mui/material/DialogTitle'
import DialogContent from '@mui/material/DialogContent'
import DialogActions from '@mui/material/DialogActions'
import MenuItem from '@mui/material/MenuItem'
import TextField from '@mui/material/TextField'
import Typography from '@mui/material/Typography'
import { apiFetch } from '../lib/api-client'
import { panelToolbarSx } from '../ui/layout'

/** Isolated preview: never writes a workspace artifact or replaces either viewer. */
function DiagnosticViewport({ url, label }: { url: string; label: string }) {
  const target = useRef<HTMLDivElement>(null)
  const [error, setError] = useState('')
  const [loading, setLoading] = useState(true)
  useEffect(() => {
    const controller = new AbortController()
    let plugin: import('molstar/lib/mol-plugin-ui/context').PluginUIContext | undefined
    const load = async () => {
      const response = await apiFetch(url, { signal: controller.signal, cache: 'no-store' })
      if (!response.ok) throw new Error(`Cannot preview structure: HTTP ${response.status}`)
      const pdb = await response.text()
      const [{ createPluginUI }, { renderReact18 }, { DefaultPluginUISpec }] = await Promise.all([
        import('molstar/lib/mol-plugin-ui'), import('molstar/lib/mol-plugin-ui/react18'), import('molstar/lib/mol-plugin-ui/spec'),
      ])
      if (controller.signal.aborted || !target.current) return
      plugin = await createPluginUI({ target: target.current, render: renderReact18, spec: {
        ...DefaultPluginUISpec(), layout: { initial: { isExpanded: false, showControls: false } },
      } })
      if (controller.signal.aborted) { plugin.dispose(); plugin = undefined; return }
      const data = await plugin.builders.data.rawData({ data: pdb, label })
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb')
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default')
      if (!controller.signal.aborted) setLoading(false)
    }
    void load().catch(reason => { if (!controller.signal.aborted) setError(String(reason)) })
    return () => { controller.abort(); plugin?.dispose() }
  }, [url, label])
  return <Box sx={{ position: 'relative', height: '55vh', minHeight: 260 }}>
    <div ref={target} style={{ position: 'absolute', inset: 0 }} />
    {loading && !error && <Typography variant="body2" sx={{ position: 'relative', p: 1 }}>Loading diagnostic structure…</Typography>}
    {error && <Alert severity="error" sx={{ position: 'relative' }}>{error}</Alert>}
  </Box>
}

export function RunStructureInspection({ workspaceId, runDirectory }: { workspaceId: string; runDirectory: string }) {
  const [open, setOpen] = useState(false)
  const [files, setFiles] = useState<string[]>([])
  const [selected, setSelected] = useState('')
  const [error, setError] = useState('')
  const [loading, setLoading] = useState(false)
  const endpoint = `/api/workspaces/${encodeURIComponent(workspaceId)}/run-inspection?summary=${encodeURIComponent(`${runDirectory}/run-summary.json`)}`
  useEffect(() => {
    if (!open) return
    const controller = new AbortController()
    apiFetch(endpoint, { signal: controller.signal }).then(async response => {
      if (!response.ok) throw new Error(`Cannot list retained structures: HTTP ${response.status}`)
      const data = await response.json() as { files: string[] }
      if (controller.signal.aborted) return
      setFiles(data.files); setSelected(data.files[0] || '')
    }).catch(reason => { if (!controller.signal.aborted) setError(String(reason)) })
      .finally(() => { if (!controller.signal.aborted) setLoading(false) })
    return () => controller.abort()
  }, [open, endpoint])
  return <>
    <Button onClick={() => { setError(''); setFiles([]); setSelected(''); setLoading(true); setOpen(true) }}>Inspect retained structures</Button>
    <Dialog open={open} onClose={() => setOpen(false)} fullWidth maxWidth="lg">
      <DialogTitle>Diagnostic structures · inspection only</DialogTitle>
      <DialogContent>
        <Alert severity="warning">Retained source, intermediate, or rejected structures—not accepted run outputs. Viewing does not change their status or the active structure. See the run logs and per-seed report for validation findings.</Alert>
        {error && <Alert severity="error">{error}</Alert>}
        <Box sx={panelToolbarSx}>
          {loading ? 'Loading retained structures…' : !files.length ? 'No retained diagnostic PDB structures available.' :
            <TextField select label={`Structure (${files.length})`} value={selected} fullWidth onChange={event => setSelected(event.target.value)}>
              {files.map(file => <MenuItem key={file} value={file}>{file}</MenuItem>)}
            </TextField>}
        </Box>
        {selected && <DiagnosticViewport key={selected} url={`${endpoint}&file=${encodeURIComponent(selected)}`} label={`Diagnostic only · ${selected}`} />}
      </DialogContent>
      <DialogActions><Button onClick={() => setOpen(false)}>Close</Button></DialogActions>
    </Dialog>
  </>
}
