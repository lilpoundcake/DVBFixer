import { useCallback, useEffect, useRef, useState } from 'react'
import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import CircularProgress from '@mui/material/CircularProgress'
import IconButton from '@mui/material/IconButton'
import List from '@mui/material/List'
import ListItemButton from '@mui/material/ListItemButton'
import ListItemText from '@mui/material/ListItemText'
import Menu from '@mui/material/Menu'
import MenuItem from '@mui/material/MenuItem'
import TextField from '@mui/material/TextField'
import ToggleButton from '@mui/material/ToggleButton'
import ToggleButtonGroup from '@mui/material/ToggleButtonGroup'
import Tooltip from '@mui/material/Tooltip'
import Typography from '@mui/material/Typography'
import AddIcon from '@mui/icons-material/Add'
import FolderOpenIcon from '@mui/icons-material/FolderOpen'
import UploadFileIcon from '@mui/icons-material/UploadFile'
import ArrowUpwardIcon from '@mui/icons-material/ArrowUpward'
import ArrowDownwardIcon from '@mui/icons-material/ArrowDownward'
import { useStructureStore, type ViewerSlot } from '../stores/structureStore'
import { useWorkspaceStore, workspaceFileUrl, type WorkspaceArtifact, type WorkspaceManifest } from '../stores/workspaceStore'
import { useSelectionStore } from '../stores/selectionStore'

export function ProjectLibrary() {
  const uploadRef = useRef<HTMLInputElement>(null)
  const workspaces = useWorkspaceStore(state => state.workspaces)
  const active = useWorkspaceStore(state => state.active)
  const loading = useWorkspaceStore(state => state.loading)
  const error = useWorkspaceStore(state => state.error)
  const refresh = useWorkspaceStore(state => state.refresh)
  const activate = useWorkspaceStore(state => state.activate)
  const createWorkspace = useWorkspaceStore(state => state.createWorkspace)
  const updateWorkspace = useWorkspaceStore(state => state.update)
  const saveWorkspace = useWorkspaceStore(state => state.save)
  const reloadWorkspace = useWorkspaceStore(state => state.reload)
  const discardActive = useWorkspaceStore(state => state.discardActive)
  const setTextPreview = useWorkspaceStore(state => state.setTextPreview)
  const [selectedArtifacts, setSelectedArtifacts] = useState<Set<string>>(new Set())
  const selectionAnchor = useRef<number | null>(null)
  const [loadTarget, setLoadTarget] = useState<ViewerSlot>('primary')
  const [contextMenu, setContextMenu] = useState<{ mouseX: number; mouseY: number; kind: 'workspace' | 'artifact'; id: string; name: string } | null>(null)
  const plugin = useStructureStore(state => state.plugin)
  const secondaryPlugin = useStructureStore(state => state.secondaryPlugin)
  const setFileName = useStructureStore(state => state.setFileName)
  const setSecondaryFileName = useStructureStore(state => state.setSecondaryFileName)
  const setMeta = useStructureStore(state => state.setMeta)
  const setLoading = useStructureStore(state => state.setLoading)
  const setError = useStructureStore(state => state.setError)
  const clearSelection = useSelectionStore(state => state.clearSelection)

  const loadArtifact = useCallback(async (workspace: WorkspaceManifest, artifact: WorkspaceArtifact, slot: ViewerSlot) => {
    if (artifact.kind !== 'structure') {
      if (/\.(fasta|fa|faa|pir|aln|json|txt|html|csv|log|dat|top|itp|mdp|md|ya?ml|toml|xml|pml|py|sh|mol2|sdf)$/i.test(artifact.file)) {
        setTextPreview({ workspaceId: workspace.id, file: artifact.file, name: artifact.name })
        window.dispatchEvent(new Event('dvbfixer:open-text-viewer'))
      } else window.open(workspaceFileUrl(workspace.id, artifact.file), '_blank', 'noopener,noreferrer')
      return
    }
    const target = slot === 'primary' ? plugin : secondaryPlugin
    if (!target) { setError(slot === 'primary' ? 'Open a 3D Structure tab first' : 'Open a 3D Structure (B) tab first'); return }
    setLoading(true); setError(null)
    if (slot === 'primary') clearSelection()
    try {
      const response = await fetch(workspaceFileUrl(workspace.id, artifact.file))
      if (!response.ok) throw new Error(`HTTP ${response.status}`)
      await target.clear()
      const data = await target.builders.data.rawData({ data: await response.text(), label: artifact.file })
      const format = /\.(cif|mmcif)$/i.test(artifact.file) ? 'mmcif' : 'pdb'
      const trajectory = await target.builders.structure.parseTrajectory(data, format)
      await target.builders.structure.hierarchy.applyPreset(trajectory, 'default')
      if (slot === 'primary') {
        setFileName(artifact.file)
        setMeta({ name: artifact.name })
        updateWorkspace({ primaryFile: artifact.file })
      } else {
        setSecondaryFileName(artifact.file)
        updateWorkspace({ secondaryFile: artifact.file })
      }
    } catch (reason: any) {
      setError(`Failed to load ${artifact.name}: ${reason.message || String(reason)}`)
    } finally { setLoading(false) }
  }, [clearSelection, plugin, secondaryPlugin, setError, setFileName, setLoading, setMeta, setSecondaryFileName, setTextPreview, updateWorkspace])

  useEffect(() => { setSelectedArtifacts(new Set()); selectionAnchor.current = null }, [active?.id])

  const openWorkspace = useCallback(async (id: string) => {
    const workspace = await activate(id)
    const primary = workspace.artifacts.find(item => item.file === workspace.primaryFile)
    const secondary = workspace.artifacts.find(item => item.file === workspace.secondaryFile)
    if (primary) await loadArtifact(workspace, primary, 'primary')
    else { setFileName(null); if (plugin) await plugin.clear() }
    if (secondary && secondaryPlugin) await loadArtifact(workspace, secondary, 'secondary')
    else { setSecondaryFileName(null); if (secondaryPlugin) await secondaryPlugin.clear() }
  }, [activate, loadArtifact, plugin, secondaryPlugin, setFileName, setSecondaryFileName])

  useEffect(() => {
    refresh().then(() => {
      const current = useWorkspaceStore.getState()
      if (!current.active && current.workspaces[0]) openWorkspace(current.workspaces[0].id).catch(() => {})
      else if (!current.active && !current.workspaces.length) createWorkspace().catch(() => {})
    }).catch(() => {})
  }, [createWorkspace, openWorkspace, refresh])

  const importFiles = async (files: FileList | null) => {
    if (!active || !files?.length) return
    setLoading(true); setError(null)
    try {
      await saveWorkspace()
      for (const file of Array.from(files)) {
        const response = await fetch(`/api/workspaces/${encodeURIComponent(active.id)}/import`, {
          method: 'POST', headers: { 'X-File-Name': encodeURIComponent(file.name) }, body: file,
        })
        const body = await response.json().catch(() => ({}))
        if (!response.ok) throw new Error(body.error || `Import failed: ${file.name}`)
      }
      await reloadWorkspace()
      await refresh()
    } catch (reason: any) { setError(reason.message || String(reason)) }
    finally { setLoading(false); if (uploadRef.current) uploadRef.current.value = '' }
  }

  const removeContextItem = async (target: NonNullable<typeof contextMenu>) => {
    setContextMenu(null)
    setLoading(true); setError(null)
    try {
      if (target.kind === 'artifact') {
        if (!active) return
        await saveWorkspace()
        const targetIds = selectedArtifacts.has(target.id) ? [...selectedArtifacts] : [target.id]
        const artifacts = active.artifacts.filter(item => targetIds.includes(item.id))
        for (const artifact of artifacts) {
          const response = await fetch(`/api/workspaces/${encodeURIComponent(active.id)}/artifacts/${encodeURIComponent(artifact.id)}`, { method: 'DELETE' })
          if (!response.ok) {
            const body = await response.json().catch(() => ({}))
            throw new Error(body.error || `Delete failed for ${artifact.name}: HTTP ${response.status}`)
          }
        }
        if (artifacts.some(artifact => artifact.file === useStructureStore.getState().fileName)) {
          if (plugin) await plugin.clear()
          setFileName(null)
          clearSelection()
        }
        if (artifacts.some(artifact => artifact.file === useStructureStore.getState().secondaryFileName)) {
          if (secondaryPlugin) await secondaryPlugin.clear()
          setSecondaryFileName(null)
        }
        setSelectedArtifacts(new Set())
        selectionAnchor.current = null
        await reloadWorkspace()
        await refresh()
      } else {
        const wasActive = active?.id === target.id
        if (wasActive) await saveWorkspace()
        const response = await fetch(`/api/workspaces/${encodeURIComponent(target.id)}`, { method: 'DELETE' })
        if (!response.ok) {
          const body = await response.json().catch(() => ({}))
          throw new Error(body.error || `Delete failed: HTTP ${response.status}`)
        }
        if (wasActive) discardActive()
        await refresh()
        if (wasActive) {
          if (plugin) await plugin.clear()
          if (secondaryPlugin) await secondaryPlugin.clear()
          setFileName(null); setSecondaryFileName(null); clearSelection()
          const next = useWorkspaceStore.getState().workspaces[0]
          if (next) await openWorkspace(next.id)
          else await createWorkspace('Untitled workspace')
        }
      }
    } catch (reason: any) { setError(reason.message || String(reason)) }
    finally { setLoading(false) }
  }

  const renameContextItem = async (target: NonNullable<typeof contextMenu>) => {
    setContextMenu(null)
    const name = window.prompt(`Rename ${target.kind}`, target.name)?.trim()
    if (!name || name === target.name) return
    try {
      await saveWorkspace()
      const url = target.kind === 'workspace'
        ? `/api/workspaces/${encodeURIComponent(target.id)}/rename`
        : `/api/workspaces/${encodeURIComponent(active?.id || '')}/artifacts/${encodeURIComponent(target.id)}/rename`
      const response = await fetch(url, { method: 'PATCH', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ name }) })
      const body = await response.json().catch(() => ({}))
      if (!response.ok) throw new Error(body.error || `Rename failed: HTTP ${response.status}`)
      if (target.kind === 'artifact') await reloadWorkspace()
      else if (active?.id === target.id) await reloadWorkspace()
      await refresh()
    } catch (reason: any) { setError(reason.message || String(reason)) }
  }

  const reorderSelected = async (direction: -1 | 1) => {
    if (!active || !selectedArtifacts.size) return
    const visible = active.artifacts.filter(item => !item.hidden)
    if (direction < 0) {
      for (let index = 1; index < visible.length; index++) if (selectedArtifacts.has(visible[index].id) && !selectedArtifacts.has(visible[index - 1].id)) [visible[index - 1], visible[index]] = [visible[index], visible[index - 1]]
    } else {
      for (let index = visible.length - 2; index >= 0; index--) if (selectedArtifacts.has(visible[index].id) && !selectedArtifacts.has(visible[index + 1].id)) [visible[index], visible[index + 1]] = [visible[index + 1], visible[index]]
    }
    let visibleIndex = 0
    updateWorkspace({ artifacts: active.artifacts.map(item => item.hidden ? item : visible[visibleIndex++]) })
    await useWorkspaceStore.getState().save()
  }

  return <Box sx={{ height: '100%', display: 'flex', flexDirection: 'column' }}>
    <Box sx={{ px: 1, py: 0.75, borderBottom: 1, borderColor: 'divider', display: 'flex', gap: 0.5, alignItems: 'center' }}>
      <Typography variant="caption" sx={{ fontWeight: 700, textTransform: 'uppercase', letterSpacing: 1 }}>Projects</Typography>
      <Tooltip title="New workspace"><IconButton size="small" onClick={() => createWorkspace()}><AddIcon fontSize="small" /></IconButton></Tooltip>
      <input ref={uploadRef} hidden multiple type="file" onChange={event => importFiles(event.target.files)} />
      <Tooltip title="Import files into active workspace"><span><IconButton size="small" disabled={!active} onClick={() => uploadRef.current?.click()}><UploadFileIcon fontSize="small" /></IconButton></span></Tooltip>
      <ToggleButtonGroup size="small" exclusive value={loadTarget} onChange={(_event, value) => value && setLoadTarget(value)} sx={{ ml: 'auto' }}>
        <ToggleButton value="primary">A</ToggleButton><ToggleButton value="secondary" disabled={!secondaryPlugin}>B</ToggleButton>
      </ToggleButtonGroup>
    </Box>
    {error && <Typography variant="caption" color="error" sx={{ px: 1, py: 0.5 }}>{error}</Typography>}
    {loading && <CircularProgress size={16} sx={{ m: 1 }} />}
    <List dense disablePadding sx={{ borderBottom: 1, borderColor: 'divider' }}>
      {workspaces.map(workspace => <ListItemButton key={workspace.id} selected={active?.id === workspace.id} onClick={() => openWorkspace(workspace.id)}
        onContextMenu={event => { event.preventDefault(); setContextMenu({ mouseX: event.clientX + 2, mouseY: event.clientY - 6, kind: 'workspace', id: workspace.id, name: workspace.name }) }}>
        <FolderOpenIcon sx={{ fontSize: 16, mr: 1 }} />
        <ListItemText primary={workspace.name} secondary={`${workspace.artifactCount} files`} />
      </ListItemButton>)}
    </List>
    {active && <>
      <Box sx={{ p: 1 }}><TextField size="small" fullWidth label="Workspace name" value={active.name} onChange={event => updateWorkspace({ name: event.target.value })} /></Box>
      <Box sx={{ px: 1, display: 'flex', alignItems: 'center' }}><Typography variant="caption" sx={{ fontWeight: 700, color: 'text.secondary' }}>FILES</Typography>
        <Tooltip title="Move selected files up"><span><IconButton size="small" disabled={!selectedArtifacts.size} onClick={() => reorderSelected(-1)}><ArrowUpwardIcon fontSize="inherit" /></IconButton></span></Tooltip>
        <Tooltip title="Move selected files down"><span><IconButton size="small" disabled={!selectedArtifacts.size} onClick={() => reorderSelected(1)}><ArrowDownwardIcon fontSize="inherit" /></IconButton></span></Tooltip>
        {!!selectedArtifacts.size && <Typography variant="caption" color="primary">{selectedArtifacts.size} selected</Typography>}
      </Box>
      <List dense disablePadding sx={{ flex: 1, overflow: 'auto' }}>
        {active.artifacts.filter(artifact => !artifact.hidden).map((artifact, index, visible) => <ListItemButton key={artifact.id} selected={selectedArtifacts.has(artifact.id)} onClick={event => {
          if (event.shiftKey && selectionAnchor.current !== null) {
            const start = Math.min(selectionAnchor.current, index); const end = Math.max(selectionAnchor.current, index)
            setSelectedArtifacts(current => new Set([...current, ...visible.slice(start, end + 1).map(item => item.id)]))
          } else if (event.ctrlKey || event.metaKey) {
            setSelectedArtifacts(current => { const next = new Set(current); if (next.has(artifact.id)) next.delete(artifact.id); else next.add(artifact.id); return next })
            selectionAnchor.current = index
          } else {
            setSelectedArtifacts(new Set([artifact.id])); selectionAnchor.current = index
            loadArtifact(active, artifact, loadTarget).catch(() => {})
          }
        }} sx={{
          pl: 2,
          bgcolor: selectedArtifacts.has(artifact.id) ? 'rgba(29,50,97,0.24)' : index % 2 ? 'rgba(29,50,97,0.045)' : 'transparent',
          '&.Mui-selected': { bgcolor: 'rgba(29,50,97,0.24)' },
          '&.Mui-selected:hover': { bgcolor: 'rgba(29,50,97,0.31)' },
        }}
          onContextMenu={event => { event.preventDefault(); setContextMenu({ mouseX: event.clientX + 2, mouseY: event.clientY - 6, kind: 'artifact', id: artifact.id, name: artifact.name }) }}>
          <ListItemText primary={artifact.name} secondary={[artifact.folder, artifact.file].filter(Boolean).join(' · ')} />
        </ListItemButton>)}
        {!active.artifacts.length && <Box sx={{ p: 2, textAlign: 'center' }}><Typography variant="caption" color="text.secondary">Import files to begin.</Typography></Box>}
      </List>
    </>}
    <Button size="small" onClick={() => refresh()} sx={{ m: 0.5 }}>Refresh</Button>
    <Menu open={!!contextMenu} onClose={() => setContextMenu(null)} anchorReference="anchorPosition"
      anchorPosition={contextMenu ? { top: contextMenu.mouseY, left: contextMenu.mouseX } : undefined}>
      <MenuItem onClick={() => contextMenu && renameContextItem(contextMenu)}>Rename</MenuItem>
      <MenuItem onClick={() => contextMenu && removeContextItem(contextMenu)}>Move {contextMenu?.kind === 'workspace' ? 'workspace' : contextMenu && selectedArtifacts.has(contextMenu.id) && selectedArtifacts.size > 1 ? `${selectedArtifacts.size} selected files` : 'file'} to trash</MenuItem>
    </Menu>
  </Box>
}
