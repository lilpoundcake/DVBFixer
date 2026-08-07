import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import Alert from '@mui/material/Alert'
import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import Chip from '@mui/material/Chip'
import CircularProgress from '@mui/material/CircularProgress'
import FormControl from '@mui/material/FormControl'
import IconButton from '@mui/material/IconButton'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import Select from '@mui/material/Select'
import Tab from '@mui/material/Tab'
import Tabs from '@mui/material/Tabs'
import TextField from '@mui/material/TextField'
import Tooltip from '@mui/material/Tooltip'
import Typography from '@mui/material/Typography'
import AddIcon from '@mui/icons-material/Add'
import DeleteIcon from '@mui/icons-material/DeleteOutlined'
import RedoIcon from '@mui/icons-material/Redo'
import UndoIcon from '@mui/icons-material/Undo'
import PlayArrowIcon from '@mui/icons-material/PlayArrow'
import AccountTreeIcon from '@mui/icons-material/AccountTree'
import { useStructureStore } from '../stores/structureStore'

interface TemplateSelection {
  id: string
  file: string
  chain: string
  targetChain: string
  label?: string
}
interface AlignmentSpan { start: number; end: number }
interface AlignmentRow {
  id: string
  kind: 'target' | 'template'
  sequence: string
  templateId?: string
}
interface AlignmentGroup {
  chainId: string
  rows: AlignmentRow[]
  masks: Record<string, AlignmentSpan[]>
}
interface HomologyProject {
  version: 1
  id: string
  name: string
  targetFasta: string
  templates: TemplateSelection[]
  engine: 'mafft' | 'muscle' | 'clustalo'
  alignmentGroups: AlignmentGroup[]
  structuralAlignment?: string
  modelOptions: Record<string, string | number | boolean>
  createdAt: string
  updatedAt: string
}
interface ProjectSummary { id: string; name: string; updatedAt: string }
interface ArtifactEntry { file?: string; name?: string; kind?: string }
interface EngineStatus { available: boolean; path: string | null }

function targetChainIds(fasta: string): string[] {
  return fasta.split(/\r?\n/)
    .filter(line => line.trim().startsWith('>'))
    .map(line => line.trim().slice(1).split(/\s+/)[0])
    .filter(Boolean)
}

function cloneGroups(groups: AlignmentGroup[]): AlignmentGroup[] {
  return JSON.parse(JSON.stringify(groups)) as AlignmentGroup[]
}

export function HomologyPanel() {
  const [project, setProject] = useState<HomologyProject | null>(null)
  const [projects, setProjects] = useState<ProjectSummary[]>([])
  const [artifacts, setArtifacts] = useState<ArtifactEntry[]>([])
  const [engines, setEngines] = useState<Record<string, EngineStatus>>({})
  const [tab, setTab] = useState(0)
  const [busy, setBusy] = useState('')
  const [error, setError] = useState('')
  const [message, setMessage] = useState('')
  const [saved, setSaved] = useState(true)
  const [selection, setSelection] = useState<Record<string, { start: number; end: number }>>({})
  const [undo, setUndo] = useState<AlignmentGroup[][]>([])
  const [redo, setRedo] = useState<AlignmentGroup[][]>([])
  const alignmentImportRef = useRef<HTMLInputElement>(null)
  const libraryVersion = useStructureStore(state => state.libraryVersion)

  const refreshProjects = useCallback(async () => {
    const response = await fetch('/api/homology/projects', { cache: 'no-store' })
    if (!response.ok) throw new Error(`Projects: HTTP ${response.status}`)
    setProjects(await response.json() as ProjectSummary[])
  }, [])

  const loadProject = useCallback(async (id: string) => {
    setBusy('Loading project')
    try {
      const response = await fetch(`/api/homology/projects/${encodeURIComponent(id)}`, { cache: 'no-store' })
      const body = await response.json() as HomologyProject & { error?: string }
      if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
      setProject(body)
      setUndo([])
      setRedo([])
      setSaved(true)
    } finally {
      setBusy('')
    }
  }, [])

  const createProject = useCallback(async () => {
    setBusy('Creating project')
    setError('')
    try {
      const response = await fetch('/api/homology/projects', { method: 'POST' })
      const body = await response.json() as HomologyProject & { error?: string }
      if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
      setProject(body)
      await refreshProjects()
      setSaved(true)
    } catch (reason: any) {
      setError(reason.message || String(reason))
    } finally {
      setBusy('')
    }
  }, [refreshProjects])

  useEffect(() => {
    refreshProjects()
      .then(async () => {
        const response = await fetch('/api/homology/projects', { cache: 'no-store' })
        const list = await response.json() as ProjectSummary[]
        if (list[0]) await loadProject(list[0].id)
        else await createProject()
      })
      .catch(reason => setError(reason.message || String(reason)))
  }, [createProject, loadProject, refreshProjects])

  useEffect(() => {
    fetch(`/structures/index.json?t=${Date.now()}`, { cache: 'no-store' })
      .then(response => response.ok ? response.json() : [])
      .then((entries: ArtifactEntry[]) => setArtifacts(entries.filter(entry => entry.file && /\.(pdb|cif|mmcif)$/i.test(entry.file))))
      .catch(() => setArtifacts([]))
  }, [libraryVersion])

  useEffect(() => {
    fetch('/api/homology/engines', { cache: 'no-store' })
      .then(response => response.ok ? response.json() : {})
      .then((status: Record<string, EngineStatus>) => setEngines(status))
      .catch(() => setEngines({}))
  }, [])

  useEffect(() => {
    if (!project || saved) return
    const timer = window.setTimeout(async () => {
      try {
        const response = await fetch(`/api/homology/projects/${encodeURIComponent(project.id)}`, {
          method: 'PUT', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(project),
        })
        if (!response.ok) throw new Error(`Autosave: HTTP ${response.status}`)
        setSaved(true)
        refreshProjects().catch(() => {})
      } catch (reason: any) {
        setError(reason.message || String(reason))
      }
    }, 700)
    return () => window.clearTimeout(timer)
  }, [project, refreshProjects, saved])

  const update = useCallback((patch: Partial<HomologyProject>) => {
    setProject(current => current ? { ...current, ...patch } : current)
    setSaved(false)
  }, [])

  const targetChains = useMemo(() => targetChainIds(project?.targetFasta || ''), [project?.targetFasta])

  const addTemplate = () => {
    if (!project || !artifacts[0]?.file) return
    const template: TemplateSelection = {
      id: crypto.randomUUID(), file: artifacts[0].file, chain: 'A', targetChain: targetChains[0] || 'A',
    }
    update({ templates: [...project.templates, template] })
  }

  const updateTemplate = (id: string, patch: Partial<TemplateSelection>) => {
    if (!project) return
    update({ templates: project.templates.map(template => template.id === id ? { ...template, ...patch } : template) })
  }

  const callProjectAction = async (action: 'align' | 'salign') => {
    if (!project) return
    setBusy(action === 'align' ? 'Generating alignment' : 'Structural alignment')
    setError('')
    setMessage('')
    try {
      const response = await fetch(`/api/homology/${action}`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(project),
      })
      const body = await response.json() as HomologyProject & { error?: string }
      if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
      setProject(body)
      setSaved(true)
      if (action === 'align') { setTab(2); setUndo([]); setRedo([]) }
      setMessage(action === 'align' ? 'Alignment generated' : `Structural alignment saved: ${body.structuralAlignment}`)
    } catch (reason: any) {
      setError(reason.message || String(reason))
    } finally {
      setBusy('')
    }
  }

  const commitGroups = (next: AlignmentGroup[]) => {
    if (!project) return
    setUndo(history => [...history.slice(-49), cloneGroups(project.alignmentGroups)])
    setRedo([])
    update({ alignmentGroups: next })
  }

  const editRow = (groupIndex: number, rowId: string, sequence: string) => {
    if (!project) return
    const original = project.alignmentGroups[groupIndex].rows.find(row => row.id === rowId)
    if (!original || sequence.replace(/-/g, '') !== original.sequence.replace(/-/g, '')) {
      setError('Alignment editing may move gaps but cannot change the underlying amino-acid sequence. Edit residues in Target instead.')
      return
    }
    setError('')
    const next = cloneGroups(project.alignmentGroups)
    next[groupIndex].rows = next[groupIndex].rows.map(row => row.id === rowId ? { ...row, sequence: sequence.toUpperCase() } : row)
    if (next[groupIndex].masks[rowId]) next[groupIndex].masks[rowId] = []
    commitGroups(next)
  }

  const insertGap = (groupIndex: number, row: AlignmentRow) => {
    const point = selection[row.id]?.start ?? row.sequence.length
    editRow(groupIndex, row.id, `${row.sequence.slice(0, point)}-${row.sequence.slice(point)}`)
  }

  const deleteGap = (groupIndex: number, row: AlignmentRow) => {
    const point = selection[row.id]?.start ?? 0
    const index = row.sequence[point] === '-' ? point : row.sequence.lastIndexOf('-', point)
    if (index >= 0) editRow(groupIndex, row.id, `${row.sequence.slice(0, index)}${row.sequence.slice(index + 1)}`)
  }

  const addMask = (groupIndex: number, row: AlignmentRow) => {
    if (!project || row.kind !== 'template') return
    const selected = selection[row.id] || { start: 0, end: row.sequence.length }
    if (selected.end <= selected.start) return
    const next = cloneGroups(project.alignmentGroups)
    const current = next[groupIndex].masks[row.id] || []
    next[groupIndex].masks[row.id] = [...current, selected].sort((a, b) => a.start - b.start)
    commitGroups(next)
  }

  const removeMask = (groupIndex: number, rowId: string, maskIndex: number) => {
    if (!project) return
    const next = cloneGroups(project.alignmentGroups)
    next[groupIndex].masks[rowId] = (next[groupIndex].masks[rowId] || []).filter((_mask, index) => index !== maskIndex)
    commitGroups(next)
  }

  const undoAlignment = () => {
    if (!project || !undo.length) return
    const previous = undo[undo.length - 1]
    setUndo(undo.slice(0, -1))
    setRedo(history => [...history, cloneGroups(project.alignmentGroups)])
    update({ alignmentGroups: previous })
  }

  const redoAlignment = () => {
    if (!project || !redo.length) return
    const next = redo[redo.length - 1]
    setRedo(redo.slice(0, -1))
    setUndo(history => [...history, cloneGroups(project.alignmentGroups)])
    update({ alignmentGroups: next })
  }

  const exportAlignment = () => {
    if (!project) return
    const text = project.alignmentGroups.flatMap(group => group.rows.map(row => (
      `>${group.chainId}|${row.id}\n${row.sequence}\n`
    ))).join('')
    const url = URL.createObjectURL(new Blob([text], { type: 'text/plain' }))
    const anchor = document.createElement('a')
    anchor.href = url
    anchor.download = `${project.name.replace(/[^a-zA-Z0-9_-]+/g, '_') || 'homology'}_alignment.fasta`
    anchor.click()
    URL.revokeObjectURL(url)
  }

  const importAlignment = async (file: File | undefined) => {
    if (!project || !file) return
    try {
      const imported = new Map<string, string>()
      let key = ''
      for (const raw of (await file.text()).split(/\r?\n/)) {
        const line = raw.trim()
        if (!line) continue
        if (line.startsWith('>')) {
          key = line.slice(1).split(/\s+/)[0]
          imported.set(key, '')
        } else if (key) imported.set(key, `${imported.get(key) || ''}${line.toUpperCase()}`)
      }
      const next = cloneGroups(project.alignmentGroups)
      for (const group of next) {
        for (const row of group.rows) {
          const sequence = imported.get(`${group.chainId}|${row.id}`)
          if (!sequence) throw new Error(`Imported alignment is missing ${group.chainId}|${row.id}`)
          if (sequence.replace(/-/g, '') !== row.sequence.replace(/-/g, '')) {
            throw new Error(`Imported alignment changes the ungapped sequence of ${row.id}`)
          }
          row.sequence = sequence
        }
        if (new Set(group.rows.map(row => row.sequence.length)).size !== 1) {
          throw new Error(`Imported rows for target chain ${group.chainId} have different lengths`)
        }
        group.masks = {}
      }
      commitGroups(next)
    } catch (reason: any) {
      setError(reason.message || String(reason))
    } finally {
      if (alignmentImportRef.current) alignmentImportRef.current.value = ''
    }
  }

  const runModel = async () => {
    if (!project) return
    setBusy('Running Modeller')
    setError('')
    setMessage('')
    try {
      const response = await fetch('/api/homology/model', {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(project),
      })
      const body = await response.json() as { error?: string; outputFile?: string; stderr?: string }
      if (!response.ok) throw new Error(body.error || body.stderr || `HTTP ${response.status}`)
      setMessage(`Model complete: ${body.outputFile}`)
      useStructureStore.getState().bumpLibraryVersion()
    } catch (reason: any) {
      setError(reason.message || String(reason))
    } finally {
      setBusy('')
    }
  }

  if (!project) return <Box sx={{ p: 2 }}><CircularProgress size={20} /> {busy || 'Loading homology workspace…'}</Box>

  return (
    <Box sx={{ height: '100%', display: 'flex', flexDirection: 'column' }}>
      <Box sx={{ px: 1, py: 0.75, display: 'flex', gap: 1, alignItems: 'center', borderBottom: 1, borderColor: 'divider' }}>
        <FormControl size="small" sx={{ minWidth: 220 }}>
          <InputLabel>Project</InputLabel>
          <Select label="Project" value={project.id} onChange={event => loadProject(event.target.value)}>
            {projects.map(item => <MenuItem key={item.id} value={item.id}>{item.name}</MenuItem>)}
          </Select>
        </FormControl>
        <TextField size="small" label="Project name" value={project.name} onChange={event => update({ name: event.target.value })} />
        <Tooltip title="New project"><IconButton onClick={createProject}><AddIcon /></IconButton></Tooltip>
        <Typography variant="caption" color="text.secondary" sx={{ ml: 'auto' }}>{saved ? 'Saved' : 'Saving…'}</Typography>
      </Box>
      <Tabs value={tab} onChange={(_event, value) => setTab(value)} sx={{ borderBottom: 1, borderColor: 'divider' }}>
        <Tab label="1 Target" /><Tab label="2 Templates" /><Tab label="3 Alignment" /><Tab label="4 Model" />
      </Tabs>
      <Box sx={{ flex: 1, overflow: 'auto', p: 1.5 }}>
        {error && <Alert severity="error" sx={{ mb: 1 }}>{error}</Alert>}
        {message && <Alert severity="success" sx={{ mb: 1 }}>{message}</Alert>}
        {busy && <Alert severity="info" icon={<CircularProgress size={16} />} sx={{ mb: 1 }}>{busy}</Alert>}

        {tab === 0 && (
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
            <Typography variant="body2">Enter one FASTA record per target chain. Sequence residues are edited here; the alignment tab edits gaps only.</Typography>
            <TextField multiline minRows={12} label="Target FASTA" value={project.targetFasta}
              onChange={event => update({ targetFasta: event.target.value, alignmentGroups: [] })}
              slotProps={{ htmlInput: { spellCheck: false } }}
              sx={{ '& textarea': { fontFamily: 'ui-monospace, monospace' } }} />
            <Typography variant="caption">Detected chains: {targetChains.join(', ') || 'none'}</Typography>
          </Box>
        )}

        {tab === 1 && (
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
            <Box sx={{ display: 'flex', gap: 1 }}>
              <Button startIcon={<AddIcon />} variant="outlined" onClick={addTemplate} disabled={!artifacts.length}>Add template</Button>
              <Button startIcon={<AccountTreeIcon />} variant="outlined" onClick={() => callProjectAction('salign')} disabled={busy !== '' || project.templates.length < 2}>Structural align</Button>
              <Button variant="contained" onClick={() => callProjectAction('align')} disabled={busy !== '' || !targetChains.length || !project.templates.length}>Generate MSA</Button>
              <FormControl size="small" sx={{ minWidth: 150 }}>
                <InputLabel>MSA engine</InputLabel>
                <Select label="MSA engine" value={project.engine} onChange={event => update({ engine: event.target.value as HomologyProject['engine'] })}>
                  <MenuItem value="mafft" disabled={engines.mafft?.available === false}>MAFFT{engines.mafft?.available === false ? ' (missing)' : ''}</MenuItem>
                  <MenuItem value="muscle" disabled={engines.muscle?.available === false}>MUSCLE 5{engines.muscle?.available === false ? ' (missing)' : ''}</MenuItem>
                  <MenuItem value="clustalo" disabled={engines.clustalo?.available === false}>Clustal Omega{engines.clustalo?.available === false ? ' (missing)' : ''}</MenuItem>
                </Select>
              </FormControl>
            </Box>
            {project.templates.map(template => (
              <Box key={template.id} sx={{ display: 'grid', gridTemplateColumns: 'minmax(280px, 1fr) 100px 160px 40px', gap: 1, alignItems: 'center' }}>
                <FormControl size="small"><InputLabel>Structure</InputLabel><Select label="Structure" value={template.file} onChange={event => updateTemplate(template.id, { file: event.target.value })}>
                  {artifacts.map(item => <MenuItem key={item.file} value={item.file}>{item.name || item.file}</MenuItem>)}
                </Select></FormControl>
                <TextField size="small" label="Chain" value={template.chain} onChange={event => updateTemplate(template.id, { chain: event.target.value })} />
                <FormControl size="small"><InputLabel>Target chain</InputLabel><Select label="Target chain" value={template.targetChain} onChange={event => updateTemplate(template.id, { targetChain: event.target.value })}>
                  {targetChains.map(chain => <MenuItem key={chain} value={chain}>{chain}</MenuItem>)}
                </Select></FormControl>
                <IconButton onClick={() => update({ templates: project.templates.filter(item => item.id !== template.id) })}><DeleteIcon /></IconButton>
              </Box>
            ))}
          </Box>
        )}

        {tab === 2 && (
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1.5 }}>
            <Box sx={{ display: 'flex', gap: 0.5 }}>
              <input ref={alignmentImportRef} hidden type="file" accept=".fasta,.fa,.faa,.aln" onChange={event => importAlignment(event.target.files?.[0])} />
              <Tooltip title="Undo alignment edit"><span><IconButton disabled={!undo.length} onClick={undoAlignment}><UndoIcon /></IconButton></span></Tooltip>
              <Tooltip title="Redo alignment edit"><span><IconButton disabled={!redo.length} onClick={redoAlignment}><RedoIcon /></IconButton></span></Tooltip>
              <Button onClick={() => callProjectAction('align')} disabled={busy !== ''}>Reset from {project.engine}</Button>
              <Button onClick={() => alignmentImportRef.current?.click()}>Import aligned FASTA</Button>
              <Button onClick={exportAlignment} disabled={!project.alignmentGroups.length}>Export alignment</Button>
            </Box>
            {!project.alignmentGroups.length && <Alert severity="info">Generate an alignment from the Templates tab.</Alert>}
            {project.alignmentGroups.map((group, groupIndex) => {
              const lengths = new Set(group.rows.map(row => row.sequence.length))
              return (
                <Box key={group.chainId} sx={{ border: 1, borderColor: lengths.size === 1 ? 'divider' : 'error.main', borderRadius: 1, p: 1 }}>
                  <Typography variant="subtitle2" sx={{ mb: 1 }}>Target chain {group.chainId} · {group.rows[0]?.sequence.length || 0} columns</Typography>
                  {lengths.size !== 1 && <Alert severity="warning" sx={{ mb: 1 }}>Rows must have equal lengths before modeling.</Alert>}
                  {group.rows.map(row => (
                    <Box key={row.id} sx={{ mb: 1 }}>
                      <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5, mb: 0.25 }}>
                        <Chip size="small" color={row.kind === 'target' ? 'primary' : 'default'} label={`${row.kind}: ${row.id}`} />
                        <Button size="small" onClick={() => insertGap(groupIndex, row)}>Insert gap</Button>
                        <Button size="small" onClick={() => deleteGap(groupIndex, row)}>Delete gap</Button>
                        {row.kind === 'template' && <Button size="small" onClick={() => addMask(groupIndex, row)}>Use selection as template span</Button>}
                      </Box>
                      <TextField fullWidth multiline maxRows={4} value={row.sequence}
                        onChange={event => editRow(groupIndex, row.id, event.target.value.replace(/\s+/g, ''))}
                        onSelect={event => {
                          const target = event.target as HTMLTextAreaElement
                          setSelection(current => ({ ...current, [row.id]: { start: target.selectionStart, end: target.selectionEnd } }))
                        }}
                        slotProps={{ htmlInput: { spellCheck: false } }}
                        sx={{ '& textarea': { fontFamily: 'ui-monospace, monospace', wordBreak: 'break-all', fontSize: '0.75rem' } }} />
                      {row.kind === 'template' && (
                        <Box sx={{ display: 'flex', gap: 0.5, mt: 0.5, flexWrap: 'wrap' }}>
                          {(group.masks[row.id] || []).map((mask, maskIndex) => (
                            <Chip key={`${mask.start}-${mask.end}-${maskIndex}`} size="small" label={`${mask.start + 1}–${mask.end}`}
                              onDelete={() => removeMask(groupIndex, row.id, maskIndex)} />
                          ))}
                          {!(group.masks[row.id] || []).length && <Typography variant="caption" color="text.secondary">Whole row selected by default</Typography>}
                        </Box>
                      )}
                    </Box>
                  ))}
                </Box>
              )
            })}
          </Box>
        )}

        {tab === 3 && (
          <Box sx={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))', gap: 1.5 }}>
            <TextField size="small" type="number" label="Models" value={project.modelOptions['--num-models'] ?? 5}
              onChange={event => update({ modelOptions: { ...project.modelOptions, '--num-models': Number(event.target.value) } })} />
            <FormControl size="small"><InputLabel>MD refinement</InputLabel><Select label="MD refinement" value={project.modelOptions['--md-level'] ?? 'fast'}
              onChange={event => update({ modelOptions: { ...project.modelOptions, '--md-level': event.target.value } })}>
              {['none', 'fast', 'slow', 'very_slow', 'slow_large'].map(level => <MenuItem key={level} value={level}>{level}</MenuItem>)}
            </Select></FormControl>
            <Button variant="contained" startIcon={busy ? <CircularProgress size={14} /> : <PlayArrowIcon />}
              onClick={runModel} disabled={busy !== '' || !project.alignmentGroups.length}>Run Modeller</Button>
            <Typography variant="caption" color="text.secondary" sx={{ gridColumn: '1 / -1' }}>
              Each painted span becomes an independent derived template fragment. Overlapping spans remain overlapping Modeller restraints. Rows without masks use their complete aligned coverage.
            </Typography>
          </Box>
        )}
      </Box>
    </Box>
  )
}
