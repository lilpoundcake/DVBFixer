import { create } from 'zustand'

export interface WorkspaceArtifact {
  id: string
  file: string
  name: string
  kind: 'structure' | 'artifact'
  folder?: string
  parent?: string
  command?: string
  artifactType?: string
  hidden?: boolean
}

export interface WorkspaceManifest {
  version: 1
  id: string
  name: string
  createdAt: string
  updatedAt: string
  primaryFile: string | null
  secondaryFile: string | null
  artifacts: WorkspaceArtifact[]
  toolState: Record<string, any>
}

interface WorkspaceSummary { id: string; name: string; updatedAt: string; artifactCount: number }

interface WorkspaceState {
  workspaces: WorkspaceSummary[]
  active: WorkspaceManifest | null
  loading: boolean
  error: string | null
  revision: number
  textPreview: { workspaceId: string; file: string; name: string } | null
  refresh: () => Promise<void>
  activate: (id: string) => Promise<WorkspaceManifest>
  reload: () => Promise<WorkspaceManifest | null>
  createWorkspace: (name?: string) => Promise<WorkspaceManifest>
  update: (patch: Partial<WorkspaceManifest>) => void
  updateToolState: (panel: string, value: unknown) => void
  save: () => Promise<void>
  discardActive: () => void
  setTextPreview: (preview: WorkspaceState['textPreview']) => void
  bumpRevision: () => void
}

let saveTimer: number | undefined

export function workspaceFileUrl(workspaceId: string, file: string): string {
  return `/api/workspaces/${encodeURIComponent(workspaceId)}/files/${file.split('/').map(encodeURIComponent).join('/')}`
}

export const useWorkspaceStore = create<WorkspaceState>((set, get) => ({
  workspaces: [], active: null, loading: false, error: null, revision: 0, textPreview: null,
  refresh: async () => {
    const response = await fetch('/api/workspaces', { cache: 'no-store' })
    if (!response.ok) throw new Error(`Workspaces: HTTP ${response.status}`)
    set({ workspaces: await response.json() as WorkspaceSummary[] })
  },
  activate: async (id) => {
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    if (get().active) await get().save()
    set({ loading: true, error: null })
    try {
      const response = await fetch(`/api/workspaces/${encodeURIComponent(id)}`, { cache: 'no-store' })
      const body = await response.json() as WorkspaceManifest & { error?: string }
      if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
      set(state => ({ active: body, loading: false, revision: state.revision + 1 }))
      return body
    } catch (reason: any) {
      set({ loading: false, error: reason.message || String(reason) })
      throw reason
    }
  },
  reload: async () => {
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    const id = get().active?.id
    if (!id) return null
    const response = await fetch(`/api/workspaces/${encodeURIComponent(id)}`, { cache: 'no-store' })
    const body = await response.json() as WorkspaceManifest & { error?: string }
    if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
    set(state => ({ active: body, revision: state.revision + 1 }))
    return body
  },
  createWorkspace: async (name = 'Untitled workspace') => {
    const response = await fetch('/api/workspaces', {
      method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ name }),
    })
    const body = await response.json() as WorkspaceManifest & { error?: string }
    if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
    await get().refresh()
    set(state => ({ active: body, revision: state.revision + 1 }))
    return body
  },
  update: (patch) => {
    set(state => ({ active: state.active ? { ...state.active, ...patch } : null }))
    if (saveTimer !== undefined) window.clearTimeout(saveTimer)
    saveTimer = window.setTimeout(() => { get().save().catch(() => {}) }, 700)
  },
  updateToolState: (panel, value) => {
    const active = get().active
    if (active) get().update({ toolState: { ...active.toolState, [panel]: value } })
  },
  save: async () => {
    const active = get().active
    if (!active) return
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    const response = await fetch(`/api/workspaces/${encodeURIComponent(active.id)}`, {
      method: 'PUT', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(active),
    })
    const body = await response.json() as WorkspaceManifest & { error?: string }
    if (!response.ok) { const error = body.error || `HTTP ${response.status}`; set({ error }); throw new Error(error) }
    set({ active: body })
    get().refresh().catch(() => {})
  },
  discardActive: () => {
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    set(state => ({ active: null, revision: state.revision + 1 }))
  },
  setTextPreview: (textPreview) => set({ textPreview }),
  bumpRevision: () => set(state => ({ revision: state.revision + 1 })),
}))
