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
  organism?: string
  method?: string
  resolution?: string | number
  description?: string
  iggSubtype?: string
  allotype?: string
  equivalentChains?: string[][]
}

export interface WorkspaceArtifactMetadataPatch {
  name?: string
  organism?: string | null
  method?: string | null
  resolution?: string | number | null
  description?: string | null
  iggSubtype?: string | null
  allotype?: string | null
  equivalentChains?: string[][] | null
}

export interface WorkspaceManifest {
  version: 1
  revision: number
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
  initialized: boolean
  revision: number
  textPreview: { workspaceId: string; file: string; name: string } | null
  refresh: () => Promise<void>
  initialize: () => Promise<void>
  activate: (id: string) => Promise<WorkspaceManifest>
  reload: () => Promise<WorkspaceManifest | null>
  createWorkspace: (name?: string) => Promise<WorkspaceManifest>
  reorderWorkspaces: (ids: string[]) => Promise<void>
  update: (patch: Partial<WorkspaceManifest>) => void
  updateToolState: (panel: string, value: unknown) => void
  save: () => Promise<void>
  updateArtifactMetadata: (artifactId: string, metadata: WorkspaceArtifactMetadataPatch) => Promise<WorkspaceManifest>
  discardActive: () => void
  setTextPreview: (preview: WorkspaceState['textPreview']) => void
  bumpRevision: () => void
}

let saveTimer: number | undefined
let saveInFlight: Promise<void> | undefined
let saveRequested = false
let metadataWriteInFlight: Promise<WorkspaceManifest> | undefined
let initializePromise: Promise<void> | undefined
let refreshGeneration = 0
let activationGeneration = 0

function staleRequest(message: string): Error {
  const error = new Error(message)
  error.name = 'AbortError'
  return error
}

export function workspaceFileUrl(workspaceId: string, file: string): string {
  return `/api/workspaces/${encodeURIComponent(workspaceId)}/files/${file.split('/').map(encodeURIComponent).join('/')}`
}

export const useWorkspaceStore = create<WorkspaceState>((set, get) => ({
  workspaces: [], active: null, loading: false, error: null, initialized: false, revision: 0, textPreview: null,
  refresh: async () => {
    const generation = ++refreshGeneration
    const response = await fetch('/api/workspaces', { cache: 'no-store' })
    if (!response.ok) throw new Error(`Workspaces: HTTP ${response.status}`)
    const workspaces = await response.json() as WorkspaceSummary[]
    if (generation === refreshGeneration) set({ workspaces })
  },
  initialize: async () => {
    if (get().initialized) return
    if (initializePromise) return initializePromise
    initializePromise = (async () => {
      set({ loading: true, error: null })
      try {
        await get().refresh()
        const current = get()
        if (!current.active) {
          if (current.workspaces[0]) await get().activate(current.workspaces[0].id)
          else await get().createWorkspace()
        }
        set({ initialized: true, loading: false })
      } catch (reason: any) {
        set({ loading: false, error: reason.message || String(reason) })
        throw reason
      } finally {
        initializePromise = undefined
      }
    })()
    return initializePromise
  },
  activate: async (id) => {
    const generation = ++activationGeneration
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    if (get().active) await get().save()
    if (generation !== activationGeneration) throw staleRequest('Workspace activation was superseded')
    set({ loading: true, error: null })
    try {
      const response = await fetch(`/api/workspaces/${encodeURIComponent(id)}`, { cache: 'no-store' })
      const body = await response.json() as WorkspaceManifest & { error?: string }
      if (!response.ok) throw new Error(body.error || `HTTP ${response.status}`)
      if (generation !== activationGeneration) throw staleRequest('Workspace activation was superseded')
      set(state => ({ active: body, loading: false, revision: state.revision + 1 }))
      return body
    } catch (reason: any) {
      if (generation === activationGeneration) set({ loading: false, error: reason.message || String(reason) })
      throw reason
    }
  },
  reload: async () => {
    // Reload is a user-visible synchronization point: flush pending local
    // edits instead of silently cancelling the debounced save.
    if (saveTimer !== undefined || saveInFlight) await get().save()
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
  reorderWorkspaces: async (ids) => {
    const byId = new Map(get().workspaces.map(workspace => [workspace.id, workspace]))
    set({ workspaces: ids.flatMap(id => byId.has(id) ? [byId.get(id)!] : []) })
    const response = await fetch('/api/workspaces/order', {
      method: 'PATCH', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ ids }),
    })
    if (!response.ok) throw new Error(`Workspace reorder: HTTP ${response.status}`)
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
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    if (!get().active) return
    if (metadataWriteInFlight) {
      await metadataWriteInFlight
      return get().save()
    }
    if (saveInFlight) {
      saveRequested = true
      return saveInFlight
    }
    saveInFlight = (async () => {
      do {
        saveRequested = false
        const active = get().active
        if (!active) return
        const response = await fetch(`/api/workspaces/${encodeURIComponent(active.id)}`, {
          method: 'PATCH',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            revision: active.revision,
            toolState: active.toolState,
            primaryFile: active.primaryFile,
            secondaryFile: active.secondaryFile,
            artifactOrder: active.artifacts.map(artifact => artifact.id),
          }),
        })
        const body = await response.json() as WorkspaceManifest & { error?: string }
        if (!response.ok) {
          const detail = body.error || `HTTP ${response.status}`
          const error = response.status === 409
            ? `Workspace save conflict. Reload the workspace before saving again. ${detail}`
            : detail
          set({ error })
          throw new Error(error)
        }
        // A save response is only an acknowledgement. Preserve any newer
        // local edits and merge back just the concurrency token/timestamp.
        set(state => state.active?.id === active.id
          ? { active: { ...state.active, revision: body.revision, updatedAt: body.updatedAt }, error: null }
          : state)
        get().refresh().catch(() => {})
      } while (saveRequested)
    })()
    try {
      await saveInFlight
    } finally {
      saveInFlight = undefined
    }
  },
  updateArtifactMetadata: async (artifactId, metadata) => {
    if (saveTimer !== undefined || saveInFlight) await get().save()
    if (metadataWriteInFlight) await metadataWriteInFlight
    const active = get().active
    if (!active) throw new Error('No active workspace')
    metadataWriteInFlight = (async () => {
      const response = await fetch(`/api/workspaces/${encodeURIComponent(active.id)}/artifacts/${encodeURIComponent(artifactId)}/metadata`, {
        method: 'PATCH',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ revision: active.revision, ...metadata }),
      })
      const body = await response.json() as WorkspaceManifest & { error?: string }
      if (!response.ok) {
        const detail = body.error || `HTTP ${response.status}`
        const error = response.status === 409
          ? `Workspace metadata conflict. Reload the workspace before saving again. ${detail}`
          : detail
        set({ error })
        throw new Error(error)
      }
      const returnedArtifact = body.artifacts.find(artifact => artifact.id === artifactId)
      set(state => state.active?.id === active.id && returnedArtifact
        ? {
            active: {
              ...state.active,
              revision: body.revision,
              updatedAt: body.updatedAt,
              artifacts: state.active.artifacts.map(artifact => artifact.id === artifactId ? returnedArtifact : artifact),
            },
            error: null,
          }
        : state)
      get().refresh().catch(() => {})
      return body
    })()
    try {
      return await metadataWriteInFlight
    } finally {
      metadataWriteInFlight = undefined
    }
  },
  discardActive: () => {
    activationGeneration++
    if (saveTimer !== undefined) { window.clearTimeout(saveTimer); saveTimer = undefined }
    set(state => ({ active: null, revision: state.revision + 1 }))
  },
  setTextPreview: (textPreview) => set({ textPreview }),
  bumpRevision: () => set(state => ({ revision: state.revision + 1 })),
}))
