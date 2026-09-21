import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { useWorkspaceStore, type WorkspaceManifest } from './workspaceStore'
import { API_TOKEN_STORAGE_KEY } from '../lib/api-client'

const workspace: WorkspaceManifest = {
  version: 2,
  revision: 0,
  id: 'workspace-a',
  name: 'Workspace A',
  createdAt: '2026-01-01T00:00:00.000Z',
  updatedAt: '2026-01-01T00:00:00.000Z',
  primaryFile: null,
  secondaryFile: null,
  artifacts: [],
  toolState: {},
  ownerPrincipalId: 'owner',
  acl: [{ principalId: 'reader', role: 'reader' }],
  effectiveRole: 'owner',
}

function jsonResponse(value: unknown): Response {
  return new Response(JSON.stringify(value), {
    status: 200,
    headers: { 'Content-Type': 'application/json' },
  })
}

function browserWindow(token: string | null = null) {
  return {
    location: { href: 'https://dvbfixer.test/' },
    sessionStorage: { getItem: (key: string) => key === API_TOKEN_STORAGE_KEY ? token : null },
    setTimeout: globalThis.setTimeout,
    clearTimeout: globalThis.clearTimeout,
  }
}

describe('workspace initialization', () => {
  it('does not queue autosaves for a read-only workspace', () => {
    const readOnly = { ...workspace, effectiveRole: 'reader' as const }
    useWorkspaceStore.setState({ active: readOnly })
    useWorkspaceStore.getState().updateToolState('dvbfixer', { inputFile: 'changed.pdb' })
    expect(useWorkspaceStore.getState().active?.toolState).toEqual({})
  })

  beforeEach(() => {
    vi.stubGlobal('window', browserWindow())
    useWorkspaceStore.setState({
      workspaces: [],
      active: null,
      loading: false,
      error: null,
      initialized: false,
      revision: 0,
      textPreview: null,
    })
  })

  afterEach(() => {
    vi.useRealTimers()
    vi.unstubAllGlobals()
  })

  it('rebases Split autosave after a job adds output without losing server state or newer edits', async () => {
    vi.useFakeTimers()
    useWorkspaceStore.setState({ active: { ...workspace, revision: 142 }, initialized: true })
    const output = { id: 'split-output', file: 'split.pdb', name: 'Split', kind: 'structure' as const }
    let latest = {
      ...workspace, revision: 143, artifacts: [output], primaryFile: output.file,
      toolState: { homology: { tab: 2 } } as WorkspaceManifest['toolState'],
    }
    const patches: any[] = []
    vi.stubGlobal('fetch', vi.fn(async (input: string, init?: RequestInit) => {
      if (input === '/api/workspaces') return jsonResponse([])
      if (init?.method === 'PATCH') {
        patches.push(JSON.parse(String(init.body)))
        if (patches.length === 1) return new Response(JSON.stringify({ error: 'revision conflict' }), { status: 409 })
        latest = { ...latest, revision: 144, toolState: patches[1].toolState }
        return jsonResponse(latest)
      }
      // Simulate another edit while the conflict recovery GET is running.
      if (latest.revision === 143) useWorkspaceStore.getState().updateToolState('dvbfixer', { command: 'split', renumber: false })
      return jsonResponse(latest)
    }))

    useWorkspaceStore.getState().updateToolState('dvbfixer', { command: 'split', renumber: true })
    await useWorkspaceStore.getState().reload()

    expect(patches).toHaveLength(2)
    expect(patches[1]).toEqual({
      revision: 143,
      toolState: { homology: { tab: 2 }, dvbfixer: { command: 'split', renumber: false } },
    })
    expect(useWorkspaceStore.getState().active).toMatchObject({
      artifacts: [output], primaryFile: output.file,
      toolState: patches[1].toolState,
    })
    expect(useWorkspaceStore.getState().error).toBeNull()
    // Flush the timer scheduled by the edit inside the recovery request.
    await useWorkspaceStore.getState().save()
  })

  it('deduplicates concurrent initialization and activates the first workspace once', async () => {
    vi.stubGlobal('window', browserWindow('workspace-token'))
    const fetchMock = vi.fn(async (input: string | URL | Request, _init?: RequestInit) => {
      const url = String(input)
      if (url === '/api/workspaces') {
        return jsonResponse([{ id: workspace.id, name: workspace.name, updatedAt: workspace.updatedAt, artifactCount: 0 }])
      }
      if (url === `/api/workspaces/${workspace.id}`) return jsonResponse(workspace)
      throw new Error(`Unexpected request: ${url}`)
    })
    vi.stubGlobal('fetch', fetchMock)

    await Promise.all([
      useWorkspaceStore.getState().initialize(),
      useWorkspaceStore.getState().initialize(),
    ])

    expect(fetchMock).toHaveBeenCalledTimes(2)
    expect(useWorkspaceStore.getState()).toMatchObject({ initialized: true, active: workspace })
    for (const [, init] of fetchMock.mock.calls) {
      expect(new Headers(init?.headers).get('Authorization')).toBe('Bearer workspace-token')
    }
  })

  it('saves a narrow revisioned patch without replacing newer local state', async () => {
    useWorkspaceStore.setState({ active: workspace, initialized: true })
    let resolveResponse!: (response: Response) => void
    const responsePromise = new Promise<Response>(resolve => { resolveResponse = resolve })
    const fetchMock = vi.fn((_input: string | URL | Request, _init?: RequestInit) => responsePromise)
    vi.stubGlobal('fetch', fetchMock)

    const saving = useWorkspaceStore.getState().save()
    useWorkspaceStore.setState({
      active: { ...workspace, toolState: { homology: { tab: 3 } } },
    })
    resolveResponse(jsonResponse({ ...workspace, revision: 1, updatedAt: '2026-01-02T00:00:00.000Z' }))
    await saving

    const [, request] = fetchMock.mock.calls[0]
    expect(request).toMatchObject({ method: 'PATCH' })
    expect(JSON.parse(String(request?.body))).toEqual({
      revision: 0,
      toolState: {},
      primaryFile: null,
      secondaryFile: null,
      artifactOrder: [],
    })
    expect(useWorkspaceStore.getState().active).toMatchObject({
      revision: 1,
      toolState: { homology: { tab: 3 } },
    })
  })

  it('surfaces revision conflicts without replacing the active workspace', async () => {
    useWorkspaceStore.setState({ active: workspace, initialized: true })
    vi.stubGlobal('fetch', vi.fn(async () => new Response(JSON.stringify({ error: 'workspace revision conflict' }), {
      status: 409,
      headers: { 'Content-Type': 'application/json' },
    })))

    await expect(useWorkspaceStore.getState().save()).rejects.toThrow(/save conflict/i)
    expect(useWorkspaceStore.getState().active).toBe(workspace)
    expect(useWorkspaceStore.getState().error).toMatch(/reload the workspace/i)
  })
})
