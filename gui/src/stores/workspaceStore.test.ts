import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import { useWorkspaceStore, type WorkspaceManifest } from './workspaceStore'

const workspace: WorkspaceManifest = {
  version: 1,
  revision: 0,
  id: 'workspace-a',
  name: 'Workspace A',
  createdAt: '2026-01-01T00:00:00.000Z',
  updatedAt: '2026-01-01T00:00:00.000Z',
  primaryFile: null,
  secondaryFile: null,
  artifacts: [],
  toolState: {},
}

function jsonResponse(value: unknown): Response {
  return new Response(JSON.stringify(value), {
    status: 200,
    headers: { 'Content-Type': 'application/json' },
  })
}

describe('workspace initialization', () => {
  beforeEach(() => {
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

  afterEach(() => vi.unstubAllGlobals())

  it('deduplicates concurrent initialization and activates the first workspace once', async () => {
    const fetchMock = vi.fn(async (input: string | URL | Request) => {
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
