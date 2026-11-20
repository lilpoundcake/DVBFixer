import { afterEach, expect, it, vi } from 'vitest'
import { runAltlocOperation } from './altloc-operations'
import { apiFetch } from './api-client'
import { useWorkspaceStore, type WorkspaceManifest } from '../stores/workspaceStore'
import { useAltlocStore } from '../stores/altlocStore'
vi.mock('./api-client', () => ({ apiFetch: vi.fn() }))

afterEach(() => {
  vi.unstubAllGlobals(); vi.clearAllMocks()
  useAltlocStore.setState({ sessions: {} }); useWorkspaceStore.setState({ active: null })
})

it.each(['analyze', 'run', 'apply'] as const)('opens review immediately on %s, never again at completion', async action => {
  const dispatchEvent = vi.fn()
  vi.stubGlobal('window', { dispatchEvent })
  useWorkspaceStore.setState({ active: { id: 'w' } as WorkspaceManifest })
  useAltlocStore.getState().patch('w', { input: 'structure' })
  let finish!: (response: Response) => void
  vi.mocked(apiFetch).mockImplementationOnce(() => new Promise(resolve => { finish = resolve }))
  const operation = runAltlocOperation('w', action)
  expect(dispatchEvent).toHaveBeenCalledTimes(1)
  expect(dispatchEvent.mock.calls[0][0].type).toBe('dvbfixer:open-altlocs')
  useWorkspaceStore.setState({ active: { id: 'other' } as WorkspaceManifest })
  if (action !== 'apply') finish(Response.json({ report: { residues: [] } }))
  if (action === 'run') vi.mocked(apiFetch).mockResolvedValueOnce(Response.json({ outputArtifact: { name: 'selected.pdb' } }))
  await operation
  expect(dispatchEvent).toHaveBeenCalledTimes(1)
  expect(useAltlocStore.getState().sessions.w.busy).toBe(false)
})
