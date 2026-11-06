import { expect, it } from 'vitest'
import { selectHomologyJob } from './homology-jobs'
import { acceptManagedJobUpdate, type ManagedJobRecord } from './managed-jobs'

const job = (id: string, projectId: string, workspaceId = 'workspace'): ManagedJobRecord => ({
  id, workspaceId, command: 'homology', status: 'running', version: 1, revision: 2,
  createdAt: '2026-01-01', startedAt: null, finishedAt: null, args: [], exitCode: null,
  stdoutLog: '', stderrLog: '', outputDir: '', outputFile: null,
  workflow: { kind: 'homology', projectId, projectName: projectId },
})

it('restores only the selected Homology project in the selected workspace', () => {
  const selected = job('selected', 'VH')
  expect(selectHomologyJob([job('lower', 'vh'), job('other', 'VH', 'elsewhere'), selected], 'workspace', 'VH')).toBe(selected)
  expect(selectHomologyJob([selected], 'workspace', 'VL')).toBeNull()
})

it('prefers active work over history and rejects stale GET/SSE revisions', () => {
  const active = job('active', 'H')
  const finished = { ...job('finished', 'H'), status: 'succeeded' as const, createdAt: '2026-02-01' }
  expect(selectHomologyJob([finished, active], 'workspace', 'H')).toBe(active)
  expect(acceptManagedJobUpdate(active, { ...active, revision: 1 })).toBe(false)
  expect(acceptManagedJobUpdate(active, { ...active, revision: 3, status: 'succeeded' })).toBe(true)
  expect(acceptManagedJobUpdate({ ...finished, revision: undefined }, { ...finished, revision: undefined, status: 'running' })).toBe(false)
})
