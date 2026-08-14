import { describe, expect, it } from 'vitest'
import { createManagedJobRequest, isActiveManagedJob, managedJobStatusLabel, selectRestoredManagedJob, type ManagedJobRecord } from './managed-jobs'

function job(status: ManagedJobRecord['status'], id = status): ManagedJobRecord {
  return {
    version: 1, id, workspaceId: 'workspace', status,
    createdAt: '2026-01-01T00:00:00Z', startedAt: null, finishedAt: null,
    command: 'prepare', args: [], exitCode: status === 'failed' ? 2 : null,
    stdoutLog: 'runs/stdout.log', stderrLog: 'runs/stderr.log',
    outputDir: 'runs/job', outputFile: status === 'succeeded' ? 'runs/job/out.pdb' : null,
  }
}

describe('managed job UI helpers', () => {
  it('includes the selected command in the managed-job request', () => {
    expect(createManagedJobRequest('prepare', {
      workspaceId: 'workspace-1', inputFile: 'input.pdb', inputs: {},
      values: {}, fastaContent: '',
    })).toMatchObject({ command: 'prepare', workspaceId: 'workspace-1' })
  })

  it('restores an active job ahead of a newer terminal record', () => {
    expect(selectRestoredManagedJob([job('succeeded'), job('running')])?.status).toBe('running')
    expect(selectRestoredManagedJob([])).toBeNull()
  })

  it('identifies only queued and running records as active', () => {
    expect(isActiveManagedJob(job('queued'))).toBe(true)
    expect(isActiveManagedJob(job('running'))).toBe(true)
    expect(isActiveManagedJob(job('cancelled'))).toBe(false)
  })

  it('formats concise terminal status', () => {
    expect(managedJobStatusLabel(job('succeeded'))).toContain('out.pdb')
    expect(managedJobStatusLabel(job('failed'))).toContain('exit 2')
  })
})
