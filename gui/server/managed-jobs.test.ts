import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import {
  cancelManagedJob, createManagedJob, getManagedJob, listManagedJobs,
  type ManagedJobRecord,
} from './managed-jobs'
import { ensureWorkspaceMigration, listWorkspaces, workspaceRoot } from './workspace-api'

const directories: string[] = []
const originalEnvironment = {
  executable: process.env.DVBFIXER_EXECUTABLE,
  args: process.env.DVBFIXER_ARGS,
}

function root(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-managed-jobs-'))
  directories.push(directory)
  ensureWorkspaceMigration(directory)
  return directory
}

function useNode(script: string): void {
  process.env.DVBFIXER_EXECUTABLE = process.execPath
  process.env.DVBFIXER_ARGS = JSON.stringify(['-e', script])
}

async function terminal(
  dataRoot: string,
  workspaceId: string,
  jobId: string,
  timeout = 3000,
): Promise<ManagedJobRecord> {
  const deadline = Date.now() + timeout
  while (Date.now() < deadline) {
    const record = getManagedJob(dataRoot, workspaceId, jobId)
    if (['succeeded', 'failed', 'cancelled'].includes(record.status)) return record
    await new Promise(resolve => setTimeout(resolve, 10))
  }
  throw new Error(`job ${jobId} did not finish`)
}

afterEach(() => {
  for (const [key, value] of Object.entries({
    DVBFIXER_EXECUTABLE: originalEnvironment.executable,
    DVBFIXER_ARGS: originalEnvironment.args,
  })) {
    if (value === undefined) delete process.env[key]
    else process.env[key] = value
  }
  for (const directory of directories.splice(0)) fs.rmSync(directory, { recursive: true, force: true })
})

describe('managed DVBfixer jobs', () => {
  it('persists a successful fake-CLI lifecycle and logs', async () => {
    const dataRoot = root()
    const workspace = listWorkspaces(dataRoot)[0]
    useNode("console.log('managed ok')")
    const created = createManagedJob(dataRoot, { workspaceId: workspace.id, command: 'doctor' })
    expect(created.status).toBe('queued')
    expect(created.id).toMatch(/^[0-9a-f-]{36}$/)

    const finished = await terminal(dataRoot, workspace.id, created.id)
    expect(finished.status).toBe('succeeded')
    expect(finished.startedAt).toBeTruthy()
    expect(finished.finishedAt).toBeTruthy()
    expect(finished.exitCode).toBe(0)
    expect(finished.args).toEqual([])
    const workspaceDirectory = workspaceRoot(dataRoot, workspace.id)
    expect(fs.readFileSync(path.join(workspaceDirectory, finished.stdoutLog), 'utf8')).toContain('managed ok')
    expect(fs.existsSync(path.join(workspaceDirectory, finished.outputDir, 'job.json'))).toBe(true)
    expect(listManagedJobs(dataRoot, workspace.id).map(job => job.id)).toContain(created.id)
  })

  it('rejects concurrent jobs in one workspace and cancels the active process', async () => {
    const dataRoot = root()
    const workspace = listWorkspaces(dataRoot)[0]
    useNode('setInterval(() => {}, 1000)')
    const created = createManagedJob(dataRoot, { workspaceId: workspace.id, command: 'doctor' })
    expect(() => createManagedJob(dataRoot, { workspaceId: workspace.id, command: 'doctor' }))
      .toThrow(/already running/)
    cancelManagedJob(dataRoot, workspace.id, created.id)
    const finished = await terminal(dataRoot, workspace.id, created.id)
    expect(finished.status).toBe('cancelled')
    expect(finished.finishedAt).toBeTruthy()
    expect(fs.readFileSync(path.join(workspaceRoot(dataRoot, workspace.id), finished.stderrLog), 'utf8'))
      .toContain('cancelled')
  })
})
