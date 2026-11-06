import { EventEmitter } from 'node:events'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import crypto from 'node:crypto'
import { Readable } from 'node:stream'
import { afterEach, describe, expect, it, vi } from 'vitest'
import { registerHomologyApi, registerHomologyJobApi, type HomologyProject } from './homology-api'
import { cancelManagedJob, createManagedJob, getManagedJob, listManagedJobs, registerManagedJobApi, type ManagedJobRecord } from './managed-jobs'
import { ensureWorkspaceMigration, listWorkspaces, loadWorkspace, saveWorkspace, workspaceRoot, createWorkspace } from './workspace-api'
import * as runStorage from './run-storage'
import { resetApiShutdown, shutdownApiWork } from './api-routes'

const roots: string[] = [], jobs: Array<{ root: string; job: ManagedJobRecord }> = []
const original = { executable: process.env.DVBFIXER_EXECUTABLE, args: process.env.DVBFIXER_ARGS }
async function terminal(root: string, job: ManagedJobRecord) {
  let record = job
  await vi.waitFor(() => {
    record = getManagedJob(root, job.workspaceId, job.id)
    expect(['succeeded', 'failed', 'cancelled']).toContain(record.status)
  }, { timeout: 5000 })
  return record
}
afterEach(async () => {
  for (const item of jobs) cancelManagedJob(item.root, item.job.workspaceId, item.job.id)
  await Promise.all(jobs.splice(0).map(item => terminal(item.root, item.job)))
  vi.restoreAllMocks()
  for (const [key, value] of Object.entries({ DVBFIXER_EXECUTABLE: original.executable, DVBFIXER_ARGS: original.args })) {
    if (value === undefined) delete process.env[key]; else process.env[key] = value
  }
  for (const root of roots.splice(0)) fs.rmSync(root, { recursive: true, force: true })
})

function fixture() {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'homology-job-')); roots.push(root)
  ensureWorkspaceMigration(root)
  const workspace = listWorkspaces(root)[0]
  const directory = workspaceRoot(root, workspace.id)
  fs.writeFileSync(path.join(directory, 'template.pdb'), 'END\n')
  const project: HomologyProject = { version: 1, id: 'project-H', name: 'test project', targetFasta: '>H\nMA\n',
    templates: [{ id: 'template-1', file: 'template.pdb', chain: 'H', targetChain: 'H' }], engine: 'mafft',
    alignmentGroups: [{ chainId: 'H', rows: [{ id: 'H', kind: 'target', sequence: 'MA' }], masks: {} }],
    modelOptions: {}, createdAt: '2026-01-01', updatedAt: '2026-01-01' }
  process.env.DVBFIXER_EXECUTABLE = process.execPath
  process.env.DVBFIXER_ARGS = JSON.stringify(['-e', `
    const fs=require('node:fs');
    const fd=Number(process.env.DVBFIXER_JOB_PROGRESS_FD);
    fs.writeSync(fd, JSON.stringify({version:1,stage:'sampling',seed:7})+'\\n');
    console.log('synthetic orchestration control; no model weights');
    const command=process.argv[1],out=process.argv[process.argv.indexOf('-o')+1];
    if(command==='homology'){
      const plan=JSON.parse(fs.readFileSync(process.argv[process.argv.indexOf('--template-plan')+1],'utf8'));
      const bytes=fs.readFileSync(plan.templates[0].path);
      setInterval(()=>{if(fs.existsSync('finish')){fs.writeFileSync(out+'.pdb',bytes);fs.writeFileSync(out+'.dat','{}');process.exit(0)}},10);
    }else setInterval(()=>{if(fs.existsSync('finish'))process.exit(0)},10);
  `])
  return { root, workspace, directory, project }
}

async function api(root: string, route: string, body?: unknown, principal = 'local', method = body ? 'POST' : 'GET', contentType = 'application/json') {
  let middleware: any
  const server = { middlewares: { use: (_prefix: string, handler: any) => { middleware = handler } } } as any
  if (route.startsWith('/legacy/')) registerHomologyApi(server, root, () => principal)
  else if (route.endsWith('/homology/jobs')) registerHomologyJobApi(server, root, () => principal, 'local', 'test-version')
  else registerManagedJobApi(server, root, () => principal)
  const request = Readable.from(body === undefined ? [] : [typeof body === 'string' ? body : JSON.stringify(body)]) as any
  request.method = method; request.url = route.replace('/legacy', ''); request.headers = { 'content-type': contentType }
  let text = ''
  const response = Object.assign(new EventEmitter(), { statusCode: 0, writableFinished: false,
    setHeader: () => {}, writeHead(status: number) { this.statusCode = status }, write(value: string) { text += value },
    end(value?: string) { if (value) text += value; this.writableFinished = true; this.emit('close') },
  })
  await middleware(request, response, () => { throw new Error(`unhandled route: ${route}`) })
  return { status: response.statusCode, body: text.startsWith('data:') ? text : text ? JSON.parse(text) : null, response }
}

async function start(root: string, workspaceId: string, project: HomologyProject) {
  const response = await api(root, `/workspaces/${workspaceId}/homology/jobs`, { project })
  expect(response.status).toBe(202)
  jobs.push({ root, job: response.body })
  return response
}

describe('persistent Homology jobs', () => {
  it('survives HTTP disconnect, freezes inputs and publishes through shared jobs', async () => {
    const { root, workspace, directory, project } = fixture()
    const accepted = await start(root, workspace.id, project)
    accepted.response.emit('close')
    const job = accepted.body as ManagedJobRecord
    await vi.waitFor(() => expect(getManagedJob(root, workspace.id, job.id)).toMatchObject({
      status: 'running', progress: { stage: 'sampling', seed: 7 }, workflow: { projectId: project.id },
    }))
    const run = path.join(directory, job.outputDir)
    const snapshot = fs.readFileSync(path.join(run, 'homology-project.json'))
    expect(job.workflow?.snapshotSha256).toBe(crypto.createHash('sha256').update(snapshot).digest('hex'))
    const plan = JSON.parse(fs.readFileSync(path.join(run, 'template-plan.json'), 'utf8'))
    expect(plan.templates[0].path).toBe(path.join(run, 'inputs/template-1.pdb'))
    fs.writeFileSync(path.join(directory, 'template.pdb'), 'CHANGED AFTER SUBMISSION\n')
    expect(fs.readFileSync(plan.templates[0].path, 'utf8')).toBe('END\n')
    fs.writeFileSync(path.join(run, 'finish'), 'done')
    const finished = await terminal(root, job)
    expect(finished).toMatchObject({ status: 'succeeded', modelStatus: 'success', workflow: { projectId: 'project-H' } })
    expect(finished.revision).toBeGreaterThan(job.revision!)
    expect(fs.readFileSync(path.join(directory, finished.outputFile!), 'utf8')).toBe('END\n')
    expect(finished.provenance?.resolvedInputs).toEqual([{ artifactId: null, file: 'template.pdb',
      sha256: crypto.createHash('sha256').update('END\n').digest('hex') }])
    expect(finished.provenance?.outputs.length).toBeGreaterThan(0)
    const artifact = loadWorkspace(root, workspace.id).artifacts.find(item => item.file === finished.outputFile)
    expect(artifact?.workflowProvenance?.request.values.homologyProject).toEqual({ id: project.id, snapshotSha256: job.workflow?.snapshotSha256 })
    const summary = JSON.parse(fs.readFileSync(path.join(run, 'run-summary.json'), 'utf8'))
    expect(summary.diagnostics.files.find((file: any) => file.file === 'homology-project.json').sha256).toBe(job.workflow?.snapshotSha256)
    expect((await api(root, `/workspaces/${workspace.id}/jobs/${job.id}`)).body).toMatchObject(finished)
    expect((await api(root, '/jobs')).body.jobs.some((item: any) => item.id === job.id)).toBe(true)
  })

  it('uses the same workspace lease for managed, legacy and ordinary jobs', async () => {
    const { root, workspace, project } = fixture()
    const { body: job } = await start(root, workspace.id, project)
    expect((await api(root, `/workspaces/${workspace.id}/homology/jobs`, { project })).status).toBe(409)
    expect(() => createManagedJob(root, { workspaceId: workspace.id, command: 'doctor' })).toThrow(/already running/)
    expect((await api(root, '/legacy/model', { ...project, workspaceId: workspace.id })).status).toBe(409)
    const stopped = await api(root, `/workspaces/${workspace.id}/jobs/${job.id}`, undefined, 'local', 'DELETE')
    expect(stopped.status).toBe(202)
    expect(stopped.body.cancellationRequestedAt).toMatch(/^\d{4}-/)
    expect((await terminal(root, job)).status).toBe('cancelled')
    const manifest = loadWorkspace(root, workspace.id)
    expect(manifest.artifacts.every(artifact => artifact.kind !== 'structure')).toBe(true)
  })

  it('cancels a queued Homology job without launching its CLI or cancelling the other workspace', async () => {
    const { root, workspace, directory, project } = fixture()
    const other = createWorkspace(root, 'other')
    const blocker = createManagedJob(root, { workspaceId: other.id, command: 'doctor' })
    jobs.push({ root, job: blocker })
    await vi.waitFor(() => expect(getManagedJob(root, other.id, blocker.id).status).toBe('running'))
    const { body: job } = await start(root, workspace.id, project)
    expect(getManagedJob(root, workspace.id, job.id).status).toBe('queued')
    cancelManagedJob(root, workspace.id, job.id)
    const cancelled = await terminal(root, job)
    expect(cancelled).toMatchObject({ status: 'cancelled', startedAt: null, outputFile: null })
    expect(fs.readFileSync(path.join(directory, cancelled.stdoutLog), 'utf8')).toBe('')
    expect(getManagedJob(root, other.id, blocker.id).status).toBe('running')
  })

  it.each(['schema', 'options', 'path', 'media', 'json', 'method'])('refuses invalid %s before reserving a job', async mode => {
    const { root, workspace, project } = fixture()
    let body: unknown = { project }, contentType = 'application/json'
    if (mode === 'schema') body = { project: { ...project, injected: true } }
    if (mode === 'options') project.modelOptions = { '--arbitrary-executable': 'bad' }
    if (mode === 'path') project.templates[0].file = '../private.pdb'
    if (mode === 'media') contentType = 'text/plain'
    if (mode === 'json') body = '{invalid'
    const response = await api(root, `/workspaces/${workspace.id}/homology/jobs`, body, 'local', mode === 'method' ? 'GET' : 'POST', contentType)
    expect(response.status).toBe(mode === 'media' ? 415 : mode === 'method' ? 405 : 400)
    expect(response.body.error).toMatchObject({ code: expect.any(String), requestId: expect.any(String) })
    expect(listManagedJobs(root, workspace.id)).toEqual([])
  })

  it('retains a failed child and its partial structure only as diagnostic evidence', async () => {
    const { root, workspace, directory, project } = fixture()
    process.env.DVBFIXER_ARGS = JSON.stringify(['-e', `const fs=require('node:fs');
      fs.writeFileSync(process.argv[process.argv.indexOf('-o')+1]+'.pdb','partial synthetic output');
      console.error('synthetic child failure');process.exit(2)`])
    const { body: job } = await start(root, workspace.id, project)
    const finished = await terminal(root, job)
    expect(finished).toMatchObject({ status: 'failed', exitCode: 2, modelStatus: 'failed', outputFile: null })
    expect(finished.error).toContain('synthetic child failure')
    const manifest = loadWorkspace(root, workspace.id)
    expect(manifest.artifacts.every(artifact => artifact.kind !== 'structure')).toBe(true)
    const summary = JSON.parse(fs.readFileSync(path.join(directory, job.outputDir, 'run-summary.json'), 'utf8'))
    expect(summary.diagnostics.files.some((file: any) => file.file.endsWith('.pdb'))).toBe(true)
  })

  it.each(['homology', 'model'])('does not publish accepted outputs when %s is stopped during packaging', async command => {
    const { root, workspace, directory, project } = fixture()
    let packed!: () => void, resume!: () => void
    const packaged = new Promise<void>(resolve => { packed = resolve })
    const paused = new Promise<void>(resolve => { resume = resolve })
    const originalCompact = runStorage.compactRun
    vi.spyOn(runStorage, 'compactRun').mockImplementation(async (...args) => {
      const summary = await originalCompact(...args)
      packed(); await paused
      return summary
    })
    let job: ManagedJobRecord
    if (command === 'homology') job = (await start(root, workspace.id, project)).body
    else {
      process.env.DVBFIXER_ARGS = JSON.stringify(['-e', `require('node:fs').writeFileSync(process.argv[process.argv.indexOf('-o')+1],'END\\n')`])
      job = createManagedJob(root, { workspaceId: workspace.id, command: 'model', inputFile: 'template.pdb' })
      jobs.push({ root, job })
    }
    try {
      if (command === 'homology') {
        await vi.waitFor(() => expect(getManagedJob(root, workspace.id, job.id).status).toBe('running'))
        fs.writeFileSync(path.join(directory, job.outputDir, 'finish'), 'done')
      }
      await packaged
      cancelManagedJob(root, workspace.id, job.id)
    } finally { resume() }
    const cancelled = await terminal(root, job)
    expect(cancelled).toMatchObject({ status: 'cancelled', outputFile: null })
    expect(loadWorkspace(root, workspace.id).artifacts.every(artifact => artifact.kind !== 'structure')).toBe(true)
    const summary = JSON.parse(fs.readFileSync(path.join(directory, job.outputDir, 'run-summary.json'), 'utf8'))
    expect(summary).toMatchObject({ status: 'cancelled', outputs: [] })
  })

  it('drains managed packaging and retained diagnostics before server shutdown completes', async () => {
    const { root, workspace, directory, project } = fixture()
    let packed!: () => void, resume!: () => void
    const packaged = new Promise<void>(resolve => { packed = resolve })
    const paused = new Promise<void>(resolve => { resume = resolve })
    const originalCompact = runStorage.compactRun
    vi.spyOn(runStorage, 'compactRun').mockImplementation(async (...args) => {
      const summary = await originalCompact(...args)
      packed(); await paused
      return summary
    })
    const { body: job } = await start(root, workspace.id, project)
    let shutdown: Promise<void> | undefined
    try {
      await vi.waitFor(() => expect(getManagedJob(root, workspace.id, job.id).status).toBe('running'))
      fs.writeFileSync(path.join(directory, job.outputDir, 'finish'), 'done')
      await packaged
      let drained = false
      shutdown = shutdownApiWork().then(() => { drained = true })
      await new Promise(resolve => setTimeout(resolve, 10))
      expect(drained).toBe(false)
      resume()
      await shutdown
      expect((await terminal(root, job)).status).toBe('cancelled')
      expect(fs.existsSync(path.join(directory, job.outputDir, 'run-summary.json'))).toBe(true)
    } finally {
      resume(); await shutdown
      resetApiShutdown()
    }
  })

  it('checks ownership before resolving inputs and prevents publication after writer revocation', async () => {
    const { root, workspace, directory, project } = fixture()
    saveWorkspace(root, { ...loadWorkspace(root, workspace.id), acl: [{ principalId: 'writer', role: 'writer' }, { principalId: 'reader', role: 'reader' }] })
    expect((await api(root, `/workspaces/${workspace.id}/homology/jobs`, { project }, 'foreign')).status).toBe(404)
    expect((await api(root, `/workspaces/${workspace.id}/homology/jobs`, { project }, 'reader')).status).toBe(403)
    const accepted = await api(root, `/workspaces/${workspace.id}/homology/jobs`, { project }, 'writer')
    expect(accepted.status).toBe(202)
    const job = accepted.body as ManagedJobRecord
    jobs.push({ root, job })
    await vi.waitFor(() => expect(getManagedJob(root, workspace.id, job.id).status).toBe('running'))
    const stream = await api(root, `/workspaces/${workspace.id}/jobs/${job.id}/events`, undefined, 'writer')
    expect(stream.response.writableFinished).toBe(false)
    saveWorkspace(root, { ...loadWorkspace(root, workspace.id), acl: [] })
    await vi.waitFor(() => expect(stream.response.writableFinished).toBe(true))
    expect((await api(root, `/workspaces/${workspace.id}/jobs/${job.id}`, undefined, 'writer')).status).toBe(404)
    fs.writeFileSync(path.join(directory, job.outputDir, 'finish'), 'done')
    const failed = await terminal(root, job)
    expect(failed.status).toBe('failed')
    expect(failed.outputFile).toBeNull()
    expect(loadWorkspace(root, workspace.id).artifacts.every(artifact => artifact.kind !== 'structure')).toBe(true)
  })
})
