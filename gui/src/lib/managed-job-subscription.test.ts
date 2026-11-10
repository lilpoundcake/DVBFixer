import { afterEach, expect, it, vi } from 'vitest'
import { watchManagedJob } from './managed-job-subscription'
import type { ManagedJobRecord } from './managed-jobs'

const { fetcher, events } = vi.hoisted(() => ({ fetcher: vi.fn(), events: vi.fn() }))
vi.mock('./api-client', () => ({ apiFetch: fetcher, apiSse: events }))
afterEach(() => { vi.useRealTimers(); fetcher.mockReset(); events.mockReset() })
const record = (revision: number): ManagedJobRecord => ({ version: 1, revision, id: 'job', workspaceId: 'workspace',
  command: 'homology', status: 'running', createdAt: '', startedAt: null, finishedAt: null,
  args: [], exitCode: null, stdoutLog: '', stderrLog: '', outputDir: '', outputFile: null })

it('disconnects the viewer without cancelling a persistent job', async () => {
  vi.useFakeTimers()
  events.mockImplementation(async function* (_url, options) {
    yield { data: JSON.stringify(record(3)) }
    await new Promise(resolve => options.signal.addEventListener('abort', resolve, { once: true }))
  })
  const onJob = vi.fn()
  const dispose = watchManagedJob('workspace', 'job', onJob, vi.fn())
  await vi.advanceTimersByTimeAsync(1)
  expect(onJob).toHaveBeenCalledWith(record(3))
  dispose()
  await vi.advanceTimersByTimeAsync(60_000)
  expect(fetcher).not.toHaveBeenCalled() // Especially no DELETE on unmount/reload.
})

it('uses authenticated GET after SSE fails and clears inaccessible job state', async () => {
  vi.useFakeTimers()
  events.mockImplementation(async function* () { yield* []; throw new Error('proxy disconnected') })
  fetcher.mockResolvedValueOnce(Response.json(record(3))).mockResolvedValueOnce(new Response('', { status: 403 }))
  const onJob = vi.fn(), denied = vi.fn()
  const dispose = watchManagedJob('workspace', 'job', onJob, denied)
  await vi.advanceTimersByTimeAsync(3000)
  expect(onJob).toHaveBeenCalledWith(record(3))
  await vi.advanceTimersByTimeAsync(3000)
  expect(denied).toHaveBeenCalledOnce()
  expect(fetcher.mock.calls.every(([, options]) => options.method !== 'DELETE')).toBe(true)
  dispose()
})

it('does not let a stale GET response replace a newer SSE record', async () => {
  vi.useFakeTimers()
  events.mockImplementation(async function* (_url, options) {
    yield { data: JSON.stringify(record(5)) }
    await new Promise(resolve => options.signal.addEventListener('abort', resolve, { once: true }))
  })
  fetcher.mockResolvedValue(Response.json(record(4)))
  const onJob = vi.fn()
  const dispose = watchManagedJob('workspace', 'job', onJob, vi.fn())
  await vi.advanceTimersByTimeAsync(3000)
  expect(onJob).toHaveBeenCalledTimes(1)
  expect(onJob).toHaveBeenCalledWith(record(5))
  dispose()
})

it('honors Retry-After during fallback without resubmitting work', async () => {
  vi.useFakeTimers()
  events.mockImplementation(async function* () { yield* []; throw new Error('no SSE') })
  fetcher.mockResolvedValueOnce(new Response('', { status: 429, headers: { 'Retry-After': '20' } }))
    .mockResolvedValueOnce(Response.json({ ...record(3), status: 'succeeded' }))
  const onJob = vi.fn()
  const dispose = watchManagedJob('workspace', 'job', onJob, vi.fn())
  await vi.advanceTimersByTimeAsync(3000)
  expect(fetcher).toHaveBeenCalledTimes(1)
  await vi.advanceTimersByTimeAsync(19_999)
  expect(fetcher).toHaveBeenCalledTimes(1)
  await vi.advanceTimersByTimeAsync(1)
  expect(onJob).toHaveBeenCalledWith({ ...record(3), status: 'succeeded' })
  expect(fetcher).toHaveBeenCalledTimes(2)
  dispose()
})
