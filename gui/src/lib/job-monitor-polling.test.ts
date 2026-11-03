import { afterEach, expect, it, vi } from 'vitest'
import { fetchJobMonitorPage, pollJobMonitor, retryAfterMs } from './job-monitor-polling'

const page = JSON.stringify({ workspaceIds: ['one'], activeJobs: [], jobs: [], nextCursor: 'older' })

afterEach(() => { vi.useRealTimers(); vi.restoreAllMocks() })

it('loads the newest page immediately and bounds refreshes without prefetching older pages', async () => {
  vi.useFakeTimers()
  const request = vi.fn(async (_url: RequestInfo | URL) => new Response(page))
  const receive = vi.fn()
  const stop = pollJobMonitor(receive, { request, paused: () => false })
  await vi.advanceTimersByTimeAsync(59_000)
  expect(request).toHaveBeenCalledTimes(12)
  expect(request.mock.calls.every(([url]) => String(url) === '/api/v1/jobs?limit=30')).toBe(true)
  stop()
  await vi.advanceTimersByTimeAsync(60_000)
  expect(request).toHaveBeenCalledTimes(12)
  await fetchJobMonitorPage('older', request)
  expect(request.mock.calls.at(-1)?.[0]).toBe('/api/v1/jobs?limit=30&before=older')
})

it('honors Retry-After without discarding known jobs or retrying other workspaces during cooldown', async () => {
  vi.useFakeTimers()
  const request = vi.fn().mockResolvedValueOnce(new Response('', { status: 429, headers: { 'Retry-After': '60' } }))
    .mockResolvedValue(new Response(page))
  const receive = vi.fn()
  const stop = pollJobMonitor(receive, { request, paused: () => false })
  await vi.advanceTimersByTimeAsync(59_999)
  expect(request).toHaveBeenCalledTimes(1)
  expect(receive.mock.calls[0][0]).toBeNull()
  await vi.advanceTimersByTimeAsync(1)
  expect(request).toHaveBeenCalledTimes(2)
  expect(request.mock.calls[1][0]).toBe('/api/v1/jobs?limit=30')
  stop()
})

it('does not poll a hidden browser tab and does not apply a late response after disposal', async () => {
  vi.useFakeTimers()
  const request = vi.fn(async () => new Response('[]'))
  const receive = vi.fn()
  const stop = pollJobMonitor(receive, { request, paused: () => true })
  await vi.advanceTimersByTimeAsync(60_000)
  expect(request).not.toHaveBeenCalled()
  stop()
  let resolve!: (response: Response) => void
  const cancel = pollJobMonitor(receive, {
    request: () => new Promise<Response>(done => { resolve = done }), paused: () => false,
  })
  cancel()
  resolve(new Response('[]'))
  await vi.advanceTimersByTimeAsync(0)
  expect(receive).not.toHaveBeenCalled()
})

it('accepts HTTP-date Retry-After and safely falls back for missing/invalid values', () => {
  expect(retryAfterMs(new Response('', { headers: { 'Retry-After': 'Thu, 01 Jan 1970 00:01:00 GMT' } }), 0)).toBe(60_000)
  expect(retryAfterMs(new Response('', { headers: { 'Retry-After': 'bad' } }))).toBe(60_000)
})
