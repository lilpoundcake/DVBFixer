import { apiFetch } from './api-client'
import type { ManagedJobRecord } from './managed-jobs'

export function retryAfterMs(response: Response, now = Date.now()): number {
  const value = response.headers.get('Retry-After')
  if (!value) return 60_000
  const delay = /^\d+$/.test(value) ? Number(value) * 1000 : Date.parse(value) - now
  return Number.isFinite(delay) && delay >= 0 ? Math.max(3000, delay) : 60_000
}

export interface JobMonitorPage {
  workspaceIds: string[]
  activeJobs: ManagedJobRecord[]
  jobs: ManagedJobRecord[]
  nextCursor: string | null
}

export class MonitorHttpError extends Error {
  readonly status: number
  readonly retryMs: number
  constructor(status: number, retryMs: number) {
    super(`HTTP ${status}`)
    this.status = status
    this.retryMs = retryMs
  }
}

export async function fetchJobMonitorPage(cursor: string | null, request = apiFetch, signal?: AbortSignal): Promise<JobMonitorPage> {
  const response = await request(`/api/v1/jobs?limit=30${cursor ? `&before=${encodeURIComponent(cursor)}` : ''}`, { signal, cache: 'no-store' })
  if (!response.ok) throw new MonitorHttpError(response.status, response.status === 429 ? retryAfterMs(response) : 5000)
  return response.json() as Promise<JobMonitorPage>
}

/** Refresh only active work and the newest page, never prefetch old history.
 * One read at a time, at most 12/minute regardless of workspace count. */
export function pollJobMonitor(
  receive: (page: JobMonitorPage | null, error: string | null) => void,
  options: { request?: typeof apiFetch; paused?: () => boolean } = {},
): () => void {
  const request = options.request || apiFetch
  const paused = options.paused || (() => typeof document !== 'undefined' && document.hidden)
  const controller = new AbortController()
  let timer: ReturnType<typeof setTimeout>
  let failures = 0
  const poll = async () => {
    if (controller.signal.aborted) return
    let delay = 5000
    if (!paused()) {
      try {
        const page = await fetchJobMonitorPage(null, request, controller.signal)
        if (controller.signal.aborted) return
        receive(page, null)
        failures = 0
      } catch (error) {
        if (controller.signal.aborted) return
        if (error instanceof MonitorHttpError && error.status === 429) {
          delay = error.retryMs
          receive(null, `API rate limit reached; monitoring will retry in ${Math.ceil(delay / 1000)}s.`)
        } else {
          delay = Math.min(60_000, 5000 * 2 ** ++failures)
          const denied = error instanceof MonitorHttpError && [401, 403].includes(error.status)
          receive(denied ? { workspaceIds: [], activeJobs: [], jobs: [], nextCursor: null } : null, (error as Error).message)
        }
      }
    }
    if (!controller.signal.aborted) timer = setTimeout(poll, delay)
  }
  void poll()
  return () => { controller.abort(); clearTimeout(timer) }
}
