import { apiFetch, apiSse } from './api-client'
import { retryAfterMs } from './job-monitor-polling'
import { acceptManagedJobUpdate, isActiveManagedJob, type ManagedJobRecord } from './managed-jobs'

/** Disposing a viewer disconnects transport only; it NEVER cancels the job. */
export function watchManagedJob(workspaceId: string, jobId: string,
  onJob: (job: ManagedJobRecord) => void, onAccessLost: (message: string) => void,
): () => void {
  const base = `/api/v1/workspaces/${encodeURIComponent(workspaceId)}/jobs/${encodeURIComponent(jobId)}`
  const controller = new AbortController()
  let disposed = false, streaming = false
  let latest: ManagedJobRecord | null = null
  let timer: ReturnType<typeof setTimeout>
  const apply = (job: ManagedJobRecord) => {
    if (disposed || job.id !== jobId || job.workspaceId !== workspaceId || !acceptManagedJobUpdate(latest, job)) return
    latest = job
    onJob(job)
    if (!isActiveManagedJob(job)) { disposed = true; controller.abort(); clearTimeout(timer) }
  }
  void (async () => {
    try {
      for await (const event of apiSse(`${base}/events`, { cache: 'no-store', signal: controller.signal })) {
        try { apply(JSON.parse(event.data) as ManagedJobRecord); streaming = true } catch { /* GET can recover. */ }
      }
    } catch { /* Authenticated GET is the reconnect/proxy fallback. */ }
    finally { streaming = false }
  })()
  const poll = async () => {
    let delay = streaming ? 30_000 : 3000
    try {
      if (typeof document === 'undefined' || !document.hidden) {
        const response = await apiFetch(base, { cache: 'no-store', signal: controller.signal })
        if (response.status === 429) delay = retryAfterMs(response)
        else if (response.status === 403 || response.status === 404 || response.status === 401) {
          disposed = true; controller.abort(); onAccessLost('Homology job is no longer accessible.'); return
        } else if (response.ok) apply(await response.json() as ManagedJobRecord)
      }
    } catch { /* Preserve the last state during temporary network failures. */ }
    if (!disposed) timer = setTimeout(poll, delay)
  }
  timer = setTimeout(poll, 3000)
  return () => { disposed = true; controller.abort(); clearTimeout(timer) }
}
