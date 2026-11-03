import { useEffect, useRef, useState } from 'react'
import Alert from '@mui/material/Alert'
import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import LinearProgress from '@mui/material/LinearProgress'
import Typography from '@mui/material/Typography'
import { apiFetch } from '../lib/api-client'
import { formatJobTimestamp, isActiveManagedJob, managedJobStatusLabel, type ManagedJobRecord } from '../lib/managed-jobs'
import { useWorkspaceStore } from '../stores/workspaceStore'
import { panelToolbarSx, workspaceRowColors } from '../ui/layout'
import { fetchJobMonitorPage, MonitorHttpError, pollJobMonitor, type JobMonitorPage } from '../lib/job-monitor-polling'

const jobKey = (job: ManagedJobRecord) => `${job.workspaceId}:${job.id}`

export function JobMonitor() {
  const workspaces = useWorkspaceStore(state => state.workspaces)
  const setTextPreview = useWorkspaceStore(state => state.setTextPreview)
  const [jobs, setJobs] = useState<ManagedJobRecord[]>([])
  const [error, setError] = useState('')
  const [loading, setLoading] = useState(true)
  const [loadingOlder, setLoadingOlder] = useState(false)
  const [nextCursor, setNextCursor] = useState<string | null>(null)
  const [olderRetryAt, setOlderRetryAt] = useState(0)
  const loadedOlder = useRef(false)
  const historyRequest = useRef<AbortController | null>(null)
  const [stopping, setStopping] = useState<Set<string>>(new Set())
  const [cancelError, setCancelError] = useState('')
  const [now, setNow] = useState(Date.now())

  const mergePage = (page: JobMonitorPage, newest: boolean) => {
    setJobs(current => {
      const completed = new Map(current.filter(job => !isActiveManagedJob(job) && page.workspaceIds.includes(job.workspaceId)).map(job => [jobKey(job), job]))
      for (const job of page.jobs) completed.set(jobKey(job), job)
      const active = (newest ? page.activeJobs : current.filter(isActiveManagedJob))
        .filter(job => page.workspaceIds.includes(job.workspaceId) && !completed.has(jobKey(job)))
      return [...active, ...completed.values()].sort((a, b) =>
        Number(isActiveManagedJob(b)) - Number(isActiveManagedJob(a)) ||
        (b.finishedAt || b.createdAt).localeCompare(a.finishedAt || a.createdAt) || jobKey(b).localeCompare(jobKey(a)))
    })
    if (!newest || !loadedOlder.current) setNextCursor(page.nextCursor)
    if (newest) setStopping(current => new Set([...current].filter(key => page.activeJobs.some(job => jobKey(job) === key))))
    setLoading(false)
    setNow(Date.now())
  }

  useEffect(() => {
    const stop = pollJobMonitor((page, message) => {
      if (page) {
        if (message) { setJobs([]); loadedOlder.current = false; setNextCursor(null) }
        else mergePage(page, true)
      }
      setError(message || '')
      setLoading(false)
    })
    return () => { stop(); historyRequest.current?.abort() }
  }, [])

  const loadOlder = async () => {
    if (!nextCursor || loadingOlder || historyRequest.current || Date.now() < olderRetryAt) return
    setLoadingOlder(true)
    const controller = new AbortController()
    historyRequest.current = controller
    try {
      const page = await fetchJobMonitorPage(nextCursor, undefined, controller.signal)
      if (controller.signal.aborted) return
      loadedOlder.current = true
      mergePage(page, false)
      setError('')
    } catch (reason) {
      if (controller.signal.aborted) return
      if (reason instanceof MonitorHttpError && reason.status === 429) setOlderRetryAt(Date.now() + reason.retryMs)
      setError((reason as Error).message)
    } finally {
      if (historyRequest.current === controller) historyRequest.current = null
      if (!controller.signal.aborted) setLoadingOlder(false)
    }
  }

  const stop = async (job: ManagedJobRecord) => {
    const key = jobKey(job)
    if (stopping.has(key)) return
    setStopping(current => new Set([...current, key]))
    setCancelError('')
    try {
      const response = await apiFetch(`/api/v1/workspaces/${encodeURIComponent(job.workspaceId)}/jobs/${encodeURIComponent(job.id)}`, { method: 'DELETE' })
      if (!response.ok) throw new Error(`HTTP ${response.status}`)
      // Keep "Stopping" until persisted state confirms a terminal outcome.
    } catch (error) {
      setCancelError(`${job.command}: ${(error as Error).message}`)
      setStopping(current => new Set([...current].filter(item => item !== key)))
    }
  }

  return <Box sx={{ display: 'flex', flexDirection: 'column', height: '100%', minHeight: 0, minWidth: 0 }}>
    <Box component="header" sx={{ ...panelToolbarSx, flexDirection: 'column', alignItems: 'flex-start', gap: 0.5 }}>
      <Typography variant="subtitle2">Job Monitor</Typography>
      <Typography variant="caption" color="text.secondary">
        Managed tasks across accessible workspaces. Progress is indeterminate when the runner supplies no work counts.
        Completion is not scientific acceptance. Existing scheduling and workspace limits still apply.
      </Typography>
    </Box>
    <Box component="section" aria-label="Tasks" sx={{ p: 1, flex: 1, minHeight: 0, overflow: 'auto' }}>
    <Typography variant="subtitle2" sx={{ mb: 1 }}>Tasks · {jobs.length} loaded</Typography>
    <Typography variant="caption" color="text.secondary">Active tasks first; completed tasks newest first. Older history loads on demand.</Typography>
    {error && <Alert severity="warning">{error}</Alert>}
    {cancelError && <Alert severity="error">{cancelError}</Alert>}
    {loading && <Typography variant="body2">Loading tasks…</Typography>}
    {!loading && !jobs.length && <Typography variant="body2" sx={{ mt: 1 }}>No managed jobs available.</Typography>}
    {jobs.map((job, index) => {
      const workspace = workspaces.find(item => item.id === job.workspaceId)
      const elapsed = Math.max(0, Math.floor(((job.finishedAt ? Date.parse(job.finishedAt) : now) - Date.parse(job.startedAt || job.createdAt)) / 1000))
      return <Box key={jobKey(job)} sx={{ px: 1.5, py: 1, borderBottom: 1, borderColor: 'divider', overflowWrap: 'anywhere',
        bgcolor: index % 2 ? workspaceRowColors.stripe : workspaceRowColors.base,
        '&:hover': { bgcolor: workspaceRowColors.hover },
      }}>
        <Typography variant="body2">{managedJobStatusLabel(job)}</Typography>
        <Typography variant="caption" display="block">{workspace?.name || job.workspaceId} · {job.id}</Typography>
        <Typography variant="caption" display="block" title={job.startedAt || undefined}>
          Started (local): {formatJobTimestamp(job.startedAt)} · Elapsed: {elapsed}s
        </Typography>
        {isActiveManagedJob(job) && <LinearProgress aria-label={`${job.command} progress (indeterminate)`} sx={{ my: 1 }} />}
        {job.error && <Typography variant="caption" color="error">{job.error}</Typography>}
        <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1, mt: 0.5 }}>
          {isActiveManagedJob(job) && <Button variant="outlined" color="error"
            disabled={stopping.has(jobKey(job)) || !workspace || workspace.effectiveRole === 'reader'}
            onClick={() => void stop(job)}>{stopping.has(jobKey(job)) ? 'Stopping…' : 'Stop'}</Button>}
          {[['Stdout', job.stdoutLog], ['Stderr', job.stderrLog]] .map(([name, file]) => <Button key={name}
            onClick={() => {
              setTextPreview({ workspaceId: job.workspaceId, file, name: `${job.command} · ${name}` })
              window.dispatchEvent(new Event('dvbfixer:open-text-viewer'))
            }}>{name}</Button>)}
        </Box>
      </Box>
    })}
    {nextCursor && <Button sx={{ mt: 1 }} variant="outlined" disabled={loadingOlder || now < olderRetryAt}
      onClick={() => void loadOlder()}>{loadingOlder ? 'Loading…' : 'Load older tasks'}</Button>}
    </Box>
  </Box>
}
