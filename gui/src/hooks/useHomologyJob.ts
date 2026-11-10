import { useEffect, useRef, useState } from 'react'
import { apiFetch } from '../lib/api-client'
import { acceptManagedJobUpdate, isActiveManagedJob, type ManagedJobRecord } from '../lib/managed-jobs'
import { selectHomologyJob } from '../lib/homology-jobs'
import { watchManagedJob } from '../lib/managed-job-subscription'

/** Shared persisted jobs, not a GUI scheduler or an HTTP-request-owned process. */
export function useHomologyJob(workspaceId: string, projectId: string) {
  const [record, setRecord] = useState<ManagedJobRecord | null>(null)
  const [error, setError] = useState('')
  const [submitting, setSubmitting] = useState(false)
  const [stopping, setStopping] = useState(false)
  const generation = useRef(0)
  const mounted = useRef(true)
  const selection = useRef({ workspaceId, projectId })
  selection.current = { workspaceId, projectId }
  const job = record?.workspaceId === workspaceId && record.workflow?.projectId === projectId ? record : null
  const activeId = job && isActiveManagedJob(job) ? job.id : null
  useEffect(() => { mounted.current = true; return () => { mounted.current = false } }, [])

  useEffect(() => {
    const request = ++generation.current
    const controller = new AbortController()
    setRecord(null); setError(''); setSubmitting(false); setStopping(false)
    if (workspaceId && projectId) void apiFetch(`/api/v1/workspaces/${encodeURIComponent(workspaceId)}/jobs`, {
      cache: 'no-store', signal: controller.signal,
    }).then(async response => {
      if (!response.ok) throw new Error(`Cannot restore Homology jobs: HTTP ${response.status}`)
      const jobs = await response.json() as ManagedJobRecord[]
      if (!controller.signal.aborted && request === generation.current) setRecord(selectHomologyJob(jobs, workspaceId, projectId))
    }).catch(reason => { if (!controller.signal.aborted && request === generation.current) setError(String(reason.message || reason)) })
    return () => { controller.abort() }
  }, [workspaceId, projectId])

  useEffect(() => {
    if (!activeId) return
    const selected = () => selection.current.workspaceId === workspaceId && selection.current.projectId === projectId
    const apply = (next: ManagedJobRecord) => {
      if (selected()) setRecord(current => current?.id === activeId && acceptManagedJobUpdate(current, next) ? next : current)
    }
    return watchManagedJob(workspaceId, activeId, apply, message => { if (selected()) { setRecord(null); setError(message) } })
  }, [workspaceId, projectId, activeId])

  useEffect(() => { if (!activeId) setStopping(false) }, [activeId])

  const start = async (project: unknown) => {
    const request = ++generation.current
    const selected = { workspaceId, projectId }
    const current = () => mounted.current && generation.current === request && selection.current.workspaceId === selected.workspaceId && selection.current.projectId === selected.projectId
    setSubmitting(true); setError('')
    try {
      const response = await apiFetch(`/api/v1/workspaces/${encodeURIComponent(workspaceId)}/homology/jobs`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ project }),
      })
      const body = await response.json() as ManagedJobRecord & { error?: { message?: string } }
      if (!current()) return
      if (!response.ok) throw new Error(body.error?.message || `Homology submit failed: HTTP ${response.status}`)
      setRecord(body)
    } catch (reason) { if (current()) setError((reason as Error).message) }
    finally { if (current()) setSubmitting(false) }
  }

  const stop = async () => {
    if (!job || !isActiveManagedJob(job) || stopping) return
    const target = job
    setStopping(true); setError('')
    try {
      const response = await apiFetch(`/api/v1/workspaces/${encodeURIComponent(target.workspaceId)}/jobs/${encodeURIComponent(target.id)}`, { method: 'DELETE' })
      if (!response.ok) throw new Error(`Cannot stop Homology: HTTP ${response.status}`)
    } catch (reason) {
      if (selection.current.workspaceId === target.workspaceId && selection.current.projectId === target.workflow?.projectId) {
        setStopping(false); setError((reason as Error).message)
      }
    }
  }
  return { job, error, submitting, stopping: stopping || Boolean(activeId && job?.cancellationRequestedAt), start, stop, active: Boolean(activeId) }
}
