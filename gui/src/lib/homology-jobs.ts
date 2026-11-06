import { isActiveManagedJob, type ManagedJobRecord } from './managed-jobs'

/** A workspace can contain several Homology projects; never borrow another run. */
export function selectHomologyJob(jobs: readonly ManagedJobRecord[], workspaceId: string, projectId: string): ManagedJobRecord | null {
  return jobs.filter(job => job.workspaceId === workspaceId && job.command === 'homology' &&
    job.workflow?.kind === 'homology' && job.workflow.projectId === projectId)
    .sort((a, b) => Number(isActiveManagedJob(b)) - Number(isActiveManagedJob(a)) || b.createdAt.localeCompare(a.createdAt))[0] || null
}
