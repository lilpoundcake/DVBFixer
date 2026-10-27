import { apiFetch } from './api-client'
import { useAltlocStore, EMPTY_ALTLOC_SESSION } from '../stores/altlocStore'
import { useWorkspaceStore } from '../stores/workspaceStore'
import { altlocSelectionComplete } from './altloc-selection'
import type { AltlocRequest, AltlocResponse } from './altloc-types'

export function openAltlocReview(workspaceId: string): void {
  if (useWorkspaceStore.getState().active?.id !== workspaceId) return
  const session = useAltlocStore.getState().sessions[workspaceId]
  const key = session?.activeKey || session?.report?.residues[0]?.key
  if (key) useAltlocStore.getState().focus(workspaceId, key)
  window.dispatchEvent(new Event('dvbfixer:open-altlocs'))
}

/** An operation survives switching dock tabs; all state remains workspace-scoped. */
export async function runAltlocOperation(workspaceId: string, action: 'analyze' | 'apply' | 'run'): Promise<void> {
  const store = useAltlocStore.getState()
  const initial = store.sessions[workspaceId] || EMPTY_ALTLOC_SESSION
  if (!initial.input || initial.busy) return
  const generation = initial.generation + 1
  store.patch(workspaceId, { generation, busy: true, error: '', result: '', reportFile: '' })
  const current = () => useAltlocStore.getState().sessions[workspaceId]?.generation === generation
  const request = async (body: AltlocRequest): Promise<AltlocResponse> => {
    const response = await apiFetch(`/api/workspaces/${encodeURIComponent(workspaceId)}/altlocs`, {
      method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
    })
    const data = await response.json()
    if (!response.ok) throw new Error(typeof data.error === 'string' ? data.error : 'Altloc operation failed')
    return data as AltlocResponse
  }
  const options = { inputArtifactId: initial.input, policy: initial.policy, fallback: initial.fallback,
    ...(initial.policy === 'prefer' ? { prefer: initial.prefer } : {}) }
  try {
    let report = initial.report
    let choices = initial.choices
    if (action !== 'apply') {
      const data = await request({ ...options, action: 'analyze' })
      if (!current()) return
      report = data.report
      choices = Object.fromEntries(report.residues.filter(row => row.selected).map(row => [row.key, row.selected!]))
      store.patch(workspaceId, { report, choices, activeKey: report.residues[0]?.key || '', focusLabel: '' })
      if (action === 'analyze' || !altlocSelectionComplete(report.residues, choices)) {
        if (initial.policy === 'manual' || action === 'run') openAltlocReview(workspaceId)
        if (action === 'run') throw new Error('Some states require manual selection in the Altlocs table before Apply.')
        return
      }
    }
    if (!report || !altlocSelectionComplete(report.residues, choices)) throw new Error('Select one state per residue in the Altlocs table.')
    const data = await request({ ...options, action: 'apply', plan: { schemaVersion: 1, sourceSha256: report.sourceSha256, selections: choices } })
    if (!current()) return
    store.patch(workspaceId, { result: `Created ${data.outputArtifact?.name}. Original input unchanged.`, reportFile: data.reportFile || '' })
    if (useWorkspaceStore.getState().active?.id === workspaceId) await useWorkspaceStore.getState().reload()
  } catch (cause) {
    if (current()) store.patch(workspaceId, { error: cause instanceof Error ? cause.message : 'Altloc operation failed' })
  } finally { if (current()) store.patch(workspaceId, { busy: false }) }
}
