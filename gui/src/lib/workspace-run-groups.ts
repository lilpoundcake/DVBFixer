import type { WorkspaceArtifact } from '../stores/workspaceStore'

export type WorkspaceRunRow =
  | { type: 'file'; artifact: WorkspaceArtifact }
  | { type: 'run'; id: string; command: string; files: WorkspaceArtifact[]; status?: string }

export function workspaceRowArtifactIds(row: WorkspaceRunRow): string[] {
  return row.type === 'file' ? [row.artifact.id] : row.files.map(file => file.id)
}

function rowKey(row: WorkspaceRunRow): string {
  return row.type === 'file' ? `file:${row.artifact.id}` : `run:${row.id}`
}

/** Native file-manager selection, including collapsed run groups. Anchors are
 * identities, not array indices that change when a group expands or collapses. */
export function selectWorkspaceRows(
  rows: readonly WorkspaceRunRow[], index: number, current: ReadonlySet<string>, anchor: string | null,
  modifiers: { shiftKey?: boolean; ctrlKey?: boolean; metaKey?: boolean } = {},
): { selected: Set<string>; anchor: string | null } {
  const row = rows[index]
  if (!row) return { selected: new Set(current), anchor }
  const ids = workspaceRowArtifactIds(row)
  if (modifiers.shiftKey && anchor !== null) {
    let start = rows.findIndex(item => rowKey(item) === anchor)
    // A previously selected file can now be inside a collapsed group.
    if (start < 0 && anchor.startsWith('file:')) {
      start = rows.findIndex(item => item.type === 'run' && item.files.some(file => `file:${file.id}` === anchor))
    }
    if (start >= 0) {
      const selected = new Set(current)
      for (const item of rows.slice(Math.min(start, index), Math.max(start, index) + 1)) {
        for (const id of workspaceRowArtifactIds(item)) selected.add(id)
      }
      return { selected, anchor }
    }
  }
  if (modifiers.ctrlKey || modifiers.metaKey) {
    const selected = new Set(current)
    const remove = ids.every(id => selected.has(id))
    for (const id of ids) { if (remove) selected.delete(id); else selected.add(id) }
    return { selected, anchor: rowKey(row) }
  }
  return { selected: new Set(ids), anchor: rowKey(row) }
}

export function selectWorkspaceContextRow(
  rows: readonly WorkspaceRunRow[], index: number, current: ReadonlySet<string>, anchor: string | null,
): { selected: Set<string>; anchor: string | null } {
  const row = rows[index]
  if (!row || workspaceRowArtifactIds(row).some(id => current.has(id))) return { selected: new Set(current), anchor }
  return selectWorkspaceRows(rows, index, current, anchor)
}

/** Group presentation only: never rewrite historical manifests or file paths. */
export function workspaceRunRows(artifacts: readonly WorkspaceArtifact[], expanded: ReadonlySet<string>): WorkspaceRunRow[] {
  const groups = new Map<string, Extract<WorkspaceRunRow, { type: 'run' }>>()
  const rows: WorkspaceRunRow[] = []
  for (const artifact of artifacts) {
    const parts = artifact.file.split('/')
    if (parts[0] !== 'runs' || parts.length < 3) {
      rows.push({ type: 'file', artifact })
      continue
    }
    const id = parts.slice(0, 2).join('/')
    let group = groups.get(id)
    if (!group) {
      group = { type: 'run', id, command: artifact.run?.command || artifact.command || 'Run', files: [],
        status: artifact.run?.status }
      groups.set(id, group)
      rows.push(group)
    }
    group.files.push(artifact)
    if (artifact.run?.status) group.status = artifact.run.status
    else if (artifact.name.startsWith('[rejected]')) group.status = 'validation_failed'
  }
  return rows.flatMap(row => row.type === 'run' && expanded.has(row.id)
    ? [row, ...row.files.map(artifact => ({ type: 'file' as const, artifact }))] : [row])
}
