import { describe, expect, it } from 'vitest'
import type { WorkspaceArtifact } from '../stores/workspaceStore'
import { selectWorkspaceContextRow, selectWorkspaceRows, workspaceRowArtifactIds, workspaceRunRows } from './workspace-run-groups'

function artifact(file: string, command = 'prepare'): WorkspaceArtifact {
  return { id: file, file, command, name: file.split('/').pop()!, kind: 'artifact' }
}

describe('Workspace run groups', () => {
  const files = [artifact('input.pdb'), artifact('runs/job_one/result.pdb'), artifact('runs/job_one/diagnostics.tar.gz'),
    artifact('runs/homology_two/result/model.pdb', 'homology'), artifact('runs/naming_three/output.pdb', 'atom-names')]
  it('shows one collapsed row per run for every command and keeps imported files independent', () => {
    const rows = workspaceRunRows(files, new Set())
    expect(rows.map(row => row.type)).toEqual(['file', 'run', 'run', 'run'])
    expect(rows[1]).toMatchObject({ id: 'runs/job_one', command: 'prepare', files: files.slice(1, 3) })
  })
  it('expands only the requested run without changing paths, ordering or input data', () => {
    const before = JSON.stringify(files)
    const rows = workspaceRunRows(files, new Set(['runs/job_one']))
    expect(rows.flatMap(row => row.type === 'file' ? [row.artifact.file] : [])).toEqual(['input.pdb', files[1].file, files[2].file])
    expect(JSON.stringify(files)).toBe(before)
  })
  it('presents rejected, failed and cancelled runs without conflating them with success', () => {
    for (const status of ['validation_failed', 'failed', 'cancelled'] as const) {
      expect(workspaceRunRows([{ ...files[1], run: { command: 'model', status } }], new Set())[0]).toMatchObject({ status })
    }
    expect(workspaceRunRows([{ ...files[1], name: '[rejected] model.pdb' }], new Set())[0]).toMatchObject({ status: 'validation_failed' })
  })

  it.each(['ctrlKey', 'metaKey'] as const)('selects and deselects entire collapsed groups with %s', key => {
    const rows = workspaceRunRows(files, new Set())
    let state = selectWorkspaceRows(rows, 1, new Set(), null)
    expect([...state.selected]).toEqual([files[1].id, files[2].id])
    state = selectWorkspaceRows(rows, 2, state.selected, state.anchor, { [key]: true })
    expect([...state.selected]).toEqual([files[1].id, files[2].id, files[3].id])
    state = selectWorkspaceRows(rows, 1, state.selected, state.anchor, { [key]: true })
    expect([...state.selected]).toEqual([files[3].id])
  })

  it('supports Shift ranges across imported files and multiple collapsed groups', () => {
    const rows = workspaceRunRows(files, new Set())
    const anchor = selectWorkspaceRows(rows, 0, new Set(), null)
    const range = selectWorkspaceRows(rows, 2, anchor.selected, anchor.anchor, { shiftKey: true })
    expect([...range.selected]).toEqual(files.slice(0, 4).map(file => file.id))
    expect(range.anchor).toBe(anchor.anchor)
  })

  it('does not change a group anchor when expansion inserts rows before the next group', () => {
    const collapsed = workspaceRunRows(files, new Set())
    const anchor = selectWorkspaceRows(collapsed, 1, new Set(), null)
    const expanded = workspaceRunRows(files, new Set(['runs/job_one']))
    const nextGroup = expanded.findIndex(row => row.type === 'run' && row.id === 'runs/homology_two')
    const selection = selectWorkspaceRows(expanded, nextGroup, anchor.selected, anchor.anchor, { shiftKey: true })
    expect([...selection.selected]).toEqual([files[1].id, files[2].id, files[3].id])
  })

  it('keeps file selection and its anchor valid after the containing group collapses', () => {
    const expanded = workspaceRunRows(files, new Set(['runs/job_one']))
    const file = expanded.findIndex(row => row.type === 'file' && row.artifact.id === files[2].id)
    const anchor = selectWorkspaceRows(expanded, file, new Set(), null)
    const collapsed = workspaceRunRows(files, new Set())
    const next = selectWorkspaceRows(collapsed, 2, anchor.selected, anchor.anchor, { shiftKey: true })
    expect(next.selected).toEqual(new Set([files[1].id, files[2].id, files[3].id]))
  })

  it('selects a range of individual files inside an expanded group without selecting unrelated runs', () => {
    const rows = workspaceRunRows(files, new Set(['runs/job_one']))
    const first = rows.findIndex(row => row.type === 'file' && row.artifact.id === files[1].id)
    const anchor = selectWorkspaceRows(rows, first, new Set(), null)
    const next = selectWorkspaceRows(rows, first + 1, anchor.selected, anchor.anchor, { shiftKey: true })
    expect(next.selected).toEqual(new Set([files[1].id, files[2].id]))
  })

  it('preserves multi-selection on right-click, including partially selected groups', () => {
    const rows = workspaceRunRows(files, new Set())
    const selected = new Set([files[0].id, files[1].id, files[3].id])
    const context = selectWorkspaceContextRow(rows, 1, selected, 'file:input.pdb')
    expect(context.selected).toEqual(selected)
    expect(context.anchor).toBe('file:input.pdb')
    expect(selectWorkspaceContextRow(rows, 3, selected, null).selected).toEqual(new Set([files[4].id]))
  })

  it('includes collapsed group members in the visible action set, without duplicate expanded files', () => {
    const rows = workspaceRunRows(files, new Set(['runs/job_one']))
    expect([...new Set(rows.flatMap(workspaceRowArtifactIds))]).toEqual(files.map(file => file.id))
  })

  it('does not reuse stale numeric positions when filtering removes the anchor', () => {
    const rows = workspaceRunRows(files.slice(1), new Set())
    const selected = selectWorkspaceRows(rows, 1, new Set(), 'file:missing', { shiftKey: true })
    expect(selected.selected).toEqual(new Set([files[3].id]))
  })
})
