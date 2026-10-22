import { renderToStaticMarkup } from 'react-dom/server'
import { useStore } from 'zustand'
import { afterEach, expect, it, vi } from 'vitest'
import { AltlocTable } from './AltlocTable'
import { useWorkspaceStore, type WorkspaceManifest } from '../stores/workspaceStore'
import { useAltlocStore } from '../stores/altlocStore'
import type { AltlocReport } from '../lib/altloc-types'

// Render client snapshots through SSR without changing the production stores.
vi.mock('../stores/workspaceStore', async importOriginal => {
  const original = await importOriginal<typeof import('../stores/workspaceStore')>()
  const store = original.useWorkspaceStore
  return { ...original, useWorkspaceStore: Object.assign(
    (selector: Parameters<typeof store>[0]) => useStore({ ...store, getInitialState: store.getState }, selector), store,
  ) }
})
vi.mock('../stores/altlocStore', async importOriginal => {
  const original = await importOriginal<typeof import('../stores/altlocStore')>()
  const store = original.useAltlocStore
  return { ...original, useAltlocStore: Object.assign(
    (selector: Parameters<typeof store>[0]) => useStore({ ...store, getInitialState: store.getState }, selector), store,
  ) }
})
afterEach(() => { vi.restoreAllMocks(); useWorkspaceStore.setState({ active: null }); useAltlocStore.setState({ sessions: {} }) })
it('contains review controls and the selection table, not a CLI launcher or embedded viewer', () => {
  const html = renderToStaticMarkup(<AltlocTable />)
  expect(html).toContain('<table')
  expect(html).toContain('Alternate-state selection')
  expect(html).toContain('Review states')
  expect(html).toContain('3D display')
  expect(html).not.toContain('<canvas')
  expect(html).not.toContain('Apply selection')
  expect(html).not.toContain('Preferred label')
  expect(html).not.toContain('Mol* preview')
})
it('renders one checked checkbox per residue and keeps drafts outside the launcher lifecycle', () => {
  useWorkspaceStore.setState({ active: { id: 'w' } as WorkspaceManifest })
  const row = { key: 'key', model: 1, chain: 'C', residueNumber: '82', insertionCode: 'A', occurrence: 1,
    selected: 'A', method: 'occupancy', sharedSerials: [], warnings: [], states: ['A', 'B'].map(label => ({ label, residueName: 'SER', occupancy: .5, serials: [1], atomNames: ['CB'], atomCount: 1 })) }
  const report: AltlocReport = { schemaVersion: 1, sourceSha256: 'source', normalizedSha256: 'normalized',
    policy: 'occupancy', prefer: null, fallback: 'error', overrides: {}, residues: [row], warnings: [],
    unresolved: [], inputAtoms: 2, removedAtoms: null, outputSha256: null, outputWritten: false,
    toolVersion: '0.9.1', chainMap: {} }
  useAltlocStore.getState().patch('w', { report, choices: { key: 'A' } })
  const first = renderToStaticMarkup(<AltlocTable />)
  expect(first.match(/checked=""/g)).toHaveLength(1)
  useAltlocStore.getState().choose('w', 'key', 'B')
  const second = renderToStaticMarkup(<AltlocTable />)
  expect(second.match(/checked=""/g)).toHaveLength(1)
  expect(useAltlocStore.getState().sessions.w.choices).toEqual({ key: 'B' })
  useAltlocStore.getState().focus('w', 'key', 'A')
  expect(useAltlocStore.getState().sessions.w.focusLabel).toBe('A')
  expect(useAltlocStore.getState().sessions.w.choices.key).toBe('B')
  useAltlocStore.getState().choose('w', 'key', 'A')
  useAltlocStore.getState().focus('w', 'key', 'B')
  useAltlocStore.getState().choose('w', 'key', 'A')
  expect(useAltlocStore.getState().sessions.w.focusLabel).toBe('B')
  expect(useAltlocStore.getState().sessions.w.choices.key).toBe('A')
  useWorkspaceStore.setState({ active: { id: 'other' } as WorkspaceManifest })
  expect(renderToStaticMarkup(<AltlocTable />)).not.toContain('82')
})
