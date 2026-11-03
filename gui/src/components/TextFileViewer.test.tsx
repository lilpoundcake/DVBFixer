import { renderToStaticMarkup } from 'react-dom/server'
import { useStore } from 'zustand'
import { afterEach, expect, it, vi } from 'vitest'
import { TextFileViewer } from './TextFileViewer'
import { useWorkspaceStore } from '../stores/workspaceStore'

vi.mock('../stores/workspaceStore', async importOriginal => {
  const original = await importOriginal<typeof import('../stores/workspaceStore')>()
  const store = original.useWorkspaceStore
  return { ...original, useWorkspaceStore: Object.assign(
    (selector: Parameters<typeof store>[0]) => useStore({ ...store, getInitialState: store.getState }, selector), store,
  ) }
})
afterEach(() => { useWorkspaceStore.setState({ textPreview: null }) })

it('wraps long text to the docking panel width while preserving whitespace', () => {
  useWorkspaceStore.setState({ textPreview: { workspaceId: 'w', file: 'stdout.log', name: 'Log' } })
  const html = renderToStaticMarkup(<TextFileViewer />)
  expect(html).toContain('white-space:pre-wrap')
  expect(html).toContain('overflow-wrap:anywhere')
  expect(html).toContain('min-width:0')
  expect(html).toContain('max-width:100%')
})
