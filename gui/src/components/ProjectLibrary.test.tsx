import { renderToStaticMarkup } from 'react-dom/server'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import { ProjectLibrary } from './ProjectLibrary'
import type { WorkspaceArtifact } from '../stores/workspaceStore'

const { workspaceState } = vi.hoisted(() => ({
  workspaceState: { active: null as null | { id: string; name: string; artifacts: WorkspaceArtifact[] }, workspaces: [], loading: false },
}))

vi.mock('../stores/workspaceStore', () => ({
  useWorkspaceStore: (selector: (state: typeof workspaceState) => unknown) => selector(workspaceState),
  workspaceFileUrl: vi.fn(),
}))
vi.mock('../stores/structureStore', () => ({
  useStructureStore: (selector: (state: object) => unknown) => selector({}),
}))
vi.mock('../stores/selectionStore', () => ({
  useSelectionStore: (selector: (state: object) => unknown) => selector({}),
}))

describe('Workspace import control', () => {
  beforeEach(() => {
    workspaceState.active = { id: 'workspace-a', name: 'My project', artifacts: [] }
    workspaceState.loading = false
  })

  it('places a labelled import button after the viewer selectors', () => {
    const markup = renderToStaticMarkup(<ProjectLibrary mode="workspace" />)
    expect(markup).toContain('Import files')
    expect(markup.indexOf('Import files</button>')).toBeGreaterThan(markup.indexOf('value="secondary"'))
    expect(markup).toContain('type="file"')
    expect(markup).toContain('multiple=""')
  })

  it.each(['no workspace', 'busy'])('disables import when %s', reason => {
    if (reason === 'no workspace') workspaceState.active = null
    else workspaceState.loading = true
    const markup = renderToStaticMarkup(<ProjectLibrary mode="workspace" />)
    expect(markup).toMatch(/<button[^>]*disabled=""[^>]*>[^]*?Import files<\/button>/)
  })

  it('keeps file import out of the workspace switcher', () => {
    expect(renderToStaticMarkup(<ProjectLibrary />)).not.toContain('Import files</button>')
  })

  it('explains that A/B selects the structure viewer, not the workspace', () => {
    const markup = renderToStaticMarkup(<ProjectLibrary mode="workspace" />)
    expect(markup).toContain('aria-label="Target 3D viewer"')
    expect(markup).toContain('title="Open selected structures in 3D viewer A"')
    expect(markup).toContain('3D viewer B (open the B tab first)')
  })

  it('renders one collapsed, keyboard-accessible group rather than every run file', () => {
    workspaceState.active!.artifacts = [
      { id: 'input', file: 'input.pdb', name: 'Source input', kind: 'structure' },
      { id: 'result', file: 'runs/job_123/result/model.pdb', name: 'Generated model', kind: 'structure', command: 'model', run: { command: 'model', status: 'validation_failed' } },
      { id: 'archive', file: 'runs/job_123/diagnostics.tar.gz', name: 'Diagnostic archive', kind: 'artifact', command: 'model' },
    ]
    const markup = renderToStaticMarkup(<ProjectLibrary mode="workspace" />)
    expect(markup).toContain('aria-label="model run job_123"')
    expect(markup).toContain('aria-label="Expand model run"')
    expect(markup).toContain('aria-expanded="false"')
    expect(markup).toContain('aria-pressed="false"')
    expect(markup).toContain('validation failed')
    expect(markup).toContain('Source input')
    expect(markup).not.toContain('Generated model')
    expect(markup).not.toContain('Diagnostic archive')
  })
})
