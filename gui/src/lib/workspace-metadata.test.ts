import { describe, expect, it } from 'vitest'
import { structureMetaFromArtifact, workspaceArtifactNameMap } from './workspace-metadata'
import type { WorkspaceArtifact } from '../stores/workspaceStore'

describe('workspace artifact metadata', () => {
  it('creates a complete structure snapshot and clears absent optional fields', () => {
    const artifact: WorkspaceArtifact = { id: 'a', file: 'files/a.pdb', name: 'A', kind: 'structure', resolution: 2.15 }
    expect(structureMetaFromArtifact(artifact)).toEqual({
      name: 'A', organism: '', method: '', resolution: '2.15', description: '',
      iggSubtype: '', allotype: '', equivalentChains: undefined,
    })
  })

  it('maps workspace-relative files to editable artifact names', () => {
    const artifacts: WorkspaceArtifact[] = [
      { id: 'a', file: 'files/a.pdb', name: 'Alpha', kind: 'structure' },
      { id: 'b', file: 'runs/b.pdb', name: 'Beta', kind: 'structure' },
    ]
    expect(workspaceArtifactNameMap(artifacts).get('runs/b.pdb')).toBe('Beta')
  })
})
