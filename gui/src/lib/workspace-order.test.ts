import { describe, expect, it } from 'vitest'
import { reorderVisibleArtifacts } from './workspace-order'

interface Artifact { id: string; file: string; hidden?: boolean }

describe('workspace artifact ordering', () => {
  it('preserves hidden and filtered-out files while moving visible PDB files', () => {
    const artifacts: Artifact[] = [
      { id: 'pdb-a', file: 'a.pdb' },
      { id: 'notes', file: 'notes.txt' },
      { id: 'internal', file: 'run.json', hidden: true },
      { id: 'pdb-b', file: 'b.pdb' },
    ]

    const reordered = reorderVisibleArtifacts(
      artifacts,
      ['pdb-a', 'pdb-b'],
      new Set(['pdb-b']),
      -1,
    )

    expect(reordered.map(item => item.id)).toEqual(['pdb-b', 'notes', 'internal', 'pdb-a'])
    expect(reordered).not.toContain(undefined)
    expect(reordered.find(item => item.id === 'notes')).toBe(artifacts[1])
    expect(reordered.find(item => item.id === 'internal')).toBe(artifacts[2])
  })

  it('moves adjacent selected rows as one stable block', () => {
    const artifacts: Artifact[] = ['a', 'b', 'c', 'd'].map(id => ({ id, file: `${id}.pdb` }))
    const reordered = reorderVisibleArtifacts(
      artifacts,
      artifacts.map(item => item.id),
      new Set(['b', 'c']),
      1,
    )
    expect(reordered.map(item => item.id)).toEqual(['a', 'd', 'b', 'c'])
  })
})
