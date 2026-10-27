import { expect, it } from 'vitest'
import { openedAltlocArtifact } from './altloc-input'

const artifacts = [
  { id: 'a', file: 'files/input.cif', kind: 'structure' },
  { id: 'b', file: 'files/other.pdb', kind: 'structure' },
  { id: 'data', file: 'files/result.json', kind: 'artifact' },
]
it('prefills and follows the opened structure until the user selects explicitly', () => {
  expect(openedAltlocArtifact(artifacts, 'files/input.cif', '', 'auto')).toBe('a')
  expect(openedAltlocArtifact(artifacts, 'files/other.pdb', 'a', 'auto')).toBe('b')
  expect(openedAltlocArtifact(artifacts, 'files/other.pdb', 'a', 'manual')).toBeNull()
  expect(openedAltlocArtifact(artifacts, 'files/input.cif', 'a', 'auto')).toBeNull()
})
it('does not import an unrelated workspace file or non-structure and recovers a removed selection', () => {
  expect(openedAltlocArtifact(artifacts, 'other-project.pdb', '', 'auto')).toBeNull()
  expect(openedAltlocArtifact(artifacts, 'files/result.json', '', 'auto')).toBeNull()
  expect(openedAltlocArtifact(artifacts, 'files/input.cif', 'removed', 'manual')).toBe('a')
})
