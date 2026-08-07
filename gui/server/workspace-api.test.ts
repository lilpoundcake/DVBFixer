import { afterEach, describe, expect, it } from 'vitest'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { ensureWorkspaceMigration, listWorkspaces, loadWorkspace, resolveWorkspaceFile, saveWorkspace } from './workspace-api'

const temporaryDirectories: string[] = []
function temp(): string {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-workspace-test-'))
  temporaryDirectories.push(directory)
  return directory
}
afterEach(() => temporaryDirectories.splice(0).forEach(directory => fs.rmSync(directory, { recursive: true, force: true })))

describe('workspace storage', () => {
  it('migrates top-level folders to projects and ungrouped files to Unsorted', () => {
    const root = temp()
    fs.mkdirSync(path.join(root, 'legacy'))
    fs.writeFileSync(path.join(root, 'legacy', 'template.pdb'), 'END\n')
    fs.writeFileSync(path.join(root, 'target.fasta'), '>A\nAGS\n')
    fs.writeFileSync(path.join(root, 'index.json'), JSON.stringify([
      { id: '__root__', kind: 'folder', name: '__root__', children: ['folder-one', 'target.fasta'] },
      { id: 'folder-one', kind: 'folder', name: 'Project One', children: ['legacy/template.pdb'] },
      { id: 'legacy/template', file: 'legacy/template.pdb', name: 'Template', kind: 'structure' },
      { id: 'target', file: 'target.fasta', name: 'Target', kind: 'artifact' },
    ]))
    ensureWorkspaceMigration(root)
    const projects = listWorkspaces(root)
    expect(projects.map(project => project.name).sort()).toEqual(['Project One', 'Unsorted'])
    const project = projects.find(item => item.name === 'Project One')!
    expect(project.artifacts[0].file).toMatch(/^files\//)
    expect(fs.existsSync(resolveWorkspaceFile(root, project.id, project.artifacts[0].file))).toBe(true)
    expect(() => resolveWorkspaceFile(root, project.id, '../index.json')).toThrow(/workspace artifact/)
  })

  it('persists versioned workflow state', () => {
    const root = temp()
    ensureWorkspaceMigration(root)
    const workspace = listWorkspaces(root)[0]
    const saved = saveWorkspace(root, { ...workspace, toolState: { dvbfixer: { inputFile: 'a.pdb' } } })
    expect(loadWorkspace(root, saved.id).toolState).toEqual({ dvbfixer: { inputFile: 'a.pdb' } })
  })
})
