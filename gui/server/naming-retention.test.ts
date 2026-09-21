import { afterEach, describe, expect, it } from 'vitest'
import crypto from 'node:crypto'
import fs from 'node:fs'
import os from 'node:os'
import path from 'node:path'
import { namingFailureDirectory, pruneNamingFailures } from './naming-retention'

const temporaryDirectories: string[] = []
const now = 1_800_000_000_000

function workspace(): string {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'dvbfixer-retention-'))
  temporaryDirectories.push(root)
  fs.mkdirSync(path.join(root, 'runs', '_failed'), { recursive: true })
  return root
}

function failure(root: string, ageHours: number): string {
  const directory = path.join(root, 'runs', '_failed', `naming_${crypto.randomUUID()}`)
  fs.mkdirSync(directory)
  fs.writeFileSync(path.join(directory, 'stderr.log'), 'diagnostics')
  const timestamp = new Date(now - ageHours * 3_600_000)
  fs.utimesSync(directory, timestamp, timestamp)
  return directory
}

afterEach(() => {
  temporaryDirectories.splice(0).forEach(root => fs.rmSync(root, { recursive: true, force: true }))
})

describe('naming failure retention', () => {
  it('keeps the newest failures subject to both age and count', () => {
    const root = workspace()
    const recent = failure(root, 1)
    const middle = failure(root, 2)
    const overCount = failure(root, 3)
    const expired = failure(root, 168)
    pruneNamingFailures(root, { maxRuns: 2, maxAgeHours: 168 }, now)
    expect(fs.existsSync(recent)).toBe(true)
    expect(fs.existsSync(middle)).toBe(true)
    expect(fs.existsSync(overCount)).toBe(false)
    expect(fs.existsSync(expired)).toBe(false)
    expect(fs.readFileSync(path.join(recent, 'stderr.log'), 'utf8')).toBe('diagnostics')
  })

  it('expires old failures even when the count limit is not reached', () => {
    const root = workspace()
    const expired = failure(root, 25)
    pruneNamingFailures(root, { maxRuns: 20, maxAgeHours: 24 }, now)
    expect(fs.existsSync(expired)).toBe(false)
  })

  it('never prunes active, successful, unrelated, or other-workspace runs', () => {
    const root = workspace()
    const otherRoot = workspace()
    const otherFailure = failure(otherRoot, 200)
    const protectedPaths = [
      path.join(root, 'runs', `naming_${crypto.randomUUID()}`),
      path.join(root, 'runs', '_failed', 'model_old'),
      path.join(root, 'runs', '_failed', 'naming_manual'),
    ]
    for (const directory of protectedPaths) {
      fs.mkdirSync(directory)
      fs.writeFileSync(path.join(directory, 'output.pdb'), 'protected')
    }
    failure(root, 200)
    pruneNamingFailures(root, { maxRuns: 1, maxAgeHours: 1 }, now)
    for (const directory of protectedPaths) {
      expect(fs.readFileSync(path.join(directory, 'output.pdb'), 'utf8')).toBe('protected')
    }
    expect(fs.existsSync(otherFailure)).toBe(true)
  })

  it('does not follow symlinks in the failure area or its nested contents', () => {
    const root = workspace()
    const outside = workspace()
    const protectedFile = path.join(outside, 'keep.txt')
    fs.writeFileSync(protectedFile, 'keep')
    const link = path.join(root, 'runs', '_failed', `naming_${crypto.randomUUID()}`)
    fs.symlinkSync(outside, link, 'dir')
    const expired = failure(root, 200)
    fs.symlinkSync(outside, path.join(expired, 'linked'), 'dir')
    fs.utimesSync(expired, new Date(0), new Date(0))
    pruneNamingFailures(root, { maxRuns: 1, maxAgeHours: 1 }, now)
    expect(fs.lstatSync(link).isSymbolicLink()).toBe(true)
    expect(fs.existsSync(expired)).toBe(false)
    expect(fs.readFileSync(protectedFile, 'utf8')).toBe('keep')
  })

  it.each(['runs', '_failed'])('rejects a symlinked %s ancestor', component => {
    const root = workspace()
    const outside = workspace()
    const target = component === 'runs' ? path.join(root, 'runs') : path.join(root, 'runs', '_failed')
    fs.rmSync(target, { recursive: true })
    fs.symlinkSync(outside, target, 'dir')
    expect(() => pruneNamingFailures(root, { maxRuns: 1, maxAgeHours: 1 }, now)).toThrow('real directory')
    expect(() => namingFailureDirectory(root)).toThrow('real directory')
  })

  it('does not create a quarantine directory just to prune it', () => {
    const root = workspace()
    const failed = path.join(root, 'runs', '_failed')
    fs.rmdirSync(failed)
    pruneNamingFailures(root, { maxRuns: 1, maxAgeHours: 1 }, now)
    expect(fs.existsSync(failed)).toBe(false)
  })
})
