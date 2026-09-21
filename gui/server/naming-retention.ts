import fs from 'node:fs'
import path from 'node:path'

export const DEFAULT_NAMING_FAILED_MAX_RUNS = 20
export const DEFAULT_NAMING_FAILED_MAX_AGE_HOURS = 7 * 24

export interface NamingFailureRetention {
  maxRuns: number
  maxAgeHours: number
}

// Only directories created by the naming API are owned by this policy.
const NAMING_RUN = /^naming_[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/

function realDirectory(directory: string): boolean {
  let stat: fs.Stats
  try { stat = fs.lstatSync(directory) } catch (error) {
    if ((error as NodeJS.ErrnoException).code === 'ENOENT') return false
    throw error
  }
  if (!stat.isDirectory() || stat.isSymbolicLink()) {
    throw new Error('naming failure storage must be a real directory')
  }
  return true
}

/** Check both ancestors before accessing or creating the quarantine directory. */
export function namingFailureDirectory(workspaceDirectory: string): string {
  const runs = path.join(workspaceDirectory, 'runs')
  if (!realDirectory(runs)) throw new Error('workspace runs directory does not exist')
  const failed = path.join(runs, '_failed')
  realDirectory(failed)
  return failed
}

/** Lazy pruning under the workspace run lock; successful and active runs are outside this tree. */
export function pruneNamingFailures(
  workspaceDirectory: string,
  settings: NamingFailureRetention,
  now = Date.now(),
): void {
  const failed = namingFailureDirectory(workspaceDirectory)
  if (!fs.existsSync(failed)) return
  const candidates = fs.readdirSync(failed, { withFileTypes: true })
    .filter(entry => entry.isDirectory() && NAMING_RUN.test(entry.name))
    .map(entry => {
      const directory = path.join(failed, entry.name)
      const stat = fs.lstatSync(directory)
      return { directory, stat }
    })
    .filter(({ stat }) => stat.isDirectory() && !stat.isSymbolicLink())
    .sort((a, b) => b.stat.mtimeMs - a.stat.mtimeMs || a.directory.localeCompare(b.directory))
  const cutoff = now - settings.maxAgeHours * 60 * 60 * 1000
  for (const [index, { directory, stat }] of candidates.entries()) {
    if (index >= settings.maxRuns || stat.mtimeMs <= cutoff) {
      fs.rmSync(directory, { recursive: true, force: true })
    }
  }
}
