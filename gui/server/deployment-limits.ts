import fs from 'node:fs'
import path from 'node:path'

export interface DeploymentResourceSettings {
  required: boolean
  dataFilesystemMaxBytes: bigint | null
  temporaryFilesystemMaxBytes: bigint | null
  cpuQuotaPercent: number | null
  memoryMaxBytes: bigint | null
  memorySwapMaxBytes: bigint | null
  tasksMax: bigint | null
  statePaths: readonly string[]
}

export interface DeploymentResourceSnapshot {
  platform: NodeJS.Platform
  cpuMax: string
  memoryMax: string
  memorySwapMax: string
  pidsMax: string
  rootReadOnly: boolean
  unexpectedWritableMounts: readonly string[]
  dataMountPoint: string
  dataCapacityBytes: bigint
  temporaryMounts: readonly { path: string; mountPoint: string; fsType: string; capacityBytes: bigint }[]
  statePathsContained: boolean
}

interface MountRecord {
  mountPoint: string
  fsType: string
  options: readonly string[]
}

function positiveBytes(name: string, value: string | undefined): bigint {
  if (!value || !/^[1-9][0-9]*$/.test(value)) throw new Error(`${name} must be a positive integer`)
  return BigInt(value)
}

function positiveSafeInteger(name: string, value: string | undefined): number {
  const parsed = Number(value)
  if (!value || !Number.isSafeInteger(parsed) || parsed <= 0) throw new Error(`${name} must be a positive safe integer`)
  return parsed
}

export function parseDeploymentResourceSettings(
  environment: NodeJS.ProcessEnv = process.env,
): DeploymentResourceSettings {
  const raw = environment.DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED
  if (raw !== undefined && raw !== '0' && raw !== '1') {
    throw new Error('DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED must be 0 or 1')
  }
  if (raw !== '1') {
    return {
      required: false, dataFilesystemMaxBytes: null, temporaryFilesystemMaxBytes: null,
      cpuQuotaPercent: null, memoryMaxBytes: null, memorySwapMaxBytes: null, tasksMax: null, statePaths: [],
    }
  }
  const statePaths = ['TMPDIR', 'HOME', 'XDG_CACHE_HOME'].map(name => {
    const value = environment[name]?.trim()
    if (!value || !path.isAbsolute(value)) throw new Error(`${name} must be an absolute path when OS limits are required`)
    return path.resolve(value)
  })
  return {
    required: true,
    dataFilesystemMaxBytes: positiveBytes(
      'DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES', environment.DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES,
    ),
    temporaryFilesystemMaxBytes: positiveBytes(
      'DVBFIXER_OS_TEMP_FILESYSTEM_MAX_BYTES', environment.DVBFIXER_OS_TEMP_FILESYSTEM_MAX_BYTES,
    ),
    cpuQuotaPercent: positiveSafeInteger(
      'DVBFIXER_OS_CPU_QUOTA_PERCENT', environment.DVBFIXER_OS_CPU_QUOTA_PERCENT,
    ),
    memoryMaxBytes: positiveBytes(
      'DVBFIXER_OS_MEMORY_MAX_BYTES', environment.DVBFIXER_OS_MEMORY_MAX_BYTES,
    ),
    memorySwapMaxBytes: environment.DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES === '0' ? 0n : positiveBytes(
      'DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES', environment.DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES,
    ),
    tasksMax: positiveBytes('DVBFIXER_OS_TASKS_MAX', environment.DVBFIXER_OS_TASKS_MAX),
    statePaths,
  }
}

function unescapeMountPath(value: string): string {
  return value.replace(/\\(040|011|012|134)/g, (_, octal: string) =>
    String.fromCharCode(Number.parseInt(octal, 8)))
}

function parseMountInfo(content: string): MountRecord[] {
  return content.trim().split('\n').flatMap(line => {
    const [left, right] = line.split(' - ', 2)
    if (!left || !right) return []
    const fields = left.split(' ')
    const rightFields = right.split(' ')
    if (!fields[4] || !rightFields[0]) return []
    return [{
      mountPoint: unescapeMountPath(fields[4]),
      fsType: rightFields[0],
      options: fields[5]?.split(',') || [],
    }]
  })
}

function containingMount(records: readonly MountRecord[], candidate: string): MountRecord {
  const matches = records.filter(record => {
    const relative = path.relative(record.mountPoint, candidate)
    return relative === '' || (!relative.startsWith('..') && !path.isAbsolute(relative))
  }).sort((left, right) => right.mountPoint.length - left.mountPoint.length)
  if (!matches[0]) throw new Error(`no mount found for ${candidate}`)
  return matches[0]
}

function filesystemCapacity(candidate: string): bigint {
  const stats = fs.statfsSync(candidate, { bigint: true })
  return stats.bsize * stats.blocks
}

function canonicalExistingPath(candidate: string): string {
  const missing: string[] = []
  let current = path.resolve(candidate)
  while (!fs.existsSync(current)) {
    const parent = path.dirname(current)
    if (parent === current) break
    missing.unshift(path.basename(current))
    current = parent
  }
  return path.join(fs.realpathSync(current), ...missing)
}

function pathContained(root: string, candidate: string): boolean {
  const relative = path.relative(root, candidate)
  return relative === '' || (!relative.startsWith('..') && !path.isAbsolute(relative))
}

export function inspectDeploymentResources(
  dataRoot: string,
  statePaths: readonly string[],
): DeploymentResourceSnapshot {
  if (process.platform !== 'linux') {
    return {
      platform: process.platform, cpuMax: '', memoryMax: '', pidsMax: '', dataMountPoint: '',
      memorySwapMax: '', rootReadOnly: false, dataCapacityBytes: 0n, temporaryMounts: [],
      unexpectedWritableMounts: [], statePathsContained: false,
    }
  }
  const cgroupLine = fs.readFileSync('/proc/self/cgroup', 'utf8').split('\n')
    .find(line => line.startsWith('0::'))
  if (!cgroupLine) throw new Error('cgroup v2 is required')
  const cgroupDirectory = path.join('/sys/fs/cgroup', cgroupLine.slice(3).replace(/^\/+/, ''))
  const mounts = parseMountInfo(fs.readFileSync('/proc/self/mountinfo', 'utf8'))
  const realDataRoot = fs.realpathSync(dataRoot)
  const dataMount = containingMount(mounts, realDataRoot)
  const rootMountLine = fs.readFileSync('/proc/self/mountinfo', 'utf8').split('\n').find(line => {
    const fields = line.split(' ')
    return fields[4] === '/'
  })
  const temporaryMounts = ['/tmp', '/var/tmp'].map(candidate => {
    const real = fs.realpathSync(candidate)
    const mount = containingMount(mounts, real)
    return { path: candidate, mountPoint: mount.mountPoint, fsType: mount.fsType, capacityBytes: filesystemCapacity(real) }
  })
  const memoryOrKernelFilesystems = new Set([
    'autofs', 'bpf', 'cgroup2', 'configfs', 'debugfs', 'devpts', 'devtmpfs', 'efivarfs',
    'fusectl', 'hugetlbfs', 'mqueue', 'nsfs', 'proc', 'pstore', 'ramfs', 'securityfs',
    'sysfs', 'tmpfs', 'tracefs',
  ])
  const unexpectedWritableMounts = mounts
    .filter(mount => mount.options.includes('rw'))
    .filter(mount => !memoryOrKernelFilesystems.has(mount.fsType))
    .filter(mount => mount.mountPoint !== dataMount.mountPoint)
    .map(mount => mount.mountPoint)
  return {
    platform: process.platform,
    cpuMax: fs.readFileSync(path.join(cgroupDirectory, 'cpu.max'), 'utf8').trim(),
    memoryMax: fs.readFileSync(path.join(cgroupDirectory, 'memory.max'), 'utf8').trim(),
    memorySwapMax: fs.readFileSync(path.join(cgroupDirectory, 'memory.swap.max'), 'utf8').trim(),
    pidsMax: fs.readFileSync(path.join(cgroupDirectory, 'pids.max'), 'utf8').trim(),
    rootReadOnly: rootMountLine?.split(' ')[5]?.split(',').includes('ro') === true,
    unexpectedWritableMounts,
    dataMountPoint: dataMount.mountPoint,
    dataCapacityBytes: filesystemCapacity(realDataRoot),
    temporaryMounts,
    statePathsContained: statePaths.every(candidate =>
      pathContained(realDataRoot, canonicalExistingPath(candidate))),
  }
}

export function validateDeploymentResources(
  snapshot: DeploymentResourceSnapshot,
  dataRoot: string,
  mutationBackupFile: string,
  settings: DeploymentResourceSettings,
): void {
  if (!settings.required) return
  if (snapshot.platform !== 'linux') throw new Error('OS resource enforcement requires Linux')
  const [cpuQuotaText, cpuPeriodText] = snapshot.cpuMax.split(/\s+/)
  const cpuQuota = BigInt(cpuQuotaText === 'max' ? '0' : cpuQuotaText)
  const cpuPeriod = BigInt(cpuPeriodText || '0')
  if (cpuQuotaText === 'max' || cpuQuota <= 0n || cpuPeriod <= 0n ||
      cpuQuota * 100n > cpuPeriod * BigInt(settings.cpuQuotaPercent!)) {
    throw new Error('the service cgroup must set CPUQuota within DVBFIXER_OS_CPU_QUOTA_PERCENT')
  }
  if (snapshot.memoryMax === 'max' || BigInt(snapshot.memoryMax) > settings.memoryMaxBytes!) {
    throw new Error('the service cgroup must set MemoryMax within DVBFIXER_OS_MEMORY_MAX_BYTES')
  }
  if (snapshot.memorySwapMax === 'max' || BigInt(snapshot.memorySwapMax) > settings.memorySwapMaxBytes!) {
    throw new Error('the service cgroup must set MemorySwapMax within DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES')
  }
  if (snapshot.pidsMax === 'max' || BigInt(snapshot.pidsMax) > settings.tasksMax!) {
    throw new Error('the service cgroup must set TasksMax within DVBFIXER_OS_TASKS_MAX')
  }
  if (!snapshot.rootReadOnly) throw new Error('the service root filesystem must be read-only')
  if (snapshot.unexpectedWritableMounts.length) {
    throw new Error(`unexpected writable persistent mounts: ${snapshot.unexpectedWritableMounts.join(', ')}`)
  }
  const realDataRoot = canonicalExistingPath(dataRoot)
  if (snapshot.dataMountPoint !== realDataRoot) {
    throw new Error('DVBFIXER_GUI_DATA_DIR must be the root of a dedicated filesystem')
  }
  if (snapshot.dataCapacityBytes > settings.dataFilesystemMaxBytes!) {
    throw new Error('workspace filesystem exceeds DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES')
  }
  if (!snapshot.statePathsContained || !pathContained(realDataRoot, canonicalExistingPath(mutationBackupFile))) {
    throw new Error('HOME, TMPDIR, XDG_CACHE_HOME, and the mutations backup must be inside workspace storage')
  }
  for (const temporary of snapshot.temporaryMounts) {
    if (temporary.mountPoint !== temporary.path || temporary.fsType !== 'tmpfs' ||
        temporary.capacityBytes > settings.temporaryFilesystemMaxBytes!) {
      throw new Error(`${temporary.path} must be a dedicated bounded tmpfs`)
    }
  }
}

export function assertDeploymentResources(
  dataRoot: string,
  mutationBackupFile: string,
  settings: DeploymentResourceSettings,
): void {
  if (!settings.required) return
  validateDeploymentResources(
    inspectDeploymentResources(dataRoot, settings.statePaths), dataRoot, mutationBackupFile, settings,
  )
}
