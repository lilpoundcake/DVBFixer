import os from 'node:os'
import path from 'node:path'
import { describe, expect, it } from 'vitest'
import {
  parseDeploymentResourceSettings, validateDeploymentResources, type DeploymentResourceSnapshot,
} from './deployment-limits'

function settings() {
  return parseDeploymentResourceSettings({
    DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED: '1',
    DVBFIXER_OS_DATA_FILESYSTEM_MAX_BYTES: '1000',
    DVBFIXER_OS_TEMP_FILESYSTEM_MAX_BYTES: '100',
    DVBFIXER_OS_CPU_QUOTA_PERCENT: '200',
    DVBFIXER_OS_MEMORY_MAX_BYTES: '1000000',
    DVBFIXER_OS_MEMORY_SWAP_MAX_BYTES: '0',
    DVBFIXER_OS_TASKS_MAX: '512',
    TMPDIR: '/srv/dvbfixer/tmp',
    HOME: '/srv/dvbfixer/home',
    XDG_CACHE_HOME: '/srv/dvbfixer/cache',
  })
}

function snapshot(overrides: Partial<DeploymentResourceSnapshot> = {}): DeploymentResourceSnapshot {
  return {
    platform: 'linux', cpuMax: '200000 100000', memoryMax: '1000000', memorySwapMax: '0',
    pidsMax: '512', rootReadOnly: true,
    unexpectedWritableMounts: [],
    dataMountPoint: '/srv/dvbfixer', dataCapacityBytes: 1000n, statePathsContained: true,
    temporaryMounts: [
      { path: '/tmp', mountPoint: '/tmp', fsType: 'tmpfs', capacityBytes: 100n },
      { path: '/var/tmp', mountPoint: '/var/tmp', fsType: 'tmpfs', capacityBytes: 100n },
    ],
    ...overrides,
  }
}

describe('deployment resource limits', () => {
  it('is opt-in and validates required configuration', () => {
    expect(parseDeploymentResourceSettings({}).required).toBe(false)
    expect(() => parseDeploymentResourceSettings({ DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED: 'yes' }))
      .toThrow(/must be 0 or 1/)
    expect(() => parseDeploymentResourceSettings({ DVBFIXER_OS_RESOURCE_LIMITS_REQUIRED: '1' }))
      .toThrow(/TMPDIR/)
  })

  it('accepts a bounded cgroup, dedicated data filesystem, and private temporary filesystems', () => {
    expect(() => validateDeploymentResources(
      snapshot(), '/srv/dvbfixer', '/srv/dvbfixer/mutations.json', settings(),
    )).not.toThrow()
  })

  it.each([
    [{ cpuMax: 'max 100000' }, /CPUQuota/],
    [{ memoryMax: 'max' }, /MemoryMax/],
    [{ memoryMax: '1000001' }, /MemoryMax/],
    [{ memorySwapMax: 'max' }, /MemorySwapMax/],
    [{ pidsMax: 'max' }, /TasksMax/],
    [{ pidsMax: '513' }, /TasksMax/],
    [{ rootReadOnly: false }, /root filesystem/],
    [{ unexpectedWritableMounts: ['/mnt/scratch'] }, /unexpected writable persistent mounts/],
    [{ dataMountPoint: '/' }, /dedicated filesystem/],
    [{ dataCapacityBytes: 1001n }, /filesystem exceeds/],
    [{ statePathsContained: false }, /must be inside/],
    [{ temporaryMounts: [{ path: '/tmp', mountPoint: '/', fsType: 'ext4', capacityBytes: 100n }] }, /bounded tmpfs/],
  ] as const)('fails closed for an incomplete resource boundary', (overrides, message) => {
    expect(() => validateDeploymentResources(
      snapshot(overrides), '/srv/dvbfixer', '/srv/dvbfixer/mutations.json', settings(),
    )).toThrow(message)
  })

  it('rejects the deployment profile on non-Linux hosts', () => {
    expect(() => validateDeploymentResources(
      snapshot({ platform: os.platform() === 'linux' ? 'darwin' : os.platform() }),
      '/srv/dvbfixer', '/srv/dvbfixer/mutations.json', settings(),
    )).toThrow(/requires Linux/)
  })

  it('requires the mutable backup to remain in the bounded data filesystem', () => {
    expect(() => validateDeploymentResources(
      snapshot(), '/srv/dvbfixer', path.resolve('/opt/dvbfixer/mutations.json'), settings(),
    )).toThrow(/mutations backup/)
  })
})
