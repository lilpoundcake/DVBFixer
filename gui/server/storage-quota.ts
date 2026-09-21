import fs from 'node:fs'
import path from 'node:path'

export const DEFAULT_MAX_UPLOAD_BYTES = 256 * 1024 * 1024
export const DEFAULT_WORKSPACE_QUOTA_BYTES = 5 * 1024 * 1024 * 1024

export interface StorageQuotaSettings {
  maxUploadBytes: number
  workspaceQuotaBytes: number
}

export interface StorageQuotaOptions {
  maxUploadBytes?: number
  workspaceQuotaBytes?: number
}

export interface WorkspaceQuotaExceededDetails {
  quotaBytes: number
  usedBytes: number
  additionalBytes: number
  projectedBytes: number
}

export class WorkspaceQuotaExceededError extends Error {
  readonly statusCode = 507
  readonly code = 'WORKSPACE_QUOTA_EXCEEDED'
  readonly details: WorkspaceQuotaExceededDetails

  constructor(details: WorkspaceQuotaExceededDetails) {
    super(`workspace quota exceeded: ${details.projectedBytes} bytes would exceed ${details.quotaBytes} bytes`)
    this.name = 'WorkspaceQuotaExceededError'
    this.details = details
  }
}

function parsePositiveSafeInteger(name: string, value: string | undefined, fallback: number, allowZero: boolean): number {
  if (value === undefined || value.trim() === '') return fallback
  const parsed = Number(value)
  if (!Number.isSafeInteger(parsed) || parsed < (allowZero ? 0 : 1)) {
    throw new Error(`${name} must be ${allowZero ? 'a non-negative' : 'a positive'} safe integer`)
  }
  return parsed
}

function validateOption(name: string, value: number, allowZero: boolean): number {
  if (!Number.isSafeInteger(value) || value < (allowZero ? 0 : 1)) {
    throw new Error(`${name} must be ${allowZero ? 'a non-negative' : 'a positive'} safe integer`)
  }
  return value
}

export function parseStorageQuotaSettings(
  environment: NodeJS.ProcessEnv = process.env,
): StorageQuotaSettings {
  return {
    maxUploadBytes: parsePositiveSafeInteger(
      'DVBFIXER_GUI_MAX_UPLOAD_BYTES',
      environment.DVBFIXER_GUI_MAX_UPLOAD_BYTES,
      DEFAULT_MAX_UPLOAD_BYTES,
      false,
    ),
    workspaceQuotaBytes: parsePositiveSafeInteger(
      'DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES',
      environment.DVBFIXER_GUI_WORKSPACE_QUOTA_BYTES,
      DEFAULT_WORKSPACE_QUOTA_BYTES,
      true,
    ),
  }
}

let settings = parseStorageQuotaSettings()
let configuredSettings: StorageQuotaSettings | null = null

function resolveStorageQuotaOptions(options: StorageQuotaOptions): StorageQuotaSettings {
  const environment = parseStorageQuotaSettings()
  return Object.freeze({
    maxUploadBytes: options.maxUploadBytes === undefined
      ? environment.maxUploadBytes
      : validateOption('maxUploadBytes', options.maxUploadBytes, false),
    workspaceQuotaBytes: options.workspaceQuotaBytes === undefined
      ? environment.workspaceQuotaBytes
      : validateOption('workspaceQuotaBytes', options.workspaceQuotaBytes, true),
  })
}

export function initializeStorageQuota(options: StorageQuotaOptions = {}): StorageQuotaSettings {
  settings = resolveStorageQuotaOptions(options)
  return settings
}

export function configureStorageQuota(options: StorageQuotaOptions = {}): StorageQuotaSettings {
  const candidate = resolveStorageQuotaOptions(options)
  if (configuredSettings && (
    configuredSettings.maxUploadBytes !== candidate.maxUploadBytes ||
    configuredSettings.workspaceQuotaBytes !== candidate.workspaceQuotaBytes
  )) {
    throw new Error('storage quota is already configured for this process')
  }
  configuredSettings = candidate
  settings = candidate
  return settings
}

export function storageQuotaSettings(): StorageQuotaSettings {
  return settings
}

export function workspaceLogicalBytes(root: string): number {
  if (!fs.existsSync(root)) return 0
  let total = 0
  const visit = (candidate: string): void => {
    const stat = fs.lstatSync(candidate)
    if (stat.isSymbolicLink()) throw new Error(`workspace storage must not contain symlinks: ${candidate}`)
    if (stat.isFile()) {
      total += stat.size
      if (!Number.isSafeInteger(total)) throw new Error('workspace logical byte count exceeds safe integer range')
      return
    }
    if (!stat.isDirectory()) return
    for (const entry of fs.readdirSync(candidate)) visit(path.join(candidate, entry))
  }
  visit(root)
  return total
}

export function availableWorkspaceBytes(root: string): number | null {
  const usedBytes = workspaceLogicalBytes(root)
  if (settings.workspaceQuotaBytes === 0) return null
  return Math.max(0, settings.workspaceQuotaBytes - usedBytes)
}

export function assertWorkspaceQuota(root: string, additionalBytes = 0): void {
  validateOption('additionalBytes', additionalBytes, true)
  const usedBytes = workspaceLogicalBytes(root)
  if (settings.workspaceQuotaBytes === 0) return
  const projectedBytes = usedBytes + additionalBytes
  if (!Number.isSafeInteger(projectedBytes) || projectedBytes > settings.workspaceQuotaBytes) {
    throw new WorkspaceQuotaExceededError({
      quotaBytes: settings.workspaceQuotaBytes,
      usedBytes,
      additionalBytes,
      projectedBytes,
    })
  }
}
