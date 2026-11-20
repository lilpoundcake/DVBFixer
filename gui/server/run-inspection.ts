/** Read-only diagnostic preview. Never extracts files or grants scientific acceptance. */
import crypto from 'node:crypto'
import fs from 'node:fs'
import path from 'node:path'
import { createGunzip } from 'node:zlib'

const MAX_PDB = 20 * 1024 * 1024
const MAX_EXPANDED = 512 * 1024 * 1024
interface Entry { file: string; bytes: number; sha256: string }
interface Inventory { status?: string; outputs?: string[]; diagnostics?: { archive: string; files: Entry[] } | null }
const safe = (file: string) => !!file && !file.includes('\\') && !path.posix.isAbsolute(file) &&
  file.split('/').every(part => part !== '.' && part !== '..' && !!part && !part.includes(':'))

export function inspectionInventory(summary: string): Inventory {
  if (fs.statSync(summary).size > 10 * 1024 * 1024) throw new Error('run summary too large')
  const data = JSON.parse(fs.readFileSync(summary, 'utf8'))
  if (data.version !== 1 || !Array.isArray(data.outputs) ||
      !data.outputs.every((file: unknown) => typeof file === 'string' && safe(file)) ||
      (data.diagnostics && (data.diagnostics.archive !== 'diagnostics.tar.gz' || !Array.isArray(data.diagnostics.files) ||
        !data.diagnostics.files.every((entry: Entry) => entry && typeof entry.file === 'string' && safe(entry.file) &&
          Number.isSafeInteger(entry.bytes) && entry.bytes >= 0 && /^[a-f0-9]{64}$/.test(entry.sha256))))) {
    throw new Error('invalid run summary')
  }
  return data
}

export function inspectionStructures(inventory: Inventory): string[] {
  // Accepted outputs are deliberately excluded. Other PDBs may be raw proposals,
  // intermediates or rejected candidates; do not claim they all failed validation.
  return [...new Set([
    ...(inventory.outputs || []).filter(file => /\.pdb$/i.test(file) &&
      (inventory.status === 'validation_failed' || file.split('/').includes('rejected'))),
    ...(inventory.diagnostics?.files || []).filter(entry => /\.pdb$/i.test(entry.file) && entry.bytes > 0 && entry.bytes <= MAX_PDB)
      .map(entry => entry.file),
  ])]
}

export async function readInspectionStructure(archive: string, inventory: Inventory, file: string): Promise<Buffer> {
  const expected = inventory.diagnostics?.files.find(entry => entry.file === file)
  if (!expected || !inspectionStructures(inventory).includes(file)) throw new Error('diagnostic structure not found')
  const source = fs.createReadStream(archive)
  const unzip = source.pipe(createGunzip())
  source.on('error', error => unzip.destroy(error))
  let buffer = Buffer.alloc(0), expanded = 0, remaining = 0, padding = 0, selected = false
  const chunks: Buffer[] = []
  try {
    for await (const chunk of unzip) {
      expanded += chunk.length
      if (expanded > MAX_EXPANDED) throw new Error('diagnostic archive exceeds preview limit')
      buffer = Buffer.concat([buffer, chunk])
      while (buffer.length) {
        if (remaining) {
          const count = Math.min(remaining, buffer.length)
          if (selected) chunks.push(Buffer.from(buffer.subarray(0, count)))
          buffer = buffer.subarray(count)
          remaining -= count
          if (remaining) break
          if (selected) {
            const result = Buffer.concat(chunks)
            if (result.length !== expected.bytes || crypto.createHash('sha256').update(result).digest('hex') !== expected.sha256) {
              throw new Error('diagnostic structure checksum mismatch')
            }
            return result
          }
        }
        if (padding) {
          const count = Math.min(padding, buffer.length)
          buffer = buffer.subarray(count); padding -= count
          if (padding) break
        }
        if (buffer.length < 512) break
        const header = buffer.subarray(0, 512)
        buffer = buffer.subarray(512)
        const field = (start: number, size: number) => header.subarray(start, start + size).toString('utf8').replace(/\0.*$/s, '')
        const name = [field(345, 155), field(0, 100)].filter(Boolean).join('/')
        const size = field(124, 12).trim()
        if (!/^[0-7]+$/.test(size)) throw new Error('invalid diagnostic archive header')
        remaining = parseInt(size, 8)
        if (!Number.isSafeInteger(remaining) || remaining > MAX_EXPANDED) throw new Error('invalid diagnostic archive size')
        padding = (512 - remaining % 512) % 512
        selected = name === file
        if (selected && (remaining !== expected.bytes || remaining > MAX_PDB || !['0', '\0'].includes(header.subarray(156, 157).toString()))) {
          throw new Error('invalid diagnostic structure entry')
        }
      }
    }
    throw new Error('diagnostic structure not found in archive')
  } finally { source.destroy(); unzip.destroy() }
}
