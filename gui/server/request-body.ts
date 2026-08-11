import type { IncomingMessage } from 'node:http'

export const MAX_JSON_BODY_BYTES = 2 * 1024 * 1024
export const MAX_UPLOAD_BODY_BYTES = 256 * 1024 * 1024

export class RequestBodyTooLargeError extends Error {
  readonly statusCode = 413

  constructor(limit: number) {
    super(`request body exceeds ${limit} bytes`)
    this.name = 'RequestBodyTooLargeError'
  }
}

export function errorStatus(error: unknown, fallback = 500): number {
  const status = (error as { statusCode?: unknown } | null)?.statusCode
  return typeof status === 'number' && status >= 400 && status <= 599 ? status : fallback
}

export function readRequestBody(
  req: IncomingMessage,
  limit = MAX_JSON_BODY_BYTES,
): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    const declaredLength = Number(req.headers['content-length'])
    if (Number.isFinite(declaredLength) && declaredLength > limit) {
      req.resume()
      reject(new RequestBodyTooLargeError(limit))
      return
    }

    const chunks: Buffer[] = []
    let size = 0
    let settled = false
    req.on('data', chunk => {
      if (settled) return
      const buffer = Buffer.from(chunk)
      size += buffer.length
      if (size > limit) {
        settled = true
        chunks.length = 0
        reject(new RequestBodyTooLargeError(limit))
        return
      }
      chunks.push(buffer)
    })
    req.on('end', () => {
      if (!settled) resolve(Buffer.concat(chunks, size))
    })
    req.on('error', error => {
      if (!settled) reject(error)
    })
  })
}
