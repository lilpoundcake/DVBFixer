import { Readable } from 'node:stream'
import { describe, expect, it } from 'vitest'
import { readRequestBody, RequestBodyTooLargeError } from './request-body'

function request(chunks: string[], contentLength?: number) {
  const stream = Readable.from(chunks) as Readable & { headers: Record<string, string> }
  stream.headers = contentLength === undefined ? {} : { 'content-length': String(contentLength) }
  return stream
}

describe('bounded request bodies', () => {
  it('returns bodies within the configured limit', async () => {
    await expect(readRequestBody(request(['ab', 'cd']) as any, 4)).resolves.toEqual(Buffer.from('abcd'))
  })

  it('rejects declared and streamed bodies beyond the configured limit', async () => {
    await expect(readRequestBody(request([] , 5) as any, 4)).rejects.toBeInstanceOf(RequestBodyTooLargeError)
    await expect(readRequestBody(request(['abc', 'de']) as any, 4)).rejects.toBeInstanceOf(RequestBodyTooLargeError)
  })
})
