import type { IncomingMessage, ServerResponse } from 'node:http'
import { ensureRequestId } from './request-observability'

export function isV1Request(request: IncomingMessage): boolean {
  return (request.url || '').startsWith('/v1/')
}

export function v1ErrorBody(
  request: IncomingMessage,
  response: ServerResponse,
  code: string,
  message: string,
): { error: string } | { error: { code: string; message: string; requestId: string } } {
  if (!isV1Request(request)) return { error: message }
  return { error: { code, message, requestId: ensureRequestId(request, response) } }
}
