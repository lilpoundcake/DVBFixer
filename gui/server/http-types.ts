import type { IncomingMessage, ServerResponse } from 'node:http'

export type ApiNext = (error?: unknown) => void
export type ApiMiddleware = (
  request: IncomingMessage,
  response: ServerResponse,
  next: ApiNext,
) => unknown | Promise<unknown>

export interface ApiMiddlewareRegistrar {
  use(path: string, middleware: ApiMiddleware): unknown
}

export interface ApiRouteHost {
  middlewares: ApiMiddlewareRegistrar
}
