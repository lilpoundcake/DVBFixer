import type { Plugin } from 'vite'
import { registerApiRoutes } from './api-routes'

export { runDvbfixer } from './dvbfixer-runner'
export { buildArgs } from './command-args'
export { getPg, sseSend, writeSSEHeaders } from './api-routes'

export function apiPlugin(): Plugin {
  return {
    name: 'tarantino-api',
    configureServer(server) {
      registerApiRoutes(server, { projectRoot: server.config.root })
    },
  }
}
