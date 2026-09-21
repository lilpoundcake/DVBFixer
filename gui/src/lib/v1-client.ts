import createClient from 'openapi-fetch'
import type { paths } from '../generated/dvbfixer-api'
import { apiFetch } from './api-client'

/** Generated OpenAPI types with the GUI's authenticated same-origin transport. */
export function createV1Client(
  baseUrl = typeof window === 'undefined' ? 'http://localhost' : window.location.origin,
) {
  return createClient<paths>({ baseUrl, fetch: apiFetch })
}

export const v1Client = createV1Client()
