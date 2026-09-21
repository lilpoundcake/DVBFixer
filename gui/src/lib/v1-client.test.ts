import { afterEach, describe, expect, it, vi } from 'vitest'
import { clearApiToken, setApiToken } from './api-client'
import { createV1Client } from './v1-client'

function memoryStorage(): Storage {
  const values = new Map<string, string>()
  return {
    get length() { return values.size }, clear: () => values.clear(),
    getItem: key => values.get(key) ?? null, key: index => [...values.keys()][index] ?? null,
    removeItem: key => { values.delete(key) }, setItem: (key, value) => { values.set(key, value) },
  }
}

afterEach(() => {
  clearApiToken()
  vi.unstubAllGlobals()
})

describe('generated V1 client', () => {
  it('uses generated paths with the authenticated GUI transport', async () => {
    vi.stubGlobal('window', {
      location: { href: 'https://dvbfixer.test/workspaces/current', origin: 'https://dvbfixer.test' },
      sessionStorage: memoryStorage(),
    })
    setApiToken('generated-client-token')
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response('[]', {
      status: 200, headers: { 'Content-Type': 'application/json' },
    }))
    vi.stubGlobal('fetch', fetchMock)
    const client = createV1Client('https://dvbfixer.test')

    const result = await client.GET('/api/v1/workspaces/{workspaceId}/jobs', {
      params: { path: { workspaceId: 'workspace-a' } },
    })

    expect(result.data).toEqual([])
    const request = fetchMock.mock.calls[0][0] as Request
    const options = fetchMock.mock.calls[0][1] as RequestInit | undefined
    expect(request.url).toContain('/api/v1/workspaces/workspace-a/jobs')
    expect(new Headers(options?.headers).get('Authorization'))
      .toBe('Bearer generated-client-token')
  })
})
