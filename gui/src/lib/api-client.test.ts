import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'
import {
  API_TOKEN_STORAGE_KEY,
  apiFetch,
  apiSse,
  clearApiToken,
  fetchApiBlob,
  getApiToken,
  setApiToken,
  subscribeToApiToken,
} from './api-client'

function memoryStorage(): Storage {
  const values = new Map<string, string>()
  return {
    get length() { return values.size },
    clear: () => values.clear(),
    getItem: key => values.get(key) ?? null,
    key: index => [...values.keys()][index] ?? null,
    removeItem: key => { values.delete(key) },
    setItem: (key, value) => { values.set(key, value) },
  }
}

describe('api client credentials', () => {
  let sessionStorage: Storage
  let localStorage: Storage

  beforeEach(() => {
    sessionStorage = memoryStorage()
    localStorage = memoryStorage()
    vi.stubGlobal('window', {
      location: { href: 'https://dvbfixer.test/workspaces/current' },
      sessionStorage,
      localStorage,
    })
  })

  afterEach(() => {
    clearApiToken()
    vi.restoreAllMocks()
    vi.unstubAllGlobals()
  })

  it('stores the token only in sessionStorage and notifies subscribers', () => {
    const subscriber = vi.fn()
    const unsubscribe = subscribeToApiToken(subscriber)

    setApiToken('  secret-value  ')

    expect(getApiToken()).toBe('secret-value')
    expect(sessionStorage.getItem(API_TOKEN_STORAGE_KEY)).toBe('secret-value')
    expect(localStorage.length).toBe(0)
    expect(subscriber).toHaveBeenCalledOnce()

    unsubscribe()
  })

  it('merges bearer authorization into same-origin API requests and preserves request options', async () => {
    setApiToken('secret-value')
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response(null, { status: 204 }))
    vi.stubGlobal('fetch', fetchMock)

    await apiFetch('/api/workspaces', {
      method: 'POST',
      cache: 'no-store',
      credentials: 'same-origin',
      body: '{"name":"test"}',
      headers: { 'Content-Type': 'application/json', 'X-Trace': 'trace-id' },
    })

    const [url, options] = fetchMock.mock.calls[0]
    expect(url).toBe('/api/workspaces')
    expect(options).toMatchObject({
      method: 'POST', cache: 'no-store', credentials: 'same-origin', body: '{"name":"test"}',
    })
    const headers = new Headers(options?.headers)
    expect(headers.get('Authorization')).toBe('Bearer secret-value')
    expect(headers.get('Content-Type')).toBe('application/json')
    expect(headers.get('X-Trace')).toBe('trace-id')
  })

  it('merges Request headers with init headers', async () => {
    setApiToken('secret-value')
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response(null, { status: 200 }))
    vi.stubGlobal('fetch', fetchMock)
    const request = new Request('https://dvbfixer.test/api/session', {
      headers: { 'X-From-Request': 'yes', 'X-Override': 'old' },
    })

    await apiFetch(request, { headers: { 'X-From-Init': 'yes', 'X-Override': 'new' } })

    const headers = new Headers(fetchMock.mock.calls[0][1]?.headers)
    expect(headers.get('X-From-Request')).toBe('yes')
    expect(headers.get('X-From-Init')).toBe('yes')
    expect(headers.get('X-Override')).toBe('new')
    expect(headers.get('Authorization')).toBe('Bearer secret-value')
  })

  it('never sends the token to non-API or cross-origin URLs', async () => {
    setApiToken('secret-value')
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response(null, { status: 200 }))
    vi.stubGlobal('fetch', fetchMock)

    await apiFetch('/structures/model.pdb')
    await apiFetch('https://files.example/api/export?token=public-query')

    for (const [url, options] of fetchMock.mock.calls) {
      expect(new Headers(options?.headers).has('Authorization'), String(url)).toBe(false)
      expect(String(url)).not.toContain('secret-value')
    }
  })

  it('clears the matching token and notifies subscribers after a 401', async () => {
    setApiToken('expired')
    const subscriber = vi.fn()
    const unsubscribe = subscribeToApiToken(subscriber)
    vi.stubGlobal('fetch', vi.fn(async () => new Response(null, { status: 401 })))

    await apiFetch('/api/session')

    expect(getApiToken()).toBeNull()
    expect(subscriber).toHaveBeenCalledOnce()
    unsubscribe()
  })

  it('does not let a stale 401 clear a newer token', async () => {
    setApiToken('old-token')
    let resolveFetch!: (response: Response) => void
    vi.stubGlobal('fetch', vi.fn(() => new Promise<Response>(resolve => { resolveFetch = resolve })))

    const pending = apiFetch('/api/workspaces')
    setApiToken('new-token')
    resolveFetch(new Response(null, { status: 401 }))
    await pending

    expect(getApiToken()).toBe('new-token')
  })

  it('fetches blobs through the authenticated transport', async () => {
    setApiToken('secret-value')
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response('PDB DATA', {
      status: 200,
      headers: { 'Content-Type': 'chemical/x-pdb' },
    }))
    vi.stubGlobal('fetch', fetchMock)

    const blob = await fetchApiBlob('/api/workspaces/a/files/model.pdb')

    expect(await blob.text()).toBe('PDB DATA')
    expect(new Headers(fetchMock.mock.calls[0][1]?.headers).get('Authorization')).toBe('Bearer secret-value')
  })

  it('parses chunked authenticated SSE responses', async () => {
    setApiToken('secret-value')
    const encoder = new TextEncoder()
    const stream = new ReadableStream<Uint8Array>({
      start(controller) {
        controller.enqueue(encoder.encode('event: progress\nid: 7\ndata: first'))
        controller.enqueue(encoder.encode('\ndata: second\n\ndata: done\n\n'))
        controller.close()
      },
    })
    const fetchMock = vi.fn(async (_input: RequestInfo | URL, _init?: RequestInit) => new Response(stream, {
      headers: { 'Content-Type': 'text/event-stream' },
    }))
    vi.stubGlobal('fetch', fetchMock)

    const events = []
    for await (const event of apiSse('/api/jobs/123/events', { method: 'POST' })) events.push(event)

    expect(events).toEqual([
      { event: 'progress', id: '7', data: 'first\nsecond' },
      { data: 'done' },
    ])
    const headers = new Headers(fetchMock.mock.calls[0][1]?.headers)
    expect(headers.get('Accept')).toBe('text/event-stream')
    expect(headers.get('Authorization')).toBe('Bearer secret-value')
  })
})
