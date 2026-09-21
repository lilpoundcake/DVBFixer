export const API_TOKEN_STORAGE_KEY = 'dvbfixer.apiToken'

type TokenSubscriber = () => void

const tokenSubscribers = new Set<TokenSubscriber>()

function getSessionStorage(): Storage | null {
  if (typeof window === 'undefined') return null
  try {
    return window.sessionStorage
  } catch {
    return null
  }
}

function notifyTokenSubscribers(): void {
  for (const subscriber of tokenSubscribers) subscriber()
}

export function getApiToken(): string | null {
  return getSessionStorage()?.getItem(API_TOKEN_STORAGE_KEY) ?? null
}

export function setApiToken(token: string): void {
  const normalized = token.trim()
  if (!normalized) {
    clearApiToken()
    return
  }
  getSessionStorage()?.setItem(API_TOKEN_STORAGE_KEY, normalized)
  notifyTokenSubscribers()
}

export function clearApiToken(): void {
  getSessionStorage()?.removeItem(API_TOKEN_STORAGE_KEY)
  notifyTokenSubscribers()
}

export function subscribeToApiToken(subscriber: TokenSubscriber): () => void {
  tokenSubscribers.add(subscriber)
  return () => tokenSubscribers.delete(subscriber)
}

function currentLocation(): URL {
  return new URL(typeof window === 'undefined' ? 'http://localhost/' : window.location.href)
}

function requestUrl(input: RequestInfo | URL): URL {
  const rawUrl = input instanceof Request ? input.url : String(input)
  return new URL(rawUrl, currentLocation())
}

function isSameOriginApiRequest(input: RequestInfo | URL): boolean {
  const location = currentLocation()
  const url = requestUrl(input)
  return url.origin === location.origin && (url.pathname === '/api' || url.pathname.startsWith('/api/'))
}

function mergedHeaders(input: RequestInfo | URL, init?: RequestInit): Headers {
  const headers = new Headers(input instanceof Request ? input.headers : undefined)
  if (init?.headers) {
    new Headers(init.headers).forEach((value, name) => headers.set(name, value))
  }
  return headers
}

/** Fetch with a session-scoped bearer credential for same-origin API routes. */
export async function apiFetch(input: RequestInfo | URL, init?: RequestInit): Promise<Response> {
  const apiRequest = isSameOriginApiRequest(input)
  const token = apiRequest ? getApiToken() : null
  const headers = mergedHeaders(input, init)
  if (token) headers.set('Authorization', `Bearer ${token}`)

  const response = await fetch(input, { ...init, headers })
  // Do not let an old in-flight request clear a credential entered more recently.
  if (apiRequest && token && response.status === 401 && getApiToken() === token) clearApiToken()
  return response
}

export async function fetchApiBlob(input: RequestInfo | URL, init?: RequestInit): Promise<Blob> {
  const response = await apiFetch(input, init)
  if (!response.ok) throw new Error(`Download failed: HTTP ${response.status}`)
  return response.blob()
}

export async function downloadApiBlob(
  input: RequestInfo | URL,
  filename: string,
  init?: RequestInit,
): Promise<void> {
  const blob = await fetchApiBlob(input, init)
  const blobUrl = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = blobUrl
  link.download = filename
  link.style.display = 'none'
  document.body.appendChild(link)
  link.click()
  link.remove()
  window.setTimeout(() => URL.revokeObjectURL(blobUrl), 0)
}

export async function openApiBlob(
  input: RequestInfo | URL,
  init?: RequestInit,
  target = '_blank',
): Promise<void> {
  // Open synchronously so popup blockers do not reject the window after fetch resolves.
  const openedWindow = window.open('', target)
  try {
    const blob = await fetchApiBlob(input, init)
    const blobUrl = URL.createObjectURL(blob)
    if (openedWindow) openedWindow.location.replace(blobUrl)
    else window.open(blobUrl, target, 'noopener,noreferrer')
    window.setTimeout(() => URL.revokeObjectURL(blobUrl), 60_000)
  } catch (error) {
    openedWindow?.close()
    throw error
  }
}

export interface ApiSseEvent {
  data: string
  event?: string
  id?: string
  retry?: number
}

function parseSseFrame(frame: string): ApiSseEvent | null {
  const event: ApiSseEvent = { data: '' }
  const data: string[] = []
  for (const line of frame.split(/\r?\n/)) {
    if (!line || line.startsWith(':')) continue
    const separator = line.indexOf(':')
    const field = separator < 0 ? line : line.slice(0, separator)
    let value = separator < 0 ? '' : line.slice(separator + 1)
    if (value.startsWith(' ')) value = value.slice(1)
    if (field === 'data') data.push(value)
    else if (field === 'event') event.event = value
    else if (field === 'id' && !value.includes('\0')) event.id = value
    else if (field === 'retry' && /^\d+$/.test(value)) event.retry = Number(value)
  }
  if (data.length === 0) return null
  event.data = data.join('\n')
  return event
}

/** Consume an authenticated SSE response, including POST-based event streams. */
export async function* apiSse(
  input: RequestInfo | URL,
  init?: RequestInit,
): AsyncGenerator<ApiSseEvent, void, undefined> {
  const headers = mergedHeaders(input, init)
  if (!headers.has('Accept')) headers.set('Accept', 'text/event-stream')
  const response = await apiFetch(input, { ...init, headers })
  if (!response.ok) throw new Error(`Event stream failed: HTTP ${response.status}`)
  if (!response.body) throw new Error('Event stream response has no body')

  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  let buffer = ''
  try {
    while (true) {
      const { done, value } = await reader.read()
      buffer += decoder.decode(value, { stream: !done })
      const frames = buffer.split(/\r?\n\r?\n/)
      buffer = frames.pop() ?? ''
      for (const frame of frames) {
        const event = parseSseFrame(frame)
        if (event) yield event
      }
      if (done) break
    }
    const finalEvent = parseSseFrame(buffer)
    if (finalEvent) yield finalEvent
  } finally {
    reader.releaseLock()
  }
}
