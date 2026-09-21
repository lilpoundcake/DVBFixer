import { useSyncExternalStore } from 'react'
import { getApiToken, subscribeToApiToken } from '../lib/api-client'

export function useApiToken(): string | null {
  return useSyncExternalStore(subscribeToApiToken, getApiToken, () => null)
}
