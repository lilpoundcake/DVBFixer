import { alignSequences } from './alignment'

/** Parse editable target FASTA into chain-id to ungapped sequence mappings. */
export function targetSequences(fasta: string): Record<string, string> {
  const records: Record<string, string> = {}
  let id = ''
  for (const raw of fasta.split(/\r?\n/)) {
    const line = raw.trim()
    if (!line) continue
    if (line.startsWith('>')) {
      id = line.slice(1).trim().split(/\s+/)[0] || ''
      if (id && records[id] === undefined) records[id] = ''
    } else if (id) {
      records[id] += line.replace(/[\s\-*?]/g, '').toUpperCase()
    }
  }
  return records
}

/** Global pairwise identity percentage, matching the Alignment panel metric. */
export function pairwiseIdentityPercent(target: string, template: string): number | null {
  if (!target || !template) return null
  const alignment = alignSequences(target, template)
  return alignment.length ? alignment.identity / alignment.length * 100 : null
}

export function chainIdentityLabel(
  id: string,
  length: number,
  target: string,
  template: string,
): string {
  const identity = pairwiseIdentityPercent(target, template)
  return `${id} · ${length} aa · ${identity === null ? '— identity' : `${identity.toFixed(1)}% identity`}`
}

export function bestMatchingChainId(
  chains: Array<{ id: string; sequence: string }>,
  target: string,
  current = '',
): string {
  if (!chains.length) return ''
  if (!target) return chains.some(chain => chain.id === current) ? current : chains[0].id
  let best = chains[0]
  let bestIdentity = pairwiseIdentityPercent(target, best.sequence) ?? -1
  for (const chain of chains.slice(1)) {
    const identity = pairwiseIdentityPercent(target, chain.sequence) ?? -1
    if (identity > bestIdentity) {
      best = chain
      bestIdentity = identity
    }
  }
  return best.id
}
