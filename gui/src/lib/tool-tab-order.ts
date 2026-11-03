/** Presentation ordering without mutating saved tool identities. */
export function orderToolTabs<T extends { name: string }>(commands: readonly T[]): T[] {
  const ordered = [...commands]
  const diagnose = ordered.findIndex(command => command.name === 'diagnose')
  if (diagnose > 0) ordered.unshift(...ordered.splice(diagnose, 1))
  const altloc = ordered.findIndex(command => command.name === 'altlocs')
  if (altloc >= 0 && ordered.some(command => command.name === 'model')) {
    const [item] = ordered.splice(altloc, 1)
    ordered.splice(ordered.findIndex(command => command.name === 'model'), 0, item)
  }
  return ordered
}
