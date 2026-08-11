export interface ConsensusRow { sequence: string }

const STRONG_GROUPS = ['STA', 'NEQK', 'NHQK', 'NDEQ', 'QHRK', 'MILV', 'MILF', 'HY', 'FYW']
const WEAK_GROUPS = ['CSA', 'ATV', 'SAG', 'STNK', 'STPA', 'SGND', 'SNDEQK', 'NDEQHK', 'NEQHRK', 'FVLIM', 'HFY']

/** Clustal-style conservation annotation for a rectangular alignment. */
export function consensusFor(rows: ConsensusRow[]): string {
  const length = rows[0]?.sequence.length || 0
  return Array.from({ length }, (_unused, column) => {
    const values = rows.map(row => row.sequence[column]).filter(value => value && value !== '-')
    if (values.length !== rows.length || !values.length) return ' '
    if (new Set(values).size === 1) return '*'
    if (STRONG_GROUPS.some(group => values.every(value => group.includes(value)))) return ':'
    if (WEAK_GROUPS.some(group => values.every(value => group.includes(value)))) return '.'
    return ' '
  }).join('')
}

/** Per-row comparison against the target/reference alignment. */
export function comparisonToReference(reference: string, sequence: string): string {
  return Array.from({ length: Math.max(reference.length, sequence.length) }, (_unused, column) => {
    const left = reference[column]
    const right = sequence[column]
    if (!left || !right || left === '-' || right === '-') return ' '
    if (left === right) return '*'
    if (STRONG_GROUPS.some(group => group.includes(left) && group.includes(right))) return ':'
    if (WEAK_GROUPS.some(group => group.includes(left) && group.includes(right))) return '.'
    return '×'
  }).join('')
}

export function updateColumnSelection(existing: Iterable<number>, sequence: string,
  column: number, anchor: number, mode: 'replace' | 'toggle' | 'extend'): number[] {
  const selected = new Set(existing)
  if (mode === 'replace') selected.clear()
  if (mode === 'toggle') {
    if (selected.has(column)) selected.delete(column)
    else if (sequence[column] !== '-') selected.add(column)
  } else if (mode === 'replace') {
    if (sequence[column] !== '-') selected.add(column)
  } else {
    for (let index = Math.min(anchor, column); index <= Math.max(anchor, column); index++) {
      if (sequence[index] !== '-') selected.add(index)
    }
  }
  return [...selected].sort((left, right) => left - right)
}
