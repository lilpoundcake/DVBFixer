/** Shared GUI/API normalization; Python remains the scientific authority. */
export function diffusionSeeds(value: unknown): number[] {
  if (value === undefined || value === null || value === '' || (typeof value === 'string' && !value.trim())) return [7]
  const tokens = Array.isArray(value) ? value : typeof value === 'string'
    ? value.trim().split(/[,\s]+/).filter(Boolean) : [value]
  if (!tokens.length) {
    if (Array.isArray(value)) return [7]
    throw new Error('Enter at least one seed, or leave the field empty for the default seed 7.')
  }
  const seeds = tokens.map(token => {
    if ((typeof token !== 'string' && typeof token !== 'number') || !/^\d+$/.test(String(token))) {
      throw new Error('Seeds must be non-negative integers separated by commas or spaces, e.g. 1, 17, 148.')
    }
    const seed = Number(token)
    if (!Number.isSafeInteger(seed)) throw new Error('Each seed must be a safely representable non-negative integer.')
    return seed
  })
  if (new Set(seeds).size !== seeds.length) throw new Error('Seeds must be unique; remove duplicate values.')
  return seeds
}

export function diffusionSeedError(value: unknown): string | undefined {
  try { diffusionSeeds(value); return undefined } catch (error) { return (error as Error).message }
}
