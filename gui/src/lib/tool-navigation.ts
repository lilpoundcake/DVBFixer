/** Restore historical group labels without changing tool IDs or panel values. */
export function resolveToolCategory(category: string, categories: readonly string[]): string {
  if (categories.includes(category)) return category
  const aliases: Record<string, string> = {
    'Glycoprotein preparation': 'Glycans',
    Analysis: 'Utilities',
  }
  const renamed = aliases[category]
  return renamed && categories.includes(renamed) ? renamed : categories[0] || ''
}
