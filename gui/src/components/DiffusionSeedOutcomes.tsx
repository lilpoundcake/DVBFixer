import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import Typography from '@mui/material/Typography'
import type { DiffusionOutcomes } from '../lib/managed-jobs'
import { inlineActionSx } from '../ui/layout'

/** Presentation only: Python owns every validation outcome and measurement. */
export function DiffusionSeedOutcomes({ diffusion, onOpenReport }: {
  diffusion?: DiffusionOutcomes
  onOpenReport?: (file: string) => void
}) {
  if (!diffusion) return null
  return <Box component="section" aria-label="Diffusion seed outcomes" sx={{ my: 1 }}>
    <Typography variant="subtitle2">Diffusion candidates · {diffusion.outcomes.length} seeds</Typography>
    <Typography variant="caption" display="block" color="text.secondary">
      Per-seed scientific outcomes, separate from process completion. Rejected coordinates are for inspection only.
    </Typography>
    {diffusion.outcomes.map(outcome => <Box key={outcome.seed} sx={{ my: 0.5, overflowWrap: 'anywhere' }}>
      <Typography variant="body2" color={outcome.status === 'accepted' ? 'success.main' : outcome.status === 'rejected' ? 'error' : 'text.secondary'}>
        Seed {outcome.seed}: {outcome.status === 'not_returned' ? 'Not returned for final validation' : outcome.status === 'accepted' ? 'Accepted' : 'Rejected'}
        {outcome.candidateId ? ` · ${outcome.candidateId}` : ''}
      </Typography>
      {outcome.findings.map((finding, index) => <Typography key={index} variant="caption" display="block">{finding}</Typography>)}
      {outcome.hardGateFailures.length > 0 && <Typography variant="caption" display="block">
        Failed gates: {outcome.hardGateFailures.join(', ')}
      </Typography>}
      {outcome.omittedFindings > 0 && <Typography variant="caption" display="block">More findings are available in the retained report.</Typography>}
    </Box>)}
    {onOpenReport && <Button sx={inlineActionSx} aria-label="Open per-seed validation report"
      onClick={() => onOpenReport(diffusion.report)}>Per-seed report</Button>}
  </Box>
}
