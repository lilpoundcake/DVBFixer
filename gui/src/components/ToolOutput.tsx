import Box from '@mui/material/Box'
import Button from '@mui/material/Button'
import Typography from '@mui/material/Typography'
import Divider from '@mui/material/Divider'
import Alert from '@mui/material/Alert'
import OpenInNewIcon from '@mui/icons-material/OpenInNew'
import type { ToolReportPreview } from '../lib/tool-output'
import { anchorRefinementNotices, diagnosticLine, outputErrors } from '../lib/tool-output'
import { inlineActionSx } from '../ui/layout'

type OutputStream = 'stdout' | 'stderr'

const scrollBlockSx = {
  maxHeight: 'min(30vh, 280px)', minHeight: 64, overflow: 'auto', bgcolor: 'background.default',
  border: 1, borderColor: 'divider', borderRadius: 1,
  '&:focus-visible': { outline: '2px solid', outlineColor: 'primary.main' },
}

function DiagnosticText({ text }: { text: string }) {
  const lines = text.split('\n')
  return <>{lines.map((raw, index) => {
    const line = diagnosticLine(raw)
    return <Box component="span" key={index} data-severity={line.severity} sx={line.severity === 'error'
      ? { color: 'error.main', fontWeight: 700 } : undefined}>
      {line.text}{index < lines.length - 1 ? '\n' : ''}
    </Box>
  })}</>
}

export function ToolOutput({ stdout, stderr, reports = [], outputFile, onOpenFile, onOpenLog }: {
  stdout: string
  stderr: string
  onOpenLog?: (stream: OutputStream) => void
  reports?: ToolReportPreview[]
  outputFile?: string
  onOpenFile?: (file: string) => void
}) {
  const reportTitle = reports.some(report => /diagnose/i.test(report.file)) ? 'Diagnostic report' : 'Output report'
  const anchorNotices = anchorRefinementNotices(stdout, stderr)
  const errors = outputErrors(stdout, stderr, ...reports.map(report => report.text))
  const streamSection = (stream: OutputStream, text: string) => text ? <Box key={stream} sx={{ p: 1 }}>
    <Box sx={{ display: 'flex', flexWrap: 'wrap', alignItems: 'center', gap: 1, mb: 0.5 }}>
      <Typography variant="caption" color="text.secondary" sx={{ fontWeight: 600 }}>
        {stream === 'stdout' ? 'Standard output' : 'Diagnostics & messages (stderr)'}
      </Typography>
      {onOpenLog && <Button size="small" endIcon={<OpenInNewIcon sx={{ fontSize: 12 }} />}
        onClick={() => onOpenLog(stream)} sx={{ ...inlineActionSx, ml: 'auto' }}>
        Open in Text Files
      </Button>}
    </Box>
    <Box component="pre" sx={{
      m: 0, color: 'text.primary', fontFamily: 'ui-monospace, monospace',
      fontSize: '0.72rem', lineHeight: 1.5, whiteSpace: 'pre-wrap', overflowWrap: 'anywhere',
    }}><DiagnosticText text={text} /></Box>
  </Box> : null

  return <Box sx={{ mt: 1 }}>
    {errors.length > 0 && <Alert severity="error" sx={{ mb: 1 }}>
      <Typography variant="body2" sx={{ fontWeight: 700 }}>Errors / validation failures</Typography>
      {errors.slice(0, 12).map(error => <Typography variant="body2" key={error} sx={{ overflowWrap: 'anywhere' }}>{error}</Typography>)}
      {errors.length > 12 && <Typography variant="caption">Additional errors are highlighted in the report and run log below.</Typography>}
    </Alert>}
    {anchorNotices.length > 0 && <Alert severity="warning" sx={{ mb: 1 }}>
      <Typography variant="body2" sx={{ fontWeight: 700 }}>Observed-coordinate diagnostics</Typography>
      {anchorNotices.map(notice => <Typography variant="body2" key={notice}>{notice}</Typography>)}
    </Alert>}
    <Typography variant="caption" sx={{ fontWeight: 700 }}>Output &amp; diagnostics</Typography>
    <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 0.5 }}>
      Reports and progress messages are not necessarily errors. Check the run status above.
    </Typography>
    <Typography variant="caption" sx={{ display: 'block', fontWeight: 700, mb: 0.5 }}>{reportTitle}</Typography>
    <Box role="region" aria-label={reportTitle} tabIndex={0} sx={scrollBlockSx}>
      {reports.map(report => <Box key={report.file} sx={{ p: 1 }}>
        <Box sx={{ display: 'flex', flexWrap: 'wrap', alignItems: 'center', gap: 1, mb: 0.5 }}>
          <Typography variant="caption" sx={{ fontWeight: 600, overflowWrap: 'anywhere' }}>
            Report · {report.file.split('/').pop()}
          </Typography>
          {onOpenFile && !report.inline && <Button size="small" onClick={() => onOpenFile(report.file)}
            endIcon={<OpenInNewIcon sx={{ fontSize: 12 }} />}
            sx={{ ...inlineActionSx, ml: 'auto' }}>
            Open in Text Files
          </Button>}
        </Box>
        {report.error ? <Typography variant="caption" color="text.secondary">
          Preview unavailable: {report.error}. Open the report in Text Files to retry.
        </Typography> : <>
          <Box component="pre" sx={{
            m: 0, color: 'text.primary', fontFamily: 'ui-monospace, monospace',
            fontSize: '0.72rem', lineHeight: 1.5, whiteSpace: 'pre-wrap', overflowWrap: 'anywhere',
          }}><DiagnosticText text={report.text || 'The report file is empty.'} /></Box>
          {report.truncated && <Typography variant="caption" color="text.secondary">
            Preview shortened. Open in Text Files to read the complete report.
          </Typography>}
        </>}
      </Box>)}
      {outputFile && !reports.some(report => report.file === outputFile) && <Box sx={{ p: 1, display: 'flex', flexWrap: 'wrap', alignItems: 'center', gap: 1 }}>
        <Typography variant="caption" sx={{ overflowWrap: 'anywhere' }}>Output file · {outputFile.split('/').pop()}</Typography>
        {onOpenFile && <Button size="small" onClick={() => onOpenFile(outputFile)} sx={{ ...inlineActionSx, ml: 'auto' }}>Open in Text Files</Button>}
      </Box>}
      {!reports.length && !outputFile && <Typography variant="caption" color="text.secondary" sx={{ display: 'block', p: 1 }}>
        This run produced no report file.
      </Typography>}
    </Box>
    <Typography variant="caption" sx={{ display: 'block', fontWeight: 700, mt: 1.5, mb: 0.5 }}>Run log</Typography>
    <Box role="region" aria-label="Run log" tabIndex={0} sx={scrollBlockSx}>
      {streamSection('stderr', stderr)}
      {stderr && stdout && <Divider />}
      {streamSection('stdout', stdout)}
      {!stdout && !stderr && <Typography variant="caption" color="text.secondary" sx={{ display: 'block', p: 1 }}>
        This run produced no console output.
      </Typography>}
    </Box>
  </Box>
}
