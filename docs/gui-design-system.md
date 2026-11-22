# GUI controls and layout

This is the GUI's control/layout contract. Apply it to new panels and when
changing existing ones. It standardizes presentation, not scientific behavior.

## Sources of truth

- [`gui/src/theme.ts`](../gui/src/theme.ts) sets application-wide MUI defaults:
  button, input, label, toggle, checkbox, tab, chip and focus styling.
- [`gui/src/ui/layout.ts`](../gui/src/ui/layout.ts) exports the dimensions and
  reusable toolbar/action styles. Import them rather than copying pixel values.
- [`layout.test.tsx`](../gui/src/ui/layout.test.tsx) checks these contracts.
- Mol* controls remain owned by Mol*. Sequence/alignment residue cells retain
  their scientific grid geometry; they are not ordinary toolbar buttons.

## Button appearance

| Element | Size and appearance | Use |
|---|---|---|
| Default button | 32 px high, at least 80 px wide, 12 px horizontal padding, 12 px sentence-case label | Ordinary actions throughout the application |
| Main workflow action | 160 × 32 px (`primaryActionSx`) | Run, Generate MSA, Parse sequence; secondary actions retain content-driven width, not a forced 160 px |
| Icon-only action | 32 × 32 px, 18 px icon, 6 px corner radius | Refresh, add, focus, remove; always a tooltip and an accessible name |
| Dense-row icon | 24 × 24 px, 14 px icon (`compactIconButtonSx`) | Elements tree, expandable data rows, alignment gap arrows, docking tab chrome |
| Inline text action | 24 px high (`inlineActionSx`), otherwise the normal label style | Open in Text Files inside report headers; not a smaller workflow toolbar |
| Toggle segment | 32 px high, at least 32 px wide | A/B viewer target, severity and numbering selectors |
| Status/legend chip | 24 px high | Passive status, counts, compact filter tags; not a substitute for Run/Import buttons |

Use the default theme size even when specifying `size="small"`. Do not add
one-off 20/22/28 px toolbar controls, shrink button labels, or adjust vertical
position with `translateY`, negative margins or arbitrary padding.

### Variants and states

- **Contained:** the main affirmative action in a group: Run, Save, Generate.
  Avoid several equally prominent contained buttons in one toolbar.
- **Outlined:** secondary actions: Import files, Add new, Structural align.
  Use a white control surface against the tinted toolbar, with a subtle hover
  tint. Unselected toggle segments share that surface; selected segments use a
  distinct blue tint. Disabled controls retain a muted surface and text.
  Outlined text/select fields share the white control surface, including inspector
  thresholds, so inputs do not blend into tinted panel headers.
- **Text:** low-emphasis actions: Reset, Cancel a dialog, open a report.
- **Error color:** destructive/cancel-running actions and actual failed-run
  status. Never make an entire stderr report red just because it is stderr.
- **Disabled/loading:** preserve dimensions and label position. Replace an icon
  with a small spinner without changing button width. Explain prerequisites in
  a tooltip; wrap disabled controls in an inline-flex span if necessary.
- **Keyboard focus:** keep the theme's visible outline. Do not disable outlines
  on application actions. A tooltip is not an accessible name: icon buttons need
  `aria-label`, and toggles should expose selection state.

## Placement and alignment

### Panel/workflow toolbars

Use `panelToolbarSx` for panel headers and `workflowToolbarSx` for Homology
workflow rows. Both use a 40 px minimum height, 4 px vertical padding,
4 px gaps, centered controls and wrapping. Panel headers have an 8 px horizontal
inset; workflow rows share the surrounding content inset. Workflow rows have an opaque
white background matching the content canvas, without a tinted band or extra
outline box. Panel headers retain `background.paper` and a bottom divider.
Rows containing labelled fields reserve equal 8 px top/bottom padding and 12 px
between wrapped rows (48 px minimum height), so controls stay vertically centered
and floating labels stay inside the colored band
without touching its top edge or the previous row. Unlabelled toolbars retain
the compact 40 px baseline. Keep at least 4 px above and below unlabelled
controls so borders do not touch panel dividers. Do not compensate with local label transforms.

All four Homology stages use `workflowContentSx` (8 px horizontal / 4 px
vertical inset) and `workflowStageSx` (4 px between blocks). Workflow toolbars
have no additional horizontal gutter. An unlabelled action row therefore has
8 px from the tab boundary to its controls and 8 px from the controls to the
next block, without stacking independent tab/toolbar/stage padding. Labelled
fields retain the shared extra label clearance. Do not shrink controls or
introduce per-stage spacing patches to compensate for nested containers.
Target's source selector uses the shared flexible command-input style: it fills
the space before Parse sequence. The selector's left edge and the action's right
edge align with the FASTA editor below, without a second toolbar inset.

Place context selectors/fields on the left, a related `actionGroupSx` on the
right, and the primary action at the end of its group. Keep contextual viewer
actions together next to their chain selectors. Labels/status are not spacers:
use flex layout and `ml: 'auto'` for the action group instead of manual offsets.

When fields have helper text, use `formActionRowSx`: align their top edges with
the action group's controls, leaving the caption below the field. Do not center
the action against the combined input-plus-helper height.
Command launchers (including Altlocs) use `commandInputSx` for flexible structure
selectors and `commandActionGroupSx` for actions. Selectors shrink toward 160 px
before wrapping, rather than reserving 320 px and prematurely pushing actions
down. When wrapping is necessary, inputs fill their row and actions stay grouped
at the right; only panels narrower than the action group split individual
actions. Preserve the standard 160×32 main action size; secondary actions such as
Analyze, Reset, Import and Export use the ordinary content-driven button width
with the same 32 px height. Do not impose a toolbar-wide button-width selector.
Short Homology fields use `workflowControlSx` (120–160 px flexible width) to avoid
prematurely orphaning Run while retaining readable labels and unchanged text.

### Fields and text

- Single-line outlined inputs/selects share the 32 px control height. Small
  floating labels and padding are defined **once in the theme**.
- `controlSx` provides 160 px selectors; `wideControlSx` provides up to 320 px
  file/name fields. Content-driven flexible fields may grow with their panel.
- Multiline FASTA/notes editors retain their own height. Do not force multiline
  inputs to 32 px. Preserve residue/grid widths used for alignment selection.
- Use body text at 12 px and supporting captions around 11 px. Do not reduce
  toolbar labels/counts to tiny 8–9 px text to force them onto one line.
- Use `minWidth: 0`, truncation with a title, or `overflowWrap: 'anywhere'` for
  long workspace/file names. Narrow panels wrap whole controls, not button text.
- Do not impose toolbar-wide child selectors that override a deliberately wider
  field. Prefer an explicit `sx={controlSx}` on the field itself.
- Repeated floating-label fields use `labelledFieldRowsSx` (12 px between rows),
  not the 4 px content-block gap. Parsed sequence rows use `sequenceRecordRowSx`:
  a 32 px checkbox column, flexible ID field and content-sized sequence length.
  Adjacent labels must not touch the preceding field's outline.
- Mixed forms use `formStackSx` (also used by `workflowStageSx`). The shared
  `UI.floatingLabelOutset` reserves 8 px before direct labelled fields and rows
  of fields, in addition to the 4 px content gap. This covers Target FASTA after
  an open/closed preview and repeated Templates rows without per-field margins.
  Toolbars advertise the 8 px clearance already supplied by their top padding
  through `--ui-label-clearance`, so the form stack does not add it twice.
  The same token drives the theme's floating-label transform.

### Inspection tables

Workspace run groups use the two-level palette in `workspaceRunColors`: run
headers alternate independently of expanded children; file stripes restart in
each run and use that run's paired background shades. A shared left edge and
indent distinguish children from parent run headers. Parent headers have bold
titles, a small gap before the run and tinted top/bottom separators; light child
stripes use `background.paper` rather than another almost-identical tint.
Expansion must not change
neighboring run bands. Selection and hover keep their shared override colors.

Interactions, Clashes and Altlocs use `InspectionTable` / `InspectionRow`
with `inspectionTableSx`: 28 px baseline rows,
4 px vertical / 8 px horizontal cell padding and alternating white / very light
blue-gray backgrounds (`inspectionRowColors`), lighter than the tinted header.
Use subtle 1 px vertical separators between columns, not around the table's
outer edges. Keep cell backgrounds transparent so row stripes and selection
remain visible across the entire row.
Hover and selected-row colors override the stripes. Rows grow for multiline
content; compact checkboxes fit without increasing the baseline height.
Keep scientific severity/type colors and output checkboxes independent of
the tapped-row highlight. `ResidueLabel` shares case-sensitive monospace identity
formatting; `InspectionBadge` supplies 18 px type/severity labels with 10 px text.
Do not override cell typography or padding locally. Altlocs owns the Review
states / 3D display toolbar; CLI input, policy and execution stay in DVBfixer.
Altlocs includes the chain in its Residue label rather than a duplicate column;
Notes uses the remaining width, with full text available on hover.

Sequence uses readable 12 px monospace letters in 12 px cells, separated by an
8 px gap every ten residues. A content-box ResizeObserver selects the number of
complete ten-residue blocks per line for each dock width (individual cells in
very narrow docks). Resizing changes wrapping, not scientific residue identity,
numbering mode or selection. Ordinary toolbar controls keep their standard size.

MSA residue cells retain their 16×22 px scientific grid. Modeling spans use the
primary fill; transient selection uses a continuous light-blue fill, or a darker
blue over an existing modeling span. Do not outline or inset-shadow each selected
letter: it creates a grid of frames across the range. Only keyboard focus uses a
cell outline (`:focus-visible`); mismatch and active-column underlines remain.
Template-row actions stay inside the 270 px sticky label column. Use the short
`Use span` label, retaining `Use selection as modeling span` as its tooltip and
accessible name. Action groups wrap whole buttons if needed and must never paint
over the residue/comparison grid or change its column alignment.

### Reports

The selected tool owns its status and output; another tool's output must never
appear during switching/loading. The diagnostic/output report and run log are
two separately labelled, keyboard-focusable scroll regions. Their Open in Text
Files actions stay at the right of each header and wrap with the header when
space is limited. Use neutral report backgrounds and text, not an error-colored
stderr container.

Explicit ERROR findings and scientific validation rejections appear in a
deduplicated error banner before warning notices, with bold/error-colored lines
in the report/log preview. Other stderr remains neutral. The frontend displays
Python's declared findings; it does not calculate or reinterpret geometry.

## Examples

### Consistent workflow row

```tsx
<Box sx={panelToolbarSx}>
  <TextField label="Input" sx={wideControlSx} />
  <Box sx={actionGroupSx}>
    <Tooltip title="Refresh inputs">
      <IconButton aria-label="Refresh inputs"><RefreshIcon /></IconButton>
    </Tooltip>
    <Button variant="contained" sx={primaryActionSx}>Run</Button>
  </Box>
</Box>
```

### Avoid local geometry patches

```tsx
// Avoid: a different button baseline and a third control height.
<Button sx={{ height: 28, transform: 'translateY(-3px)' }}>Run</Button>

// Use: shared geometry; fix the parent row, not the button's position.
<Button variant="contained" sx={primaryActionSx}>Run</Button>
```

## Verification checklist

1. Compare normal, disabled, loading, selected and keyboard-focused states.
2. Check narrow docked panels, long filenames/workspace names and toolbar wrapping.
3. Confirm labelled inputs, buttons and icon actions share a baseline; helper text
   and multiline editors must not displace neighboring controls.
4. Confirm tool switching and report navigation still work; leave molecular
   selection, masks, chain identity and scientific calculations unchanged.
5. Run `npm run lint`, `npm run typecheck`, `npm test -- --run` and `npm run build`
   from `gui/` for an application-wide styling change. Use narrow tests while
   iterating; do not silently trigger a full scientific test run.

## Workspace file identity and filtering

[`ArtifactFileIdentity`](../gui/src/components/ArtifactFileIdentity.tsx) renders
one compact row cell: category icon (18 px), flexible ellipsized display name,
and a persistent outlined format chip (24 px high, 52–88 px wide). Category uses
recorded artifact metadata when available and recognizable path extensions as a
fallback. Format always comes from the actual path, never the editable name.
Icon shape, accessible category label and badge text must remain useful without
color; do not use error/warning colors simply to distinguish file types.

The tooltip preserves the complete case-sensitive workspace path and recorded
command/source provenance. Unknown provenance stays explicitly unknown. A format
badge is presentation, not proof that every viewer/backend accepts that format.

**Structures only** is a stable-label outlined toggle action with pressed-state
styling and a tooltip explaining whether clicking shows all files or structures.
The filter and the selection count/ordering arrows form two non-wrapping groups.
There is no redundant FILES heading or vertical divider. The toolbar wraps these groups intact without orphaning
arrows or splitting the selection count. Ordinary controls remain 32 px high. The
structure-only filter includes CIF/mmCIF and other recognized structure formats,
not just `.pdb`. Rendering, range selection and arrow reordering share the same
visible list; switching the filter must not reorder files or unhide internal
bookkeeping artifacts. Classification/filtering policy lives in
[`artifact-presentation.ts`](../gui/src/lib/artifact-presentation.ts).

See also the [GUI overview](../gui/README.md) and
[Homology workspace guide](gui-homology.md).
