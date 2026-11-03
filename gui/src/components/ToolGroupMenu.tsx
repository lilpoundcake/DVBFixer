import { useEffect, useRef, useState } from 'react'
import Button from '@mui/material/Button'
import Menu from '@mui/material/Menu'
import type { PopoverActions } from '@mui/material/Popover'
import MenuItem from '@mui/material/MenuItem'
import MenuList from '@mui/material/MenuList'
import Paper from '@mui/material/Paper'
import Popper from '@mui/material/Popper'
import ChevronRightIcon from '@mui/icons-material/ChevronRight'
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown'
import { UI } from '../ui/layout'
import { orderToolTabs } from '../lib/tool-tab-order'

interface Tool { name: string; label: string; category: string }

/** One keyboard/touch route and a pointer-hover shortcut to the same tools. */
export function ToolGroupMenu({ commands, category, activeTool, onSelect }: {
  commands: readonly Tool[]; category: string; activeTool?: string
  onSelect: (category: string, tool: string) => void
}) {
  const [anchor, setAnchor] = useState<HTMLElement | null>(null)
  const [group, setGroup] = useState<{ name: string; anchor: HTMLElement } | null>(null)
  const [keyboard, setKeyboard] = useState(false)
  const trigger = useRef<HTMLButtonElement>(null)
  const menuPosition = useRef<PopoverActions | null>(null)
  const updateSubmenu = useRef<() => void>(() => {})
  const [menuWidth, setMenuWidth] = useState<number>()
  useEffect(() => {
    const button = trigger.current
    if (!button || typeof ResizeObserver === 'undefined') return
    const observer = new ResizeObserver(() => {
      setMenuWidth(button.getBoundingClientRect().width)
      menuPosition.current?.updatePosition()
      updateSubmenu.current()
    })
    observer.observe(button)
    return () => observer.disconnect()
  }, [])
  const close = () => { setAnchor(null); setGroup(null); setKeyboard(false) }
  const groups = [...new Set(commands.map(command => command.category))]
  return <>
    <Button ref={trigger} variant="outlined" aria-label={`Workflow group: ${category || 'Select group'}`} aria-haspopup="menu"
      title="Workflow group" endIcon={<ArrowDropDownIcon />}
      sx={{ flex: '0 1 auto', width: 'max-content', minWidth: UI.fieldWidth, maxWidth: '100%', justifyContent: 'space-between' }}
      aria-expanded={!!anchor} aria-controls={anchor ? 'workflow-groups' : undefined}
      onClick={event => setAnchor(event.currentTarget)}>
      <span style={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={category}>
        {category || 'Select group'}
      </span>
    </Button>
    <Menu id="workflow-groups" action={menuPosition} anchorEl={anchor} open={!!anchor} onClose={close}
      disableEnforceFocus slotProps={{ paper: { sx: {
        width: 'max-content', minWidth: menuWidth, maxWidth: 'calc(100vw - 32px)', overflowX: 'auto',
      } } }}>
      {groups.map(name => <MenuItem key={name} selected={category === name}
        aria-haspopup="menu" aria-expanded={group?.name === name}
        onMouseEnter={event => { setKeyboard(false); setGroup({ name, anchor: event.currentTarget }) }}
        onClick={event => { setKeyboard(true); setGroup({ name, anchor: event.currentTarget }) }}
        onKeyDown={event => {
          if (event.key === 'ArrowRight') {
            event.preventDefault(); setKeyboard(true); setGroup({ name, anchor: event.currentTarget })
          }
        }}>
        <span style={{ whiteSpace: 'nowrap' }}>{name}</span><ChevronRightIcon sx={{ ml: 'auto', flexShrink: 0 }} />
      </MenuItem>)}
      <Popper open={!!group} anchorEl={group?.anchor} placement="right-start"
        modifiers={[
          { name: 'flip', options: { fallbackPlacements: ['bottom-start', 'top-start'], padding: 8 } },
          { name: 'preventOverflow', options: { padding: 8 } },
        ]}
        popperRef={instance => { updateSubmenu.current = () => instance?.forceUpdate() }}
        sx={{ zIndex: theme => theme.zIndex.modal + 1 }}>
        <Paper elevation={8} sx={{ maxHeight: '70vh', overflow: 'auto', width: 'max-content', maxWidth: 'calc(100vw - 32px)' }}>
          <MenuList key={`${group?.name}-${keyboard}`} aria-label={`${group?.name} tools`} autoFocusItem={keyboard}
            onKeyDown={event => {
              event.stopPropagation()
              if (event.key === 'ArrowLeft') {
                event.preventDefault(); event.stopPropagation(); group?.anchor.focus(); setGroup(null)
              }
              if (event.key === 'Escape') { event.stopPropagation(); close(); trigger.current?.focus() }
            }}>
            {orderToolTabs(commands.filter(command => command.category === group?.name)).map(command => <MenuItem
              key={command.name} selected={command.name === activeTool}
              onClick={() => { onSelect(command.category, command.name); close(); trigger.current?.focus() }}>
              <span style={{ whiteSpace: 'normal', overflowWrap: 'anywhere' }}>{command.label}</span>
            </MenuItem>)}
          </MenuList>
        </Paper>
      </Popper>
    </Menu>
  </>
}
