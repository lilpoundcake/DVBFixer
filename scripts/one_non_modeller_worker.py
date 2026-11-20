"""Ordinary one-proposal MODELLER with the existing explicit multichain PIR."""
from __future__ import annotations

import sys
from pathlib import Path

from scripts import run_mac_modeller_proposal as proposal
from scripts.run_retained_loop_modeller_controls import multichain_alignment

if __name__ == "__main__":
    request, source, output = map(lambda p: Path(p).resolve(), sys.argv[1:4])
    # Same native initial/fast-loop protocol; only identify the complete template
    # chain range rather than MODELLER's first-chain wildcard interpretation.
    proposal.modeller_alignment = multichain_alignment
    proposal.run(request, source, output, 101)
