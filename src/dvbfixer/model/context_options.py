"""Public context choices, with an unadvertised legacy replay spelling."""

class _ContextChoices(tuple):
    def __contains__(self, value):
        return value == "none" or super().__contains__(value)


def add_context_argument(group):
    group.add_argument(
        "--diffusion-context", choices=_ContextChoices(("auto", "local", "full")), default="auto",
        help="Protenix input: auto uses the full structure within the supported 512-residue budget, "
             "otherwise a rule-selected motif; local forces the motif, full requires the complete input. "
             "Protpardelle always uses its rule-selected motif (default: auto).",
    )
