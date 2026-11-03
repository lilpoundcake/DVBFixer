"""Human-readable output from independently validated candidate summaries."""

from pathlib import Path

from dvbfixer.model.diffusion.contract import DiffusionResult


def rejection_output(result: DiffusionResult, bundle: Path) -> tuple[str, ...]:
    """Render only final failures, preserving per-candidate identity and gates."""
    lines: list[str] = []
    for candidate, summary in zip(result.candidates, result.validation_summaries, strict=True):
        if summary.passed:
            continue
        lines.append(
            f"ERROR: Diffusion candidate REJECTED: {candidate.candidate_id} "
            f"(seed {candidate.seed}; final selected coordinates)."
        )
        details = tuple(dict.fromkeys(
            warning.removeprefix("Final validation: ")
            for warning in summary.warnings if warning.startswith("Final validation: ")
        ))
        lines.extend(f"ERROR: {detail}" for detail in details)
        lines.append("ERROR: Failed gates: " + ", ".join(summary.hard_gate_failures))
        if not details:
            lines.append("No residue-level detail is available for these failures; inspect the validation report.")
        lines.append(f"Validation report: {bundle / (candidate.candidate_id + '.diffusion.json')}")
        lines.append("Rejected coordinates are saved for inspection only, not accepted for use.")
    return tuple(lines)
