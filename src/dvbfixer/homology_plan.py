"""Materialize GUI/CLI template-part plans for comparative modeling."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from dvbfixer.homology import AA3TO1
from dvbfixer.model.diffusion.contract import ResidueIdentity, TemplateOwnership
from dvbfixer.salign import run_biopython_superposition


def _residues(pdb: Path, chain: str) -> list[dict]:
    records: list[dict] = []
    current = None
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[21].strip() != chain:
            continue
        aa = AA3TO1.get(line[17:20].strip().upper())
        if not aa:
            continue
        key = line[22:27]
        if current is None or current["key"] != key:
            current = {"key": key, "aa": aa, "lines": []}
            records.append(current)
        current["lines"].append(line)
    if not records:
        raise ValueError(f"no protein residues found for chain {chain} in {pdb}")
    return records


def materialize_template_plan(plan_path: Path, workdir: Path,
                              msa_engine: str = "auto",
                              verbose: bool = False) -> tuple[list[str], str]:
    """Fit selected templates and create one coordinate-preserving known.

    The JSON plan contains ``templates`` and ``alignmentGroups`` using the
    documented Homology workspace schema. Alignment ranges are zero-based,
    half-open column intervals. Earlier template rows win overlap columns.
    """
    plan = json.loads(plan_path.read_text())
    workdir.mkdir(parents=True, exist_ok=True)
    templates = [dict(item) for item in plan.get("templates", [])]
    groups = plan.get("alignmentGroups", [])
    if not templates or not groups:
        raise ValueError("template plan requires templates and alignmentGroups")
    by_id = {item["id"]: item for item in templates}
    global_reference = Path(templates[0]["path"]).resolve()
    fitted_by_id: dict[str, Path] = {}

    for group in groups:
        members = [item for item in templates if item["targetChain"] == group["chainId"]]
        reference = next((item for item in members
                          if Path(item["path"]).resolve() == global_reference), None)
        if len(groups) > 1 and reference is None:
            raise ValueError(
                f"target chain {group['chainId']} has no chain from common "
                f"reference structure {global_reference}"
            )
        if reference:
            members = [reference] + [item for item in members if item["id"] != reference["id"]]
        if len(members) == 1:
            fitted_by_id[members[0]["id"]] = Path(members[0]["path"]).resolve()
            continue
        fit_dir = workdir / "fitted" / str(group["chainId"])
        specs = [f"{Path(item['path']).resolve()}:{item['chain']}" for item in members]
        fitted = run_biopython_superposition(
            specs, workdir / f"structural_alignment_{group['chainId']}.pir",
            fit_dir, msa_engine=msa_engine, verbose=verbose,
        )
        for item, fitted_path in zip(members, fitted):
            fitted_by_id[item["id"]] = fitted_path

    blocks: list[str] = []
    coverage_groups: list[dict[str, Any]] = []
    pdb_lines: list[str] = []
    serial = 1
    first = last = None
    used_chains: set[str] = set()
    group_chains: dict[str, str] = {}
    available = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    for group in groups:
        label = str(group["chainId"])
        preferred = label.upper()[-1] if label.upper() in {"VH", "VL"} else (label if len(label) == 1 else label[0])
        if len(preferred) != 1 or preferred in used_chains:
            preferred = next(candidate for candidate in available if candidate not in used_chains)
        used_chains.add(preferred)
        group_chains[label] = preferred

    for group in groups:
        rows = group.get("rows", [])
        target = next((row for row in rows if row.get("kind") == "target"), None)
        if target is None:
            raise ValueError(f"target row missing for chain {group['chainId']}")
        length = len(target["sequence"])
        if any(len(row["sequence"]) != length for row in rows):
            raise ValueError(f"alignment rows for target chain {group['chainId']} differ in length")
        chosen: list[tuple[dict[str, Any], str] | None] = [None] * length
        for row in (row for row in rows if row.get("kind") == "template"):
            template_id = row.get("templateId", row["id"])
            template = by_id.get(template_id)
            if template is None:
                raise ValueError(f"template metadata missing for {row['id']}")
            residues = _residues(fitted_by_id.get(template_id, Path(template["path"])), template["chain"])
            if row["sequence"].replace("-", "") != "".join(item["aa"] for item in residues):
                raise ValueError(f"alignment row {row['id']} does not match {template['path']}:{template['chain']}")
            by_column = []
            index = 0
            for character in row["sequence"]:
                by_column.append(None if character == "-" else residues[index])
                if character != "-":
                    index += 1
            mode = group.get("maskModes", {}).get(row["id"],
                    "ranges" if group.get("masks", {}).get(row["id"]) else "all")
            spans = [] if mode == "none" else group.get("masks", {}).get(row["id"], []) \
                if mode == "ranges" else [{"start": 0, "end": length}]
            for span in spans:
                for column in range(max(0, span["start"]), min(length, span["end"])):
                    if target["sequence"][column] != "-" and chosen[column] is None:
                        residue = by_column[column]
                        if residue is not None:
                            chosen[column] = (residue, template_id)

        ordinal = 0
        block = []
        coverage: list[dict[str, Any]] = []
        chain = group_chains[str(group["chainId"])]
        for column, target_aa in enumerate(target["sequence"]):
            if target_aa != "-":
                ordinal += 1
            selected = chosen[column]
            residue = selected[0] if selected is not None else None
            owner = selected[1] if selected is not None else None
            block.append(residue["aa"] if residue else "-")
            if target_aa != "-":
                coverage.append({
                    "targetIndex": ordinal - 1,
                    "alignmentColumn": column,
                    "targetResidue": target_aa,
                    "covered": residue is not None,
                    "templateId": owner,
                    "templateResidue": residue["aa"] if residue else None,
                    "pdbResidue": {
                        "chain": chain,
                        "residueNumber": str(ordinal),
                        "insertionCode": "",
                    },
                })
            if residue is None:
                continue
            first = first or (ordinal, chain)
            last = (ordinal, chain)
            for raw in residue["lines"]:
                line = raw.ljust(80)
                pdb_lines.append(
                    f"ATOM  {serial:5d}{line[11:21]}{chain}{ordinal:4d} {line[27:]}\n"
                )
                serial += 1
        blocks.append("".join(block))
        coverage_groups.append({
            "targetChain": str(group["chainId"]),
            "pdbChain": chain,
            "targetSequence": target["sequence"].replace("-", ""),
            "coverage": coverage,
        })
        pdb_lines.append("TER\n")

    if first is None or last is None:
        raise ValueError("template plan selects no target-aligned residues")
    code = "selected_template_mosaic"
    pdb = workdir / f"{code}.pdb"
    pdb.write_text("".join(pdb_lines) + "END\n")
    coverage_path = workdir / "selected_template_mosaic.coverage.json"
    coverage_path.write_text(
        json.dumps(
            {
                "schemaVersion": 1,
                "sourcePlan": str(plan_path.resolve()),
                "mosaic": pdb.name,
                "maskConvention": "zero-based-half-open",
                "overlapPolicy": "earlier-template-wins",
                "groups": coverage_groups,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    target_blocks = [next(row for row in group["rows"] if row.get("kind") == "target")["sequence"]
                     for group in groups]
    pir = workdir / "selected_template_mosaic.pir"
    pir.write_text(
        f">P1;{code}\nstructureX:{code}:{first[0]}:{first[1]}:{last[0]}:{last[1]}::::\n"
        f"{'/'.join(blocks)}*\n>P1;target\nsequence:target::::::::\n{'/'.join(target_blocks)}*\n"
    )
    return [str(pdb)], str(pir)


def prepare_mosaic_diffusion_inputs(
    coverage_path: Path,
    target_chains: list[tuple[str, str]],
    workdir: Path,
) -> tuple[Path, tuple[TemplateOwnership, ...]]:
    """Validate template ownership and write a PDB-chain-keyed target FASTA."""
    raw = json.loads(coverage_path.read_text(encoding="utf-8"))
    if raw.get("schemaVersion") != 1:
        raise ValueError("unsupported mosaic coverage schema")
    targets = dict(target_chains)
    groups = raw.get("groups")
    if not isinstance(groups, list) or not groups:
        raise ValueError("mosaic coverage contains no target groups")

    fasta_lines: list[str] = []
    ownership: list[TemplateOwnership] = []
    seen_pdb_chains: set[str] = set()
    uncovered = 0
    for group in groups:
        logical_chain = str(group.get("targetChain", ""))
        pdb_chain = str(group.get("pdbChain", ""))
        expected_sequence = str(group.get("targetSequence", ""))
        if logical_chain not in targets:
            raise ValueError(
                f"mosaic target chain {logical_chain!r} is absent from target FASTA"
            )
        if targets[logical_chain] != expected_sequence:
            raise ValueError(
                f"mosaic target sequence for chain {logical_chain!r} does not match FASTA"
            )
        if len(pdb_chain) != 1 or pdb_chain in seen_pdb_chains:
            raise ValueError("mosaic PDB chain mapping must be unique and one character")
        seen_pdb_chains.add(pdb_chain)
        coverage = group.get("coverage")
        if not isinstance(coverage, list) or len(coverage) != len(expected_sequence):
            raise ValueError("mosaic coverage length does not match its target sequence")
        for expected_index, item in enumerate(coverage):
            if item.get("targetIndex") != expected_index:
                raise ValueError("mosaic coverage target indices must be contiguous")
            if not item.get("covered"):
                uncovered += 1
                continue
            if item.get("templateResidue") != item.get("targetResidue"):
                raise ValueError(
                    "mosaic diffusion does not yet support template-covered substitutions"
                )
            residue = item.get("pdbResidue", {})
            if residue.get("chain") != pdb_chain:
                raise ValueError("mosaic coverage residue chain mapping is inconsistent")
            ownership.append(TemplateOwnership(
                chain=pdb_chain,
                target_index=expected_index,
                residue=ResidueIdentity(
                    pdb_chain,
                    str(residue.get("residueNumber", "")),
                    str(residue.get("insertionCode", "")),
                ),
                template_id=str(item.get("templateId", "")),
            ))
        fasta_lines.extend((f">{pdb_chain}", expected_sequence))

    extra_targets = set(targets) - {
        str(group.get("targetChain", "")) for group in groups
    }
    if extra_targets:
        raise ValueError(
            "target FASTA contains chains absent from mosaic plan: "
            + ", ".join(sorted(extra_targets))
        )
    if not uncovered:
        raise ValueError("mosaic diffusion requires at least one uncovered target residue")
    fasta_path = workdir / "selected_template_mosaic.target.fasta"
    fasta_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    return fasta_path, tuple(ownership)
