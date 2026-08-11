"""Materialize GUI/CLI template-part plans for comparative modeling."""

from __future__ import annotations

import json
from pathlib import Path

from dvbfixer.homology import AA3TO1
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
        chosen = [None] * length
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
                        chosen[column] = by_column[column]

        ordinal = 0
        block = []
        chain = group_chains[str(group["chainId"])]
        for column, target_aa in enumerate(target["sequence"]):
            if target_aa != "-":
                ordinal += 1
            residue = chosen[column]
            block.append(residue["aa"] if residue else "-")
            if residue is None:
                continue
            first = first or (ordinal, chain)
            last = (ordinal, chain)
            for raw in residue["lines"]:
                line = raw.ljust(80)
                pdb_lines.append(
                    f"{line[:6]}{serial:5d}{line[11:21]}{chain}{ordinal:4d} {line[27:]}\n"
                )
                serial += 1
        blocks.append("".join(block))
        pdb_lines.append("TER\n")

    if first is None or last is None:
        raise ValueError("template plan selects no target-aligned residues")
    code = "selected_template_mosaic"
    pdb = workdir / f"{code}.pdb"
    pdb.write_text("".join(pdb_lines) + "END\n")
    target_blocks = [next(row for row in group["rows"] if row.get("kind") == "target")["sequence"]
                     for group in groups]
    pir = workdir / "selected_template_mosaic.pir"
    pir.write_text(
        f">P1;{code}\nstructureX:{code}:{first[0]}:{first[1]}:{last[0]}:{last[1]}::::\n"
        f"{'/'.join(blocks)}*\n>P1;target\nsequence:target::::::::\n{'/'.join(target_blocks)}*\n"
    )
    return [str(pdb)], str(pir)
