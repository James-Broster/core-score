"""Orchestrate all retained pair comparisons and export their measurements."""

import json
from pathlib import Path
from core_score.results import clean, file_hash, write_csv, write_json
from core_score.stages.utils.comparison_utils import prepare_pair_tasks, run_pair_tasks


def compare_groups(groups, output, workers=1, *, ost="ost", chemlib=None,
                   timeout=None, radius_fallback=True, provenance=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    tasks, rows, details = prepare_pair_tasks(groups, ost=ost, chemlib=chemlib, timeout=timeout)
    for row, detail in run_pair_tasks(
        tasks, workers, ost=ost, chemlib=chemlib,
        timeout=timeout, radius_fallback=radius_fallback,
    ):
        rows.append(row)
        details.append(detail)
    rows.sort(key=lambda row: (row["ensemble"], row["site"], row["pose_i"], row["pose_j"]))
    write_pair_results(groups, rows, details, output, provenance)
    return rows


def write_pair_results(groups, rows, details, output, provenance):
    output = Path(output)
    for group in groups:
        group["input_fingerprints"] = {
            p["pose_id"]: {
                f: file_hash(f)
                for f in (
                    [p["prepared"]["path"]] + p["prepared"].get("ligand_sdfs", [])
                )
            }
            for p in group["poses"]
            if "prepared" in p
        }
    write_csv(output / "pairs.csv", rows)
    with (output / "ligand_details.jsonl").open("w") as handle:
        for detail in details:
            handle.write(json.dumps(clean(detail), allow_nan=False) + "\n")
    write_json(
        output / "analysis_manifest.json",
        dict(schema="core-score.analysis.v1", groups=groups, provenance=provenance),
    )
