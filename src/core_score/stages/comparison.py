"""Orchestrate all retained pair comparisons and export their measurements."""

from pathlib import Path
from core_score.stages.utils.comparison_utils import prepare_pair_tasks, run_pair_tasks, write_pair_results


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
