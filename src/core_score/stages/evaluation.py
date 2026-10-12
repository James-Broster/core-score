"""Orchestrate optional ground-truth comparisons and export their results."""

from pathlib import Path
from core_score.stages.utils.comparison_utils import run_pair_tasks
from core_score.stages.utils.evaluation_utils import prepare_reference_tasks, label_reference_results, summarize_reference
from core_score.results import write_json, write_csv


def evaluate_groups(groups, output, thresholds=(2.0,), workers=1, *, ost="ost",
                    chemlib=None, timeout=None, radius_fallback=True):
    tasks, records = prepare_reference_tasks(groups, output, ost=ost, chemlib=chemlib, timeout=timeout)
    results = run_pair_tasks(tasks, workers, ost=ost, chemlib=chemlib,
                             timeout=timeout, radius_fallback=radius_fallback)
    for task, (row, detail) in zip(tasks, results):
        row.update(pose_id=task[1]["pose_id"], reference=task[2]["path"], detail=detail)
        records.append(row)
    labels = label_reference_results(records, thresholds)
    write_json(Path(output) / "reference_details.json", records)
    write_csv(Path(output) / "reference_evaluation.csv", labels)
    return labels
