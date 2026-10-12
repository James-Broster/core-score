"""Orchestrate optional ground-truth comparisons and export their results."""

from pathlib import Path
from core_score.stages.utils.comparison_utils import run_pair_tasks
from core_score.stages.utils.evaluation_utils import prepare_reference_tasks, label_reference_results
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


def summarize_reference(scores, labels, output):
    """Descriptive evaluation only; never fit a CORE cutoff on these same observations."""
    summaries = []
    for score in scores:
        if score["status"] != "complete":
            continue
        for threshold in sorted({r["threshold_angstrom"] for r in labels}):
            matching = {
                r["pose_id"]: r["success"]
                for r in labels
                if r["ensemble"] == score["ensemble"]
                and r["threshold_angstrom"] == threshold
            }
            members = score["group_pose_ids"]
            known = [matching.get(p) for p in members]
            complete = all(v is not None for v in known)
            summaries.append(
                dict(
                    ensemble=score["ensemble"],
                    site=score["site"],
                    core_score=score["S"],
                    threshold_angstrom=threshold,
                    inequality="rmsd < threshold",
                    winning_group_size=len(members),
                    evaluated_members=sum(v is not None for v in known),
                    winning_group_success_fraction=(
                        sum(known) / len(known) if complete and known else None
                    ),
                    note="Reference assignment is OST-native; multicopy labels describe whole poses, not site-specific correctness",
                )
            )
    write_csv(Path(output) / "reference_summary.csv", summaries)
    return summaries
