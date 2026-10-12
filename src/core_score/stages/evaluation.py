"""Optional comparison of every prediction against its ensemble's ground-truth structure."""

import math
from core_score.results import fingerprint, write_csv, write_json
from core_score.stages.comparison import run_pair_tasks
from core_score.stages.preparation import prepare


def evaluate_groups(groups, scores, output, thresholds=(2.0,), *, workers=1, timeout=None, radius_fallback=True):
    """Ligand RMSD of each prepared prediction to the reference, labelled at each threshold.
    Reference results never enter CORE."""
    tasks, records = [], []
    for group in groups:
        destination = output / "references" / fingerprint(group["ensemble"])[:20] / "reference.cif"
        reference = prepare(group["reference"], destination, timeout=timeout)
        for pose in group["poses"]:
            error = reference.get("error") or pose.get("error")
            if error:
                records.append(dict(ensemble=group["ensemble"], pose_id=pose["pose_id"], rmsd=None, error=error))
            else:
                tasks.append((dict(group, site="reference"), pose, reference))
    for row, detail in run_pair_tasks(tasks, workers, timeout=timeout, radius_fallback=radius_fallback):
        records.append(dict(row, pose_id=row["pose_i"], detail=detail))
    labels = label_reference_results(records, thresholds)
    write_json(output / "reference_details.json", records)
    write_csv(output / "reference_evaluation.csv", labels)
    write_csv(output / "reference_summary.csv", summarize_reference(scores, labels))
    return labels


def label_reference_results(records, thresholds):
    """Success means RMSD < threshold; a failed or ambiguous comparison is unknown, not a failure."""
    labels = []
    for record in records:
        value = record.get("rmsd")
        known = not record.get("error") and value is not None and math.isfinite(value)
        for threshold in thresholds:
            labels.append(dict(
                ensemble=record["ensemble"],
                pose_id=record["pose_id"],
                threshold_angstrom=threshold,
                inequality="rmsd < threshold",
                rmsd=value,
                success=bool(value < threshold) if known else None,
                radius_fallback_used=record.get("radius_fallback_used", False),
            ))
    return labels


def summarize_reference(scores, labels):
    """Each ensemble's CORE Score beside the fraction of its poses within the same threshold.

    Descriptive only; never fit a CORE cutoff on these same observations. Labels describe
    whole poses, so every site of a multi-copy ensemble shares them.
    """
    summaries = []
    for score in scores:
        if score["core_score"] is None:
            continue
        success = [label["success"] for label in labels if label["ensemble"] == score["ensemble"]
                   and label["threshold_angstrom"] == score["threshold_angstrom"]]
        known = [value for value in success if value is not None]
        summaries.append(dict(
            ensemble=score["ensemble"],
            site=score["site"],
            threshold_angstrom=score["threshold_angstrom"],
            core_score=score["core_score"],
            pose_count=len(success),
            evaluated_count=len(known),
            correct_count=sum(known),
            correct_fraction=sum(known) / len(success) if success and len(known) == len(success) else None,
        ))
    return summaries
