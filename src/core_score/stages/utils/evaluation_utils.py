"""Turn reference RMSDs into success labels."""

import math


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
