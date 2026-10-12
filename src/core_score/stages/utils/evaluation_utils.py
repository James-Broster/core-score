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
