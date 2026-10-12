"""Optional comparison of every prediction against its ensemble's ground-truth structure."""

import math
from core_score.results import fingerprint, write_json
from core_score.stages.comparison import run_pair_tasks
from core_score.stages.preparation import prepare


def evaluate_groups(groups, output, thresholds=(2.0,), *, workers=1, timeout=None, radius_fallback=True):
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
    write_json(output / "reference_details.json", records)
    keys = ((group["ensemble"], pose["pose_id"]) for group in groups for pose in group["poses"])
    order = {key: i for i, key in enumerate(keys)}
    labels = label_reference_results(records, thresholds)
    return sorted(labels, key=lambda label: order[label["ensemble"], label["pose_id"]])


def correct_column(threshold):
    """The reference_evaluation.csv column saying whether a pose's RMSD is below `threshold`."""
    return f"rmsd_below_{threshold:g}A"


def label_reference_results(records, thresholds):
    """One row per pose. A pose is correct at T when RMSD < T. A failed or ambiguous comparison
    is unknown (blank RMSD and labels), not a failure."""
    labels = []
    for record in records:
        value = record.get("rmsd")
        known = not record.get("error") and value is not None and math.isfinite(value)
        label = dict(ensemble=record["ensemble"], pose_id=record["pose_id"], rmsd=value if known else None,
                     lddt_pli=record.get("lddt_pli") if known else None)
        for threshold in thresholds:
            label[correct_column(threshold)] = bool(value < threshold) if known else None
        label["radius_fallback_used"] = record.get("radius_fallback_used", False)
        error = record.get("error")
        label["error"] = (error["error"] if isinstance(error, dict) else error) or None  # message only; details in JSON
        labels.append(label)
    return labels
