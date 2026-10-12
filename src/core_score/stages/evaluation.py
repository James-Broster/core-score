"""Optional comparison of every prediction against its ensemble's ground-truth structure."""

from core_score.results import fingerprint, write_csv, write_json
from core_score.stages.preparation import prepare
from core_score.stages.utils.comparison_utils import run_pair_tasks
from core_score.stages.utils.evaluation_utils import label_reference_results, summarize_reference


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
