"""Optional comparison of every prediction against its ensemble's ground-truth structure."""

from core_score.results import fingerprint, write_csv, write_json
from core_score.stages.preparation import prepare
from core_score.stages.utils.comparison_utils import run_pair_tasks
from core_score.stages.utils.evaluation_utils import label_reference_results


def evaluate_groups(groups, output, thresholds=(2.0,), *, workers=1, timeout=None, radius_fallback=True):
    """Ligand RMSD of each prepared prediction to the reference, labelled at each threshold."""
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
    write_csv(output / "reference_summary.csv", summaries)
    return summaries
