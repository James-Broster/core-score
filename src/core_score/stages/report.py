"""Write the files a user reads: every complex's CORE Score with its verdict at each cutoff,
a summary per cutoff, PoseBusters validity for every pose and, with references, every
pose's result against its reference. Everything else is in details/."""

from core_score.results import write_csv
from core_score.stages.evaluation import correct_column
from core_score.stages.posebusters import summary_row


def write_report(groups, scores, pb_records, labels, output, cutoffs=(0.5, 0.6, 0.7, 0.8, 0.9)):
    validity = None
    if pb_records is not None:
        validity = {(r["ensemble"], r["pose_id"]): r["pb_valid"] for r in pb_records}
    ranked = any("confidence" in group["poses"][0] for group in groups)
    rows = core_table(groups, scores, validity, labels, cutoffs, ranked)
    write_csv(output / "core_scores.csv", rows)
    judged = bool(labels) and ranked and validity is not None  # top poses exist and can be checked
    write_csv(output / "cutoff_summary.csv", summarize_cutoffs(rows, cutoffs, judged))
    if pb_records is not None:
        write_csv(output / "posebusters.csv", [summary_row(record) for record in pb_records])
    if labels:
        write_csv(output / "reference_evaluation.csv", labels)
    return rows


def top_pose(group, validity):
    """The PoseBusters-valid pose with the highest confidence, the paper's top-ranked prediction.

    Returns the pose, or None with the reason there is none. `validity` maps (ensemble,
    pose ID) to pb_valid, or is None when PoseBusters was skipped. Ties go to the earlier pose.
    """
    if "confidence" not in group["poses"][0]:
        return None, "no confidences given"
    if validity is None:
        return None, "PoseBusters was skipped"
    valid = [pose for pose in group["poses"] if validity.get((pose["ensemble"], pose["pose_id"])) is True]
    if not valid:
        return None, "no PoseBusters-valid pose"
    return max(valid, key=lambda pose: pose["confidence"]), None


def core_table(groups, scores, validity, labels, cutoffs, ranked):
    """One row per complex (ensemble and site) and threshold, highest CORE first.

    `above_<cutoff>` is True when CORE >= cutoff. Top-pose columns appear when confidences
    were given, and reference columns when references were evaluated.
    """
    by_name = {group["ensemble"]: group for group in groups}
    by_pose = {(label["ensemble"], label["pose_id"]): label for label in labels}

    def correct(pose, threshold):
        return by_pose.get((pose["ensemble"], pose["pose_id"]), {}).get(correct_column(threshold))

    rows = []
    for score in scores:
        group, core, threshold = by_name[score["ensemble"]], score["core_score"], score["threshold_angstrom"]
        row = dict(ensemble=score["ensemble"], site=score["site"], threshold_angstrom=threshold, core_score=core)
        for cutoff in cutoffs:
            row[f"above_{cutoff:g}"] = None if core is None else core >= cutoff
        row["pose_count"] = score["pose_count"]
        if validity is not None:
            row["valid_pose_count"] = sum(validity.get((p["ensemble"], p["pose_id"])) is True for p in group["poses"])
        note = score["error"]
        if ranked:
            pose, reason = top_pose(group, validity)
            row.update(top_pose=pose["pose_id"] if pose else None, top_pose_confidence=pose["confidence"] if pose else None)
            note = note or reason
        if labels:
            known = [correct(p, threshold) for p in group["poses"]]
            row["fraction_correct"] = sum(known) / len(known) if all(k is not None for k in known) else None
            if ranked:
                row["top_pose_correct"] = correct(pose, threshold) if pose else None
        row["note"] = note
        rows.append(row)
    rows.sort(key=lambda r: (r["threshold_angstrom"], r["core_score"] is None, -(r["core_score"] or 0)))
    return rows


def summarize_cutoffs(rows, cutoffs, judged=False):
    """Per threshold and cutoff: how many complexes are above and below it and, with references
    and confidences, how often the top pose is correct on each side."""
    summaries = []
    for threshold in dict.fromkeys(row["threshold_angstrom"] for row in rows):
        scored = [row for row in rows if row["threshold_angstrom"] == threshold and row["core_score"] is not None]
        for cutoff in cutoffs:
            above = [row for row in scored if row[f"above_{cutoff:g}"]]
            below = [row for row in scored if not row[f"above_{cutoff:g}"]]
            summary = dict(threshold_angstrom=threshold, core_cutoff=cutoff, complexes=len(scored),
                           above=len(above), below=len(below),
                           fraction_above=len(above) / len(scored) if scored else None)
            if judged:
                for side, side_rows in (("above", above), ("below", below)):
                    known = [row["top_pose_correct"] for row in side_rows if row["top_pose_correct"] is not None]
                    summary[f"top_poses_known_{side}"] = len(known)
                    summary[f"top_poses_correct_{side}"] = sum(known)
                    summary[f"accuracy_{side}"] = sum(known) / len(known) if known else None
            summaries.append(summary)
    return summaries
