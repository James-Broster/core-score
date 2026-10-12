"""Write the two files a user reads: core_scores.csv, one row per complex, and poses.csv,
one row per pose. Everything behind them is in details/."""

from core_score.results import write_csv
from core_score.stages.evaluation import correct_column


def write_report(groups, scores, pb_records, labels, output, cutoffs=(0.5, 0.6, 0.7, 0.8, 0.9)):
    validity = None
    if pb_records is not None:
        validity = {(r["ensemble"], r["pose_id"]): r["pb_valid"] for r in pb_records}
    ranked = any("confidence" in group["poses"][0] for group in groups)
    rows = core_table(groups, scores, validity, labels, cutoffs, ranked)
    write_csv(output / "core_scores.csv", rows)
    write_csv(output / "poses.csv", pose_table(groups, pb_records, labels, ranked))
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



def pose_table(groups, pb_records, labels, ranked):
    """One row per pose: its confidence, whether PoseBusters passed it and, with references,
    its RMSD to the reference.

    Columns that do not apply are left out, and a blank is unknown. The reference columns
    are whichever `label_reference_results()` produced, so there is one per threshold.
    """
    checks = {(r["ensemble"], r["pose_id"]): r for r in pb_records or []}
    measured = {(label["ensemble"], label["pose_id"]): label for label in labels}
    rows = []
    for group in groups:
        for pose in group["poses"]:
            key = (pose["ensemble"], pose["pose_id"])
            row = dict(ensemble=pose["ensemble"], pose_id=pose["pose_id"])
            if ranked:
                row["confidence"] = pose.get("confidence")
            record = checks.get(key)
            if pb_records is not None:
                # True, False, or blank when unknown. Which checks failed is in details/posebusters.json.
                row["pb_valid"] = record["pb_valid"] if record else None
            label = measured.get(key, {})
            for name, value in label.items():
                if name not in ("ensemble", "pose_id", "error", "radius_fallback_used"):
                    row[name] = value
            row["note"] = message(pose.get("error")) or (record or {}).get("error") or message(label.get("error"))
            rows.append(row)
    return rows


def message(error):
    """The one-line message of a diagnostic record, or of a plain error string."""
    return (error["error"] if isinstance(error, dict) else error) or None
