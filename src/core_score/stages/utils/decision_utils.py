"""Whether CORE accepts each ensemble at each cutoff, and each ensemble's top-ranked pose."""


def top_pose(group, validity):
    """The PoseBusters-valid pose with the highest confidence, the paper's top-ranked prediction.

    Returns the pose, or None with the reason there is none. `validity` maps pose ID to
    pb_valid, or is None when PoseBusters was skipped. Ties go to the earlier pose.
    """
    if "confidence" not in group["poses"][0]:
        return None, "no confidences given"
    if validity is None:
        return None, "PoseBusters was skipped"
    valid = [pose for pose in group["poses"] if validity.get(pose["pose_id"]) is True]
    if not valid:
        return None, "no PoseBusters-valid pose"
    return max(valid, key=lambda pose: pose["confidence"]), None


def decide(groups, scores, validity, labels, cutoffs):
    """One row per scored site ensemble, threshold and CORE cutoff."""
    tops = {g["ensemble"]: top_pose(g, None if validity is None else validity.get(g["ensemble"], {})) for g in groups}
    success = {(l["ensemble"], l["pose_id"], l["threshold_angstrom"]): l["success"] for l in labels}
    rows = []
    for score in scores:
        pose, reason = tops[score["ensemble"]]
        for cutoff in cutoffs:
            rows.append(dict(
                ensemble=score["ensemble"],
                site=score["site"],
                threshold_angstrom=score["threshold_angstrom"],
                core_cutoff=cutoff,
                core_score=score["core_score"],
                likely_correct=None if score["core_score"] is None else score["core_score"] >= cutoff,
                top_pose=pose["pose_id"] if pose else None,
                top_pose_confidence=pose["confidence"] if pose else None,
                top_pose_correct=success.get((score["ensemble"], pose["pose_id"], score["threshold_angstrom"])) if pose else None,
                note=score["error"] or reason,
            ))
    return rows


def summarize_cutoffs(decisions):
    """Per threshold and cutoff: how many ensembles CORE accepts and, with references, how
    many accepted top-ranked poses are correct."""
    summaries = []
    for key in dict.fromkeys((d["threshold_angstrom"], d["core_cutoff"]) for d in decisions):
        rows = [d for d in decisions if (d["threshold_angstrom"], d["core_cutoff"]) == key]
        scored = [d for d in rows if d["likely_correct"] is not None]
        accepted = [d for d in scored if d["likely_correct"]]
        judged = [d["top_pose_correct"] for d in accepted if d["top_pose_correct"] is not None]
        summaries.append(dict(
            threshold_angstrom=key[0],
            core_cutoff=key[1],
            scored=len(scored),
            accepted=len(accepted),
            coverage=len(accepted) / len(scored) if scored else None,
            accepted_with_known_top_pose=len(judged),
            accepted_top_pose_correct=sum(judged),
            accuracy=sum(judged) / len(judged) if judged else None,
        ))
    return summaries
