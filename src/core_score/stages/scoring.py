"""Build each site ensemble's RMSD matrix, calculate CORE, and save the scores."""

from core_score.results import diagnostic, write_csv, write_json
from core_score.stages.utils.score_utils import build_matrix, core_score


def score_groups(groups, rows, output):
    scores = []
    for group in groups:
        ids = sorted(pose["pose_id"] for pose in group["poses"])
        result = dict(ensemble=group["ensemble"], site=group["site"], pose_count=len(ids))
        error = group.get("error")
        if not error:
            try:
                selected = [r for r in rows if (r["ensemble"], r["site"]) == (group["ensemble"], group["site"])]
                score = core_score(build_matrix(ids, selected))
                result.update(score, status="complete", group_pose_ids=[ids[i] for i in score["group_members"]], pose_ids=ids)
            except ValueError as exc:
                error = diagnostic("scoring", exc)
        if error:
            result.update(status="incomplete", error=error["error"])
        scores.append(result)
    write_json(output / "scores.json", scores)
    write_csv(output / "scores.csv", scores)
    return scores
