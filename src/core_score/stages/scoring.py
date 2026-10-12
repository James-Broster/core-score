"""Build each ensemble matrix, calculate CORE, and export score summaries."""

from pathlib import Path
from core_score.stages.utils.score_utils import build_matrix, core_score
from core_score.results import diagnostic, write_json, write_csv


def score_groups(groups, rows, output):
    scores = []
    for group in groups:
        ids = sorted(p["pose_id"] for p in group["poses"])
        result = dict(
            ensemble=group["ensemble"],
            site=group["site"],
            original_pose_count=group.get("original_pose_count", len(ids)),
            retained_pose_count=len(ids),
            filter_pb_invalid=group.get("filter_pb_invalid", False),
        )
        try:
            if group.get("stage_error"):
                raise ValueError(group["stage_error"]["error"])
            selected = [
                r
                for r in rows
                if r["ensemble"] == group["ensemble"]
                and str(r["site"]) == str(group["site"])
            ]
            matrix = build_matrix(ids, selected)
            score = core_score(matrix)
            result.update(
                score,
                status="complete",
                group_pose_ids=[ids[i] for i in score["group_members"]],
                pose_ids=ids,
            )
        except Exception as exc:
            result.update(
                status="incomplete", diagnostic=diagnostic("matrix_or_score", exc)
            )
        scores.append(result)
    write_json(Path(output) / "scores.json", scores)
    write_csv(Path(output) / "scores.csv", scores)
    return scores
