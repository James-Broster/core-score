"""Build each site ensemble's RMSD matrix, calculate CORE at each threshold, and save the scores."""

from core_score.results import diagnostic, write_csv, write_json
from core_score.stages.utils.score_utils import build_matrix, core_score


def score_groups(groups, rows, output, thresholds=(2.0,)):
    """One score per site ensemble and threshold T. An ensemble with a missing or failed pair is not scored."""
    scores = []
    for group in groups:
        ids = [pose["pose_id"] for pose in group["poses"]]
        values, error = dict.fromkeys(thresholds), group.get("error")
        if not error:
            try:
                key = (group["ensemble"], group["site"])
                matrix = build_matrix(ids, [r for r in rows if (r["ensemble"], r["site"]) == key])
                values = {threshold: core_score(matrix, threshold) for threshold in thresholds}
            except ValueError as exc:
                error = diagnostic("scoring", exc)
        scores += [
            dict(ensemble=group["ensemble"], site=group["site"], pose_count=len(ids), threshold_angstrom=threshold,
                 core_score=value, error=error["error"] if error else None)
            for threshold, value in values.items()
        ]
    write_json(output / "scores.json", scores)
    write_csv(output / "scores.csv", scores)
    return scores
