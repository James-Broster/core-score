"""Run PoseBusters dock checks on every pose and save the results."""

from core_score.results import write_csv, write_json
from core_score.stages.utils.posebusters_utils import check_pose, dock_checks


def run_pb(groups, output, skip=False):
    """Exclude poses failing a dock check; an unknown result blocks the ensemble."""
    if skip:
        return groups
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock", max_workers=0)
    checks = dock_checks(buster)
    records, filtered = [], []
    for group in groups:
        results = [check_pose(buster, checks, pose) for pose in group["poses"]]
        records += results
        kept = [pose for pose, result in zip(group["poses"], results) if result["pb_valid"] is not False]
        filtered.append(dict(group, poses=kept))
        if any(result["pb_valid"] is None for result in results):
            filtered[-1]["error"] = dict(stage="posebusters", error="One or more poses have unknown PoseBusters status")
    write_json(output / "posebusters.json", records)
    write_csv(output / "posebusters.csv", [
        dict(ensemble=r["ensemble"], pose_id=r["pose_id"], pb_valid=r["pb_valid"],
             failed_checks=sorted({c for copy in r.get("copies", []) for c in copy.get("failed_checks", [])}),
             error=r.get("error") or next((c["error"] for c in r.get("copies", []) if c.get("error")), None))
        for r in records
    ])
    return filtered
