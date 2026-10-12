"""PoseBusters dock checks for every pose. The results are saved for use downstream, such as
choosing the top-ranked valid pose. CORE uses every pose whatever its result."""

from core_score.results import write_csv, write_json
from core_score.stages.utils.posebusters_utils import check_pose, dock_checks, summary_row


def run_posebusters(groups, output):
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock", max_workers=0)
    checks = dock_checks(buster)
    records = [check_pose(buster, checks, pose) for group in groups for pose in group["poses"]]
    write_json(output / "posebusters.json", records)
    write_csv(output / "posebusters.csv", [summary_row(record) for record in records])
    return records
