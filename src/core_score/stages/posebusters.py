"""Orchestrate PoseBusters checking, filtering and result export."""

from pathlib import Path
from core_score.stages.utils.posebusters_utils import load_results, pose_result, filter_group
from core_score.results import write_json, write_csv


def run_pb(groups, output, skip=False, filter_invalid=False, supplied=None):
    if skip and (filter_invalid or supplied):
        raise ValueError("--skip-posebusters cannot be combined with filtering or supplied PB results")
    expected = [(g["ensemble"], p["pose_id"]) for g in groups for p in g["poses"]]
    supplied_results = load_results(supplied, expected) if supplied else None
    records, filtered = [], []
    for group in groups:
        group_results = [
            pose_result(pose, group["ensemble"], skip=skip, supplied=supplied_results)
            for pose in group["poses"]
        ]
        filtered.append(filter_group(group, group_results, filter_invalid))
        records.extend(group_results)
    write_json(Path(output) / "posebusters.json", records)
    write_csv(Path(output) / "posebusters.csv", records)
    return filtered
