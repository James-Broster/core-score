"""Compare every pair of poses within each site ensemble and save the measurements."""

from core_score.results import write_csv, write_jsonl
from core_score.stages.utils.comparison_utils import pair_tasks, run_pair_tasks


def compare_groups(groups, output, *, workers=1, timeout=None, radius_fallback=True):
    results = run_pair_tasks(pair_tasks(groups), workers, timeout=timeout, radius_fallback=radius_fallback)
    results.sort(key=lambda result: tuple(result[0][k] for k in ("ensemble", "site", "pose_i", "pose_j")))
    rows = [row for row, _ in results]
    write_csv(output / "pairs.csv", rows)
    write_jsonl(output / "ligand_details.jsonl", [detail for _, detail in results])
    return rows
