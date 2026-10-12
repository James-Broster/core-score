"""Turn CORE Scores into decisions: which complexes CORE trusts at each cutoff, and, when
confidences are given, the top-ranked pose of each and whether it is correct."""

from core_score.results import write_csv
from core_score.stages.utils.decision_utils import decide, summarize_cutoffs


def decide_groups(groups, scores, pb_records, labels, output, cutoffs=(0.5, 0.6, 0.7, 0.8, 0.9)):
    validity = None
    if pb_records is not None:
        validity = {}
        for record in pb_records:
            validity.setdefault(record["ensemble"], {})[record["pose_id"]] = record["pb_valid"]
    decisions = decide(groups, scores, validity, labels, cutoffs)
    write_csv(output / "decisions.csv", decisions)
    write_csv(output / "cutoff_summary.csv", summarize_cutoffs(decisions))
    return decisions
