"""Prepare reference comparisons and interpret their measurements."""

from pathlib import Path
import numpy as np
from core_score.stages.utils.comparison_utils import probe_prepared
from core_score.stages.utils.prepare_utils import prepare_pose
from core_score.results import diagnostic, fingerprint, write_csv


def prepare_reference_tasks(groups, output, *, ost="ost", chemlib=None, timeout=None):
    records = []
    tasks = []
    for group in groups:
        if not group.get("reference"):
            raise ValueError(
                f'Reference evaluation requested but no reference supplied for {group["ensemble"]}'
            )
        reference = group["reference"]
        try:
            dest = (
                Path(output)
                / "references"
                / fingerprint(group["ensemble"])[:20]
                / "reference.cif"
            )
            prepared = prepare_pose(reference, dest, chemistry=True)
            prepared = probe_prepared(prepared, ost=ost, chemlib=chemlib, timeout=timeout)
        except Exception as exc:
            records.extend(
                dict(
                    ensemble=group["ensemble"],
                    pose_id=p["pose_id"],
                    rmsd=None,
                    error=diagnostic("reference_preparation", exc),
                )
                for p in group["poses"]
            )
            continue
        for pose in group["poses"]:
            try:
                if "prepared" not in pose:
                    raise ValueError("Prediction failed preparation")
                model = probe_prepared(pose["prepared"], ost=ost, chemlib=chemlib, timeout=timeout)
                # Multiple reference copies use OST's native assignment; scalar correctness
                # remains undefined unless exactly one pair was assigned.
                tasks.append(
                    (
                        dict(group, site="reference"),
                        pose,
                        reference,
                        model,
                        prepared,
                    )
                )
            except Exception as exc:
                records.append(
                    dict(
                        ensemble=group["ensemble"],
                        pose_id=pose["pose_id"],
                        rmsd=None,
                        error=diagnostic("reference_comparison", exc),
                    )
                )

    return tasks, records


def label_reference_results(records, thresholds):
    labels = []
    for record in records:
        for threshold in thresholds:
            value = record.get("rmsd")
            known = not record.get("error") and value is not None and np.isfinite(value)
            labels.append(
                dict(
                    ensemble=record["ensemble"],
                    pose_id=record["pose_id"],
                    threshold_angstrom=threshold,
                    inequality="rmsd < threshold",
                    rmsd=value,
                    success=bool(value < threshold) if known else None,
                    radius_fallback_used=record.get("radius_fallback_used", False),
                )
            )
    return labels


def summarize_reference(scores, labels, output):
    """Descriptive evaluation only; never fit a CORE cutoff on these same observations."""
    summaries = []
    for score in scores:
        if score["status"] != "complete":
            continue
        for threshold in sorted({r["threshold_angstrom"] for r in labels}):
            matching = {
                r["pose_id"]: r["success"]
                for r in labels
                if r["ensemble"] == score["ensemble"]
                and r["threshold_angstrom"] == threshold
            }
            members = score["group_pose_ids"]
            known = [matching.get(p) for p in members]
            complete = all(v is not None for v in known)
            summaries.append(
                dict(
                    ensemble=score["ensemble"],
                    site=score["site"],
                    core_score=score["S"],
                    threshold_angstrom=threshold,
                    inequality="rmsd < threshold",
                    winning_group_size=len(members),
                    evaluated_members=sum(v is not None for v in known),
                    winning_group_success_fraction=(
                        sum(known) / len(known) if complete and known else None
                    ),
                    note="Reference assignment is OST-native; multicopy labels describe whole poses, not site-specific correctness",
                )
            )
    write_csv(Path(output) / "reference_summary.csv", summaries)
    return summaries
