"""Prepare reference comparisons and interpret their measurements."""

from pathlib import Path
import numpy as np
from core_score.stages.utils.comparison_utils import probe_prepared
from core_score.stages.utils.preparation.prepare import prepare_pose
from core_score.results import diagnostic, fingerprint


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
            prepared = prepare_pose(reference, dest)
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
