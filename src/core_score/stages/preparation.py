"""Orchestrate structure/chemistry preparation for each input ensemble."""

from pathlib import Path
from core_score.stages.utils.prepare_utils import prepare_pose
from core_score.results import diagnostic, fingerprint, write_json


def prepare_groups(groups, output):
    output = Path(output)
    cache = {}
    prepared_groups = []
    for group in groups:
        poses = []
        for pose in group["poses"]:
            if pose.get("prepared"):
                poses.append(pose)
                continue
            record = dict(pose)
            destination = (
                output
                / "prepared"
                / fingerprint((group["ensemble"], pose["pose_id"]))[:20]
                / "pose.cif"
            )
            try:
                record["prepared"] = prepare_pose(
                    pose, destination, chemistry=True, cache=cache
                )
            except Exception as exc:
                record["preparation_error"] = diagnostic(
                    "preparation",
                    exc,
                    structure=pose["path"],
                    ensemble=group["ensemble"],
                    pose_id=pose["pose_id"],
                )
            poses.append(record)
        # Retry chemistry-only failures once after peer graphs have become available.
        for i, pose in enumerate(poses):
            prepared = pose.get("prepared")
            if (
                prepared
                and any(c.get("pb_error") for c in prepared.get("copy_records", []))
                and not pose.get("ligand_sdfs")
            ):
                try:
                    poses[i] = dict(
                        pose,
                        prepared=prepare_pose(
                            pose, Path(prepared["path"]), chemistry=True, cache=cache
                        ),
                    )
                except Exception:
                    pass  # First preparation and its original diagnostic remain authoritative.
        prepared_groups.append(
            dict(
                group,
                poses=poses,
                original_pose_count=group.get("original_pose_count", len(poses)),
            )
        )
    write_json(
        output / "prepared_manifest.json",
        dict(schema="core-score.prepared.v1", groups=prepared_groups),
    )
    return prepared_groups

