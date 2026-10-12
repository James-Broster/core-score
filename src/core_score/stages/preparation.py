"""Prepare every pose's structure and ligand chemistry, and check OST can read the ligand."""

from core_score.results import diagnostic, write_json
from core_score.stages.utils import ost
from core_score.stages.utils.preparation.prepare import pose_directory, prepare_pose


def prepare_groups(groups, output, *, timeout=None):
    cache = {}
    prepared = []
    for group in groups:
        poses = [prepare(pose, pose_directory(output, pose) / "pose.cif", cache, timeout) for pose in group["poses"]]
        # A ligand graph that failed can borrow one from a peer prepared later in the same ensemble.
        for i, pose in enumerate(poses):
            if "prepared" in pose and any(copy["pb_error"] for copy in pose["prepared"]["copy_records"]):
                retry = prepare(group["poses"][i], pose_directory(output, pose) / "pose.cif", cache, timeout)
                if "prepared" in retry:
                    poses[i] = retry
        prepared.append(dict(group, poses=poses))
    write_json(output / "prepared_manifest.json", dict(groups=prepared))
    return prepared


def prepare(pose, destination, cache=None, timeout=None):
    """The pose with its prepared record, or with the error that stopped preparation."""
    try:
        prepared = prepare_pose(pose, destination, cache)
        if not prepared["ligand_sdfs"] and not ost.detects_ligands(prepared, timeout):
            # OST cannot read this ligand from the structure, so it gets the connectivity SDFs.
            prepared.update(ligand_sdfs=[copy["ost_sdf"] for copy in prepared["copy_records"]],
                            route="explicit_connectivity_sdf")
        return dict(pose, prepared=prepared)
    except Exception as exc:
        return dict(pose, error=diagnostic("preparation", exc, structure=pose["path"]))
