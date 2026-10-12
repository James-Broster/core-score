"""Align ligand-copy centres and write per-site prepared structures."""

from pathlib import Path
import numpy as np
from core_score.stages.utils.preparation.cif import read_document, filter_records
from core_score.stages.utils.site_algorithms import rigid_transform, apply_transform


def align_copy_centres(poses, copies, ost, timeout):
    anchor = poses[0]
    points, transforms = {}, {}
    for pose in poses:
        prepared = pose["prepared"]
        if len(prepared["copy_records"]) != copies:
            raise ValueError("Missing copy coordinates for site correspondence")
        if pose is anchor:
            transform = np.eye(4)
        else:
            fit = rigid_transform(
                ost,
                Path(prepared["path"]),
                Path(anchor["prepared"]["path"]),
                timeout,
            )
            transform = fit.transform
            transforms[pose["pose_id"]] = dict(
                rmsd=fit.rmsd,
                chain_mapping=fit.chain_mapping,
                transform=transform,
            )
        coordinates = apply_transform(
            [c["centre"] for c in prepared["copy_records"]], transform
        )
        points[pose["pose_id"]] = {
            c["copy_id"]: xyz
            for c, xyz in zip(prepared["copy_records"], coordinates)
        }
    return points, transforms


def prepare_site_pose(pose, copy, destination):
    prepared = pose["prepared"]
    if copy.get("residue"):
        filter_records(
            read_document(prepared["path"]), {tuple(copy["residue"])}, destination
        )
        site_prepared = dict(
            prepared,
            path=str(destination.resolve()),
            copies=1,
            ligand_sdfs=[],
            copy_records=[copy],
        )
    else:
        site_prepared = dict(
            prepared,
            copies=1,
            ligand_sdfs=[copy["ost_sdf"]],
            copy_records=[copy],
        )
    return dict(pose, prepared=site_prepared)
