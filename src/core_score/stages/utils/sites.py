"""Give each ligand copy a site label that is consistent across an ensemble."""

import dataclasses
from pathlib import Path
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.optimize import linear_sum_assignment
from core_score.stages.utils.preparation.cif import read_document, filter_records
from core_score.stages.utils.ost import rigid_transform


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


def apply_transform(coords: np.ndarray, transform: np.ndarray) -> np.ndarray:
    """Apply a 4x4 homogeneous transform to an (N, 3) array of coordinates. Pure arithmetic:
    the transform itself always comes from rigid_transform(), i.e. from OST."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(
            f"expected an (N, 3) coordinate array, got shape {coords.shape}"
        )
    homogeneous = np.hstack([coords, np.ones((len(coords), 1))])
    return (homogeneous @ transform.T)[:, :3]


@dataclasses.dataclass(frozen=True)
class SiteAssignment:
    """One pose's copies, each pinned to exactly one site label (0..k-1)."""

    pose_id: object
    chain_to_site: dict  # {chain_label: site_index}
    chain_to_cost: dict  # {chain_label: squared distance to that site's centroid}


def assign_sites(
    pose_points: dict, k: int, seed: int = 0, max_iter: int = 50
) -> tuple[list[SiteAssignment], np.ndarray]:
    """Assign each pose's k ligand-copy coordinates to k persistent site labels, consistent
    across all poses, with an exact one-copy-per-pose-per-site guarantee.

    pose_points: {pose_id: {chain_label: xyz}} -- xyz already in one shared coordinate frame
                 (i.e. already passed through apply_transform against a common reference).
                 Every pose must supply exactly k points; a pose with a different count is
                 the caller's responsibility to exclude beforehand and report, not silently
                 drop here (this function raises if that invariant is violated).
    k:           the group's expected number of identical copies.

    Returns (assignments, final_centroids). The site *labels* (0..k-1) are arbitrary --
    only their consistency across poses is meaningful (see site_label_permutation_distance
    for comparing two labellings up to a global relabelling).

    Method: k-means gives an initial guess for the k centroids (not a final assignment --
    plain k-means has no per-pose one-to-one constraint and is not trusted for that). Each
    pose's points are then optimally assigned to those k centroids one-to-one via the
    Hungarian algorithm (scipy.optimize.linear_sum_assignment, cost = squared distance),
    which guarantees exactly one copy per pose per site by construction, however displaced
    a copy is -- there is no distance threshold that drops a point. Centroids are then
    re-estimated from the new assignment and the process repeats to convergence (or
    max_iter), a fixed-point refinement of k-means-with-a-hard-per-pose-constraint.
    """
    pose_ids = list(pose_points)
    for pid in pose_ids:
        n = len(pose_points[pid])
        if n != k:
            raise ValueError(
                f"pose {pid!r} has {n} copies, expected exactly k={k}; caller must exclude/report this pose before calling assign_sites"
            )

    all_points = np.vstack(
        [np.asarray(list(pose_points[pid].values())) for pid in pose_ids]
    )
    centroids, _ = kmeans2(all_points, k, seed=seed, minit="++")

    assignment_of = (
        {}
    )  # pose_id -> (ordered chain labels, site index per chain, cost per chain)
    for _ in range(max_iter):
        changed = False
        for pid in pose_ids:
            chains = list(pose_points[pid].keys())
            pts = np.asarray([pose_points[pid][c] for c in chains])
            cost = ((pts[:, None, :] - centroids[None, :, :]) ** 2).sum(
                axis=2
            )  # (k, k)
            row_ind, col_ind = linear_sum_assignment(cost)
            new = tuple(col_ind[row_ind.argsort()])
            if assignment_of.get(pid, (None,))[1:2] != (new,):
                changed = True
            assignment_of[pid] = (chains, new, cost[row_ind, col_ind])
        if not changed and _ > 0:
            break
        for s in range(k):
            pts_s = []
            for pid in pose_ids:
                chains, sites, _ = assignment_of[pid]
                for c, site in zip(chains, sites):
                    if site == s:
                        pts_s.append(pose_points[pid][c])
            if pts_s:
                centroids[s] = np.mean(pts_s, axis=0)
    else:
        raise ValueError(
            "Site correspondence did not converge; no site assignments accepted"
        )

    results = []
    for pid in pose_ids:
        chains, sites, costs = assignment_of[pid]
        results.append(
            SiteAssignment(
                pose_id=pid,
                chain_to_site=dict(zip(chains, (int(s) for s in sites))),
                chain_to_cost=dict(zip(chains, (float(c) for c in costs))),
            )
        )
    return results, centroids
