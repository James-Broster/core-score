"""OST rigid superposition and the historical k-means/Hungarian site assignment."""

import dataclasses
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.optimize import linear_sum_assignment


class RigidSuperpositionError(RuntimeError):
    """Raised when `ost compare-structures --rigid-scores` cannot superpose two poses."""

    def __init__(self, message, evidence=None):
        super().__init__(message)
        self.evidence = evidence or {}


@dataclasses.dataclass(frozen=True)
class RigidResult:
    transform: (
        np.ndarray
    )  # (4, 4) homogeneous transform: superposes `model` onto `reference`
    rmsd: float  # protein backbone RMSD after superposition (informational)
    chain_mapping: dict  # OST's own reference-chain -> model-chain mapping


def rigid_transform(
    ost_bin: str, model_cif: Path, reference_cif: Path, timeout: float | None = None
) -> RigidResult:
    """Superpose `model_cif` onto `reference_cif` using OST's own rigid Kabsch superposition
    of matched protein C-alpha atoms (`ost compare-structures --rigid-scores`).

    Raises RigidSuperpositionError on any failure (process error, non-SUCCESS status, or a
    missing transform) rather than returning a sentinel -- callers must handle this
    explicitly, not have it silently absorbed.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_json = f"{tmp}/rigid.json"
        cmd = [
            ost_bin,
            "compare-structures",
            "-m",
            str(model_cif),
            "-r",
            str(reference_cif),
            "--rigid-scores",
            "-o",
            out_json,
        ]
        evidence = dict(
            command=cmd,
            timeout=timeout,
            model=str(model_cif),
            reference=str(reference_cif),
        )
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise RigidSuperpositionError(
                f"ost compare-structures timed out after {timeout}s", evidence
            ) from exc
        evidence.update(
            returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr
        )
        if proc.returncode != 0:
            raise RigidSuperpositionError(
                proc.stderr.strip()
                or proc.stdout.strip()
                or f"exit code {proc.returncode}",
                evidence,
            )
        data = json.loads(Path(out_json).read_text())
        if data.get("status") != "SUCCESS":
            raise RigidSuperpositionError(
                f"OST status={data.get('status')}", dict(evidence, ost_results=data)
            )
        transform = data.get("transform")
        if transform is None:
            raise RigidSuperpositionError(
                "OST returned no transform (chains likely could not be mapped)",
                dict(evidence, ost_results=data),
            )
        return RigidResult(
            transform=np.asarray(transform, dtype=float),
            rmsd=data.get("rmsd"),
            chain_mapping=data.get("chain_mapping") or {},
        )


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
