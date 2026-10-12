"""Prediction-only site grouping, separate from OST scoring.

Protein atom correspondence must be supplied explicitly. These functions never
infer chemical identity, change structures, or choose an OST chain mapping.
Stability of the output does not prove a biological correspondence.
"""

from dataclasses import dataclass

import numpy as np


def _vectors(value):
    array = np.asarray(value, dtype=float)
    if array.ndim != 2 or array.shape[1] < 1 or not np.isfinite(array).all():
        raise ValueError("expected finite N x D vectors")
    return array


def _points(value):
    array = _vectors(value)
    if array.shape[1] != 3:
        raise ValueError("expected finite N x 3 coordinates")
    return array


def protein_transform(moving, fixed):
    """OST SuperposeSVD fit of explicitly corresponding protein atoms.

    Row-vector convention: transformed = coordinates @ rotation + translation.
    No ligand coordinates participate in the fit. Inputs are never modified.
    Requires the OST Python environment; there is no custom fitting fallback.
    """
    moving, fixed = _points(moving), _points(fixed)
    if moving.shape != fixed.shape or len(moving) < 3:
        raise ValueError("at least three corresponding protein atoms required")
    x, y = moving - moving.mean(0), fixed - fixed.mean(0)
    if min(np.linalg.matrix_rank(x), np.linalg.matrix_rank(y)) < 2:
        raise ValueError("collinear protein atoms do not determine a rotation")
    from ost import geom
    from ost.mol.alg import SuperposeSVD

    model, reference = geom.Vec3List(), geom.Vec3List()
    for a, b in zip(moving, fixed):
        model.append(geom.Vec3(*map(float, a)))
        reference.append(geom.Vec3(*map(float, b)))
    result = SuperposeSVD(model, reference)
    transform = result.transformation
    # Adapt OST's column-vector transform to this API's row-vector convention.
    rotation = np.array([[transform[j, i] for j in range(3)] for i in range(3)])
    translation = np.array([transform[i, 3] for i in range(3)])
    return rotation, translation, float(result.rmsd)


@dataclass(frozen=True)
class Assignment:
    sites: np.ndarray
    distances: np.ndarray
    alternative_gap: float


def assign_centres(points, centres):
    """Assign each observed copy once; unoccupied sites stay unoccupied.

    The alternative gap is the total squared-distance difference to the best
    different assignment, found by excluding each chosen edge in turn. Zero
    means ambiguity; it must not be interpreted as a verified site label.
    """
    from scipy.optimize import linear_sum_assignment

    points, centres = _vectors(points), _vectors(centres)
    if points.shape[1] != centres.shape[1]:
        raise ValueError("copy and centre feature dimensions differ")
    if not len(centres) or len(points) > len(centres):
        raise ValueError("need at least as many sites as observed copies")
    if not len(points):
        return Assignment(np.array([], dtype=int), np.array([]), float("inf"))
    cost = np.sum((points[:, None, :] - centres[None, :, :]) ** 2, axis=2)
    rows, columns = linear_sum_assignment(cost)
    best = cost[rows, columns].sum()
    alternative = float("inf")
    for i, j in zip(rows, columns):
        forbidden = cost.copy()
        forbidden[i, j] = np.inf
        try:
            rr, cc = linear_sum_assignment(forbidden)
        except ValueError:
            continue
        alternative = min(alternative, float(forbidden[rr, cc].sum()))
    return Assignment(
        columns, np.sqrt(cost[rows, columns]), max(0.0, alternative - best)
    )


def group_sites(poses, initial_centres, max_iterations=100):
    """Constrained k-means on aligned centroids or corresponding feature vectors.

    poses maps prediction IDs to arrays. Empty/partial predictions are retained;
    this does not invent a missing ligand. Run multiple initialisations and
    protein anchors externally before accepting a grouping. Returns all copies,
    including distant outliers, with their assignment distances and gaps.
    """
    if not poses:
        raise ValueError("no predictions")
    points = {key: _vectors(value) for key, value in poses.items()}
    centres = _vectors(initial_centres).copy()
    previous = None
    for _ in range(max_iterations):
        assignments = {
            key: assign_centres(value, centres) for key, value in points.items()
        }
        signature = tuple(
            (key, tuple(assignments[key].sites)) for key in sorted(points)
        )
        if signature == previous:
            return assignments, centres
        updated = centres.copy()
        for site in range(len(centres)):
            members = [
                points[key][i]
                for key, a in assignments.items()
                for i, label in enumerate(a.sites)
                if label == site
            ]
            if not members:
                raise ValueError("site has no observed copies")
            updated[site] = np.mean(members, axis=0)
        centres, previous = updated, signature
    raise ValueError("site grouping did not converge")


def same_partition(first, second):
    """Compare copy->site maps, allowing only one global renaming of sites."""
    if first.keys() != second.keys():
        return False
    forward, reverse = {}, {}
    for copy in first:
        a, b = first[copy], second[copy]
        if (a in forward and forward[a] != b) or (b in reverse and reverse[b] != a):
            return False
        forward[a], reverse[b] = b, a
    return True


@dataclass(frozen=True)
class MultiStart:
    assignments: dict
    centres: np.ndarray
    objective: float
    alternative_objective_gap: float
    starts: int


def group_sites_multistart(poses, initialisations):
    """Choose the lowest-cost converged solution among supplied starts.

    Unlike rejecting every local minimum, this follows the clustering objective.
    The gap to a different partition is reported; zero indicates equally good
    competing partitions. A positive gap is not proof of biological correctness
    or of global optimality. Nonconvergence remains explicit.
    """
    solutions = []
    for initial in initialisations:
        assignments, centres = group_sites(poses, initial)
        objective = sum(float(np.sum(a.distances**2)) for a in assignments.values())
        labels = {
            (p, i): int(site)
            for p, a in assignments.items()
            for i, site in enumerate(a.sites)
        }
        solutions.append((objective, assignments, centres, labels))
    if not solutions:
        raise ValueError("at least one initialisation required")
    solutions.sort(key=lambda solution: solution[0])
    best, assignments, centres, labels = solutions[0]
    alternative = min(
        (cost for cost, _, _, other in solutions if not same_partition(labels, other)),
        default=float("inf"),
    )
    return MultiStart(
        assignments, centres, best, max(0.0, alternative - best), len(solutions)
    )


# Historical production site correspondence (unchanged clustering policy).
import dataclasses
import json
import subprocess
import tempfile
from pathlib import Path
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
