"""CORE Score, implemented exactly per the paper's Methods (complete-linkage clustering;
every group with >=2 poses formed during clustering is a candidate; D = min over candidates
of mean-internal-RMSD / sqrt(fraction-of-ensemble); S = 1/(1 + D/r0), r0 = 1 angstrom)."""

from __future__ import annotations

import itertools
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

R0 = 1.0


def _all_candidate_clusters(n: int, Z: np.ndarray) -> list[list[int]]:
    """Every cluster (>=2 leaves) formed during agglomerative merging, from the linkage matrix."""
    members = {i: [i] for i in range(n)}
    clusters = []
    for i, (a, b, _, _) in enumerate(Z):
        merged = members[int(a)] + members[int(b)]
        members[n + i] = merged
        clusters.append(merged)
    return clusters


def validate_pairwise_matrix(rmsd_matrix: np.ndarray) -> None:
    """Reject an incomplete/malformed matrix outright rather than silently tolerating it. A missing
    (NaN) off-diagonal entry must never be treated as agreement (nansum) or as disagreement (skipped
    via NaN comparisons always being False) -- both hide a candidate group that could have won. The
    diagonal is not checked: it is never read by anything in core_score()."""
    n = rmsd_matrix.shape[0]
    off_diag = ~np.eye(n, dtype=bool)
    vals = rmsd_matrix[off_diag]
    if not np.all(np.isfinite(vals)):
        n_bad = int((~np.isfinite(vals)).sum())
        raise ValueError(
            f"rmsd_matrix has {n_bad} non-finite off-diagonal entries (missing/failed pairwise comparisons); "
            "core_score() requires a complete matrix, not a partially-filled one"
        )
    if np.any(vals < 0):
        raise ValueError(
            "rmsd_matrix has negative off-diagonal entries; RMSD cannot be negative"
        )
    # Off-diagonal only, matching "the diagonal is never read": np.allclose defaults to NaN != NaN,
    # so comparing the whole matrix (diagonal included) would wrongly fail here whenever the diagonal
    # is NaN, even though the actual (off-diagonal) values are symmetric.
    if not np.allclose(rmsd_matrix[off_diag], rmsd_matrix.T[off_diag]):
        raise ValueError("rmsd_matrix's off-diagonal entries are not symmetric")


def core_score(rmsd_matrix: np.ndarray) -> dict:
    """rmsd_matrix: (N, N) symmetric, pairwise ligand RMSD between N poses. The diagonal is never
    read (see validate_pairwise_matrix / the masked mean below) -- only off-diagonal, distinct-pose
    comparisons participate, matching the paper's "mean internal pairwise RMSD" definition exactly.
    Returns dict with S (the CORE Score), D, N, and the winning candidate group's members/size.
    """
    rmsd_matrix = np.asarray(rmsd_matrix, dtype=float)
    if rmsd_matrix.ndim != 2 or rmsd_matrix.shape[0] != rmsd_matrix.shape[1]:
        raise ValueError("RMSD matrix must be square")
    n = rmsd_matrix.shape[0]
    if n < 2:
        raise ValueError("CORE Score requires at least 2 poses")
    rmsd_matrix = np.asarray(rmsd_matrix, dtype=float)
    validate_pairwise_matrix(rmsd_matrix)
    Z = linkage(squareform(rmsd_matrix, checks=False), method="complete")
    best = None
    for members in _all_candidate_clusters(n, Z):
        idx = np.array(members)
        sub = rmsd_matrix[np.ix_(idx, idx)]
        k = len(idx)
        mask = ~np.eye(k, dtype=bool)
        mean_d = sub[mask].sum() / (
            k * (k - 1)
        )  # mean of the k*(k-1) off-diagonal (pairwise) entries only
        d_val = mean_d / np.sqrt(k / n)
        if best is None or d_val < best[0]:
            best = (d_val, members)
    D, members = best
    S = 1.0 / (1.0 + D / R0)
    return dict(S=S, D=D, N=n, group_size=len(members), group_members=members)


def build_matrix(pose_ids, rows):
    """Pairwise ligand-RMSD matrix. Raises unless every pair has exactly one usable RMSD;
    a missing or failed comparison is never filled in."""
    index = {pose_id: i for i, pose_id in enumerate(pose_ids)}
    matrix = np.full((len(pose_ids), len(pose_ids)), np.nan)
    np.fill_diagonal(matrix, 0.0)
    for row in rows:
        if row["n_assigned"] == 1 and not row["error"]:
            i, j = index[row["pose_i"]], index[row["pose_j"]]
            matrix[i, j] = matrix[j, i] = row["rmsd"]
    missing = [(a, b) for (i, a), (j, b) in itertools.combinations(enumerate(pose_ids), 2)
               if not np.isfinite(matrix[i, j])]
    if missing:
        raise ValueError(f"No usable ligand RMSD for pairs {missing}; CORE needs every pair")
    return matrix
