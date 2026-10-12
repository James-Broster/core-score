"""Build complete matrices; never impute missing or failed comparisons."""

import itertools
import math
import numpy as np


def build_matrix(pose_ids, rows, metric="rmsd"):
    if len(pose_ids) != len(set(pose_ids)):
        raise ValueError("Duplicate pose identifiers")
    index = {p: i for i, p in enumerate(pose_ids)}
    matrix = np.full((len(index), len(index)), np.nan)
    np.fill_diagonal(matrix, 0.0)
    seen, invalid = set(), []
    for row in rows:
        a, b = row["pose_i"], row["pose_j"]
        if a not in index or b not in index or a == b:
            raise ValueError(f"Unexpected comparison: {a}, {b}")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError(f"Duplicate comparison: {key}")
        seen.add(key)
        try:
            value = float(row[metric])
        except (KeyError, TypeError, ValueError):
            value = float("nan")
        if (
            row.get("error")
            or not math.isfinite(value)
            or value < 0
            or int(row.get("n_assigned", 0)) != 1
        ):
            invalid.append(key)
            continue
        matrix[index[a], index[b]] = matrix[index[b], index[a]] = value
    missing = sorted(
        set(tuple(sorted(p)) for p in itertools.combinations(pose_ids, 2)) - seen
    )
    if missing or invalid:
        raise ValueError(
            f"Incomplete pairwise matrix; missing={missing}; failed={invalid}"
        )
    return matrix
