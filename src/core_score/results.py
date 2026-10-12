"""Portable, atomic output and diagnostic records."""

import csv
import hashlib
import json
import math
import traceback
from pathlib import Path


def clean(value):
    if type(value).__name__ in ("NAType", "NaTType"):
        return None
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return clean(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: (
                        json.dumps(clean(v))
                        if isinstance(v, (dict, list, tuple))
                        else clean(v)
                    )
                    for k, v in row.items()
                }
            )
    tmp.replace(path)


def fingerprint(value):
    return hashlib.sha256(json.dumps(clean(value), sort_keys=True).encode()).hexdigest()


def file_hash(path):
    with Path(path).open("rb") as handle:
        digest = hashlib.sha256()
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def diagnostic(stage, exc, **context):
    hint = None
    if isinstance(exc, FileNotFoundError):
        hint = "Check the recorded input or executable path and activate the supplied environment."
    elif "ligand" in str(exc).lower() and (
        "ambiguous" in str(exc).lower() or "absent" in str(exc).lower()
    ):
        hint = "Inspect the reported candidates; provide a ligand chain and, where needed, residue identifier."
    elif "chemical identity" in str(exc).lower():
        hint = "Supply the original ligand SDF or an explicit atom-name mapping; do not repair prediction coordinates."
    return dict(
        stage=stage,
        error=f"{type(exc).__name__}: {exc}",
        suggestion=hint,
        traceback=traceback.format_exc(),
        evidence=getattr(exc, "evidence", None),
        **context,
    )


def summarize_run(scores, rows, groups, labels):
    """Small completion report; detailed measurements are saved separately."""
    failed_comparisons = sum(bool(row.get("error")) for row in rows)
    return {
        "comparison_count": len(rows),
        "completed_comparisons": len(rows) - failed_comparisons,
        "score_count": sum(score["status"] == "complete" for score in scores),
        "status": {
            "incomplete_scores": sum(score["status"] != "complete" for score in scores),
            "failed_comparisons": failed_comparisons,
            "unresolved_ensembles": sum(bool(group.get("stage_error")) for group in groups),
            "unknown_reference_labels": sum(label["success"] is None for label in labels),
        },
    }
