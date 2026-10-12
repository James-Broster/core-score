"""Run one OST process per pair and retain OST's own scores and assignments."""

import json
import itertools
from functools import partial
from concurrent.futures import ProcessPoolExecutor
import os
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from core_score.stages.utils.radius import score_with_radius_fallback
from core_score.results import diagnostic, write_json, write_csv, file_hash


def run_ost(model, reference, ost="ost", timeout=None, chemlib=None, radius=4.0):
    env = dict(
        os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    if chemlib:
        env["OST_COMPOUNDS_CHEMLIB"] = str(Path(chemlib).resolve())
    with tempfile.TemporaryDirectory() as directory:
        output = str(Path(directory) / "ost.json")
        command = [ost, "compare-ligand-structures", "-m", model["path"]]
        if model.get("ligand_sdfs"):
            command += ["-ml", *model["ligand_sdfs"]]
        command += ["-r", reference["path"]]
        if reference.get("ligand_sdfs"):
            command += ["-rl", *reference["ligand_sdfs"]]
        command += [
            "--fault-tolerant",
            "--rmsd",
            "--full-bs-search",
            "--full-results",
            "--lddt-pli",
            "--lddt-pli-add-mdl-contacts",
            "-o",
            output,
        ]
        if radius != 4.0:
            command += ["--radius", str(radius)]
        evidence = dict(
            command=command,
            timeout=timeout,
            radius=radius,
            chemlib=env.get("OST_COMPOUNDS_CHEMLIB"),
            model=model,
            reference=reference,
        )
        try:
            process = subprocess.run(
                command, env=env, capture_output=True, text=True, timeout=timeout
            )
            evidence.update(
                returncode=process.returncode,
                stdout=process.stdout,
                stderr=process.stderr,
            )
            evidence["raw_output"] = (
                Path(output).read_text() if Path(output).exists() else None
            )
            result = json.loads(Path(output).read_text())
            if process.returncode and result.get("status") != "FAILURE":
                result = dict(
                    result, status="FAILURE", exception=process.stderr[-2000:]
                )
            if result.get("status") not in ("SUCCESS", "FAILURE"):
                result.update(
                    status="FAILURE", exception="OST did not report a recognised status"
                )
            result["_execution"] = evidence
            return result
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            if isinstance(exc, subprocess.TimeoutExpired):
                evidence.update(
                    stdout=(
                        exc.stdout.decode(errors="replace")
                        if isinstance(exc.stdout, bytes)
                        else exc.stdout
                    ),
                    stderr=(
                        exc.stderr.decode(errors="replace")
                        if isinstance(exc.stderr, bytes)
                        else exc.stderr
                    ),
                )
            report = diagnostic("comparison", exc, **evidence)
            if isinstance(exc, subprocess.TimeoutExpired):
                report["suggestion"] = (
                    "Retry with a larger --timeout or omit --timeout for no limit."
                )
            return dict(
                status="FAILURE",
                exception=f"{type(exc).__name__}: {exc}",
                _execution=report,
            )


def summarize(group, pose_i, pose_j, model, reference, result):
    copies = group.get("n_input_copies")
    row = dict(
        method=group["method"],
        target=group["target"],
        ligand_key=group["ligand_key"],
        pose_i_seed=pose_i["seed"],
        pose_i_sample=pose_i["sample"],
        pose_j_seed=pose_j["seed"],
        pose_j_sample=pose_j["sample"],
        n_copies_i=model["copies"],
        n_copies_j=reference["copies"],
        n_input_copies=copies,
        multi_copy_group=bool(copies is not None and copies >= 2),
        copies_consistent=bool(
            copies is not None and model["copies"] == copies == reference["copies"]
        ),
    )
    failed = result.get("status") == "FAILURE"
    assigned = (
        []
        if failed
        else [
            r
            for r in result.get("rmsd", {}).get("assigned_scores", [])
            if r.get("score") is not None
        ]
    )
    contacts = (
        []
        if failed
        else [
            r["score"]
            for r in result.get("lddt_pli", {}).get("assigned_scores", [])
            if r.get("score") is not None
        ]
    )
    values = {
        "rmsd": [r["score"] for r in assigned],
        "lddt_pli": contacts,
        "bb_rmsd": [r["bb_rmsd"] for r in assigned if r.get("bb_rmsd") is not None],
        "lddt_lp": [r["lddt_lp"] for r in assigned if r.get("lddt_lp") is not None],
    }
    row["n_assigned"] = len(assigned)
    for metric, scores in values.items():
        if metric in ("rmsd", "lddt_pli"):
            row[metric] = scores[0] if len(scores) == 1 else float("nan")
        row[metric + "_mean"] = float(np.mean(scores)) if scores else float("nan")
        row[metric + "_min"] = min(scores) if scores else float("nan")
        row[metric + "_max"] = max(scores) if scores else float("nan")
    row["error"] = (
        str(result.get("exception") or result.get("traceback") or "OST failed")
        if failed
        else ""
    )
    return row


def compare_pair(task, *, ost="ost", chemlib=None, timeout=None, radius_fallback=True):
    group, pose_i, pose_j, model, reference = task
    fallback = (6.0, 8.0) if radius_fallback else ()
    if group.get("site") == "reference" and (model["copies"] != 1 or reference["copies"] != 1):
        fallback = ()
    if fallback and (model["copies"] != 1 or reference["copies"] != 1):
        raise ValueError(
            "Radius fallback requires one explicit ligand copy per side; split sites first"
        )
    attempts = []

    def evaluate(radius):
        result = run_ost(model, reference, ost=ost, chemlib=chemlib, timeout=timeout, radius=radius)
        attempts.append(result)
        assigned = (
            result.get("rmsd", {}).get("assigned_scores", [])
            if result.get("status") != "FAILURE"
            else []
        )
        return dict(
            score=assigned[0].get("score") if len(assigned) == 1 else None,
            reason=(
                result.get("rmsd", {})
                if result.get("status") != "FAILURE"
                else result.get("exception", "OST failure")
            ),
        )

    _, metadata = score_with_radius_fallback(evaluate, fallback)
    result = dict(attempts[0])
    if metadata["radius_fallback_used"]:
        # Only binding-site-derived RMSD metrics use the enlarged site.
        # Keep reference-independent lDDT-PLI from the original call.
        result["rmsd"] = attempts[-1]["rmsd"]
    row = summarize(group, pose_i, pose_j, model, reference, result)
    if row["n_assigned"] == 0 and not row["error"]:
        row["error"] = (
            "No usable ligand assignment; see OST reasons in ligand_details.jsonl"
        )
    row.update(
        ensemble=group.get("ensemble", group["target"]),
        site=group.get("site", "1"),
        pose_i=pose_i.get("pose_id", str(pose_i["sample"])),
        pose_j=pose_j.get("pose_id", str(pose_j["sample"])),
    )
    row.update(metadata)
    detail = {
        k: row[k]
        for k in (
            "method",
            "target",
            "ligand_key",
            "pose_i_seed",
            "pose_i_sample",
            "pose_j_seed",
            "pose_j_sample",
            "n_copies_i",
            "n_copies_j",
            "error",
        )
    }
    detail.update(
        source_pose_i=pose_i["path"],
        source_pose_j=pose_j["path"],
        ost_results=attempts,
        radius_policy=metadata,
        cif_record_filter_version=1,
    )
    detail.update(
        ensemble=row["ensemble"],
        site=row["site"],
        pose_i=row["pose_i"],
        pose_j=row["pose_j"],
        stage="comparison",
        suggestion=(
            None
            if not row["error"]
            else "Inspect ost_results for native assignment or chemistry reasons; do not substitute zero for a missing distance."
        ),
    )
    return row, detail


def probe_prepared(prepared, *, ost="ost", chemlib=None, timeout=None):
    """Detect missing native ligand support using OST, not a model-name guess."""
    if prepared.get("ligand_sdfs"):
        return prepared
    result = run_ost(
        prepared,
        prepared,
        ost=ost, chemlib=chemlib, timeout=timeout,
    )
    prepared = dict(prepared, detection_probe=result)
    if result.get("status") == "FAILURE":
        raise RuntimeError("OST input probe failed: " + json.dumps(result))
    if not result.get("model_ligands"):
        copies = prepared.get("copy_records", [])
        if not copies or any(not r.get("ost_sdf") for r in copies):
            raise ValueError(
                "OST could not detect the ligand; supply original ligand SDFs"
            )
        prepared.update(
            ligand_sdfs=[r["ost_sdf"] for r in copies],
            route="explicit_connectivity_sdf",
        )
    return prepared


def pair_job(task, *, ost="ost", chemlib=None, timeout=None, radius_fallback=True):
    try:
        return compare_pair(task, ost=ost, chemlib=chemlib, timeout=timeout, radius_fallback=radius_fallback)
    except Exception as exc:
        group, a, b, _, _ = task
        report = diagnostic(
            "comparison",
            exc,
            ensemble=group["ensemble"],
            site=group["site"],
            pose_i=a["pose_id"],
            pose_j=b["pose_id"],
            inputs=[a["path"], b["path"]],
        )
        return (
            dict(
                ensemble=group["ensemble"],
                site=group["site"],
                pose_i=a["pose_id"],
                pose_j=b["pose_id"],
                rmsd=None,
                n_assigned=0,
                error=report["error"],
            ),
            report,
        )


def prepare_pair_tasks(groups, *, ost="ost", chemlib=None, timeout=None):
    tasks, rows, details = [], [], []
    for group in groups:
        if group.get("stage_error"):
            continue
        poses = sorted(
            group["poses"], key=lambda p: (p["seed"], p["sample"], p["pose_id"])
        )
        for pose in poses:
            try:
                pose["prepared"] = probe_prepared(pose["prepared"], ost=ost, chemlib=chemlib, timeout=timeout)
            except Exception as exc:
                pose["probe_error"] = diagnostic(
                    "ligand_detection", exc, structure=pose["path"]
                )
        for a, b in itertools.combinations(poses, 2):
            if a.get("probe_error") or b.get("probe_error"):
                detail = dict(
                    ensemble=group["ensemble"],
                    site=group["site"],
                    pose_i=a["pose_id"],
                    pose_j=b["pose_id"],
                    stage="ligand_detection",
                    errors=[a.get("probe_error"), b.get("probe_error")],
                )
                rows.append(
                    dict(
                        ensemble=group["ensemble"],
                        site=group["site"],
                        pose_i=a["pose_id"],
                        pose_j=b["pose_id"],
                        rmsd=None,
                        n_assigned=0,
                        error="Input probe failed; inspect ligand_details.jsonl",
                    )
                )
                details.append(detail)
                continue
            task = (group, a, b, a["prepared"], b["prepared"])
            tasks.append(task)
    return tasks, rows, details


def write_pair_results(groups, rows, details, output, provenance):
    output = Path(output)
    for group in groups:
        group["input_fingerprints"] = {
            p["pose_id"]: {
                f: file_hash(f)
                for f in (
                    [p["prepared"]["path"]] + p["prepared"].get("ligand_sdfs", [])
                )
            }
            for p in group["poses"]
            if "prepared" in p
        }
    write_csv(output / "pairs.csv", rows)
    with (output / "ligand_details.jsonl").open("w") as handle:
        from core_score.results import clean

        for detail in details:
            handle.write(json.dumps(clean(detail), allow_nan=False) + "\n")
    write_json(
        output / "analysis_manifest.json",
        dict(schema="core-score.analysis.v1", groups=groups, provenance=provenance),
    )


def run_pair_tasks(tasks, workers=1, *, ost="ost", chemlib=None, timeout=None, radius_fallback=True):
    """Yield pair results in task order, serially or using local worker processes."""
    worker = partial(pair_job, ost=ost, chemlib=chemlib, timeout=timeout, radius_fallback=radius_fallback)
    if workers == 1:
        for task in tasks:
            yield worker(task)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            yield from pool.map(worker, tasks)
