"""Run the analysis stages in order and save every dataset."""

from pathlib import Path
from core_score import __version__
from core_score.results import write_json
from core_score.stages.comparison import compare_groups
from core_score.stages.evaluation import evaluate_groups, summarize_reference
from core_score.stages.posebusters import run_pb
from core_score.stages.preparation import prepare_groups
from core_score.stages.scoring import score_groups
from core_score.stages.site_assignment import prepare_site_groups


def execute(
    groups,
    output,
    *,
    workers=1,
    timeout=None,
    radius_fallback=True,
    skip_posebusters=False,
    evaluate_reference=False,
    thresholds=(2.0,),
    plots=False,
):
    """Run the full analysis on loaded ensembles, save its datasets, and return completion counts."""
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a new or empty output directory")
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "run.json", dict(
        package_version=__version__, timeout=timeout, radius_fallback=radius_fallback,
        workers=workers, skip_posebusters=skip_posebusters, thresholds=thresholds,
        evaluate_reference=evaluate_reference, plots=plots,
    ))

    print("Preparing structures and ligand chemistry", flush=True)
    prepared = prepare_groups(groups, output, timeout=timeout)
    print("Applying PoseBusters checks and filtering settings", flush=True)
    filtered = run_pb(prepared, output, skip=skip_posebusters)
    site_groups = prepare_site_groups(filtered, output, timeout=timeout)
    print(f"Comparing poses with {workers} simultaneous OST worker(s)", flush=True)
    rows = compare_groups(site_groups, output, workers=workers, timeout=timeout, radius_fallback=radius_fallback)
    scores = score_groups(site_groups, rows, output)
    labels = []
    if evaluate_reference:
        labels = evaluate_groups(filtered, output, thresholds, workers=workers, timeout=timeout,
                                 radius_fallback=radius_fallback)
        summarize_reference(scores, labels, output)

    summary = summarize_run(scores, rows, site_groups, labels)
    write_json(output / "status.json", summary["status"])
    if plots:
        from core_score.stages.plotting import plot_results
        plot_results(output, output / "figures")
    return summary


def summarize_run(scores, rows, groups, labels):
    """Completion counts; the measurements themselves are in the saved datasets."""
    failed_comparisons = sum(bool(row["error"]) for row in rows)
    return {
        "comparison_count": len(rows),
        "completed_comparisons": len(rows) - failed_comparisons,
        "score_count": sum(score["status"] == "complete" for score in scores),
        "status": {
            "incomplete_scores": sum(score["status"] != "complete" for score in scores),
            "failed_comparisons": failed_comparisons,
            "unresolved_ensembles": sum(bool(group.get("error")) for group in groups),
            "unknown_reference_labels": sum(label["success"] is None for label in labels),
        },
    }
