"""Run the analysis stages in order and save every dataset."""

from importlib.metadata import version
from pathlib import Path
from core_score import __version__
from core_score.results import write_json
from core_score.stages.comparison import compare_groups
from core_score.stages.evaluation import evaluate_groups
from core_score.stages.posebusters import run_posebusters
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
        package_version=__version__, posebusters_version=version("posebusters"), timeout=timeout,
        radius_fallback=radius_fallback, workers=workers, skip_posebusters=skip_posebusters,
        thresholds=thresholds, evaluate_reference=evaluate_reference, plots=plots,
    ))

    print("Preparing structures and ligand chemistry", flush=True)
    prepared = prepare_groups(groups, output, timeout=timeout)
    if not skip_posebusters:
        print("Running PoseBusters dock checks", flush=True)
        run_posebusters(prepared, output)
    sites = prepare_site_groups(prepared, output, timeout=timeout)
    print(f"Comparing poses with {workers} simultaneous OST worker(s)", flush=True)
    rows = compare_groups(sites, output, workers=workers, timeout=timeout, radius_fallback=radius_fallback)
    scores = score_groups(sites, rows, output, thresholds)
    labels = []
    if evaluate_reference:
        print("Comparing poses with their references", flush=True)
        labels = evaluate_groups(prepared, scores, output, thresholds, workers=workers, timeout=timeout,
                                 radius_fallback=radius_fallback)

    summary = summarize_run(scores, rows, sites, labels)
    write_json(output / "status.json", summary["status"])
    if plots:
        from core_score.stages.plotting import plot_results
        plot_results(output, output / "figures")
    return summary


def summarize_run(scores, rows, groups, labels):
    """Completion counts; the measurements themselves are in the saved datasets."""
    ensembles = {(score["ensemble"], score["site"]) for score in scores}
    scored = {(score["ensemble"], score["site"]) for score in scores if score["core_score"] is not None}
    failed_comparisons = sum(bool(row["error"]) for row in rows)
    return {
        "comparison_count": len(rows),
        "completed_comparisons": len(rows) - failed_comparisons,
        "score_count": len(scored),
        "status": {
            "incomplete_scores": len(ensembles - scored),
            "failed_comparisons": failed_comparisons,
            "unresolved_ensembles": sum(bool(group.get("error")) for group in groups),
            "unknown_reference_labels": sum(label["success"] is None for label in labels),
        },
    }
