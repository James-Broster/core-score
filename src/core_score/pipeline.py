"""Orchestrate the analysis stages in execution order."""

from pathlib import Path
from core_score import __version__
from core_score.stages.preparation import prepare_groups
from core_score.stages.site_assignment import prepare_site_groups
from core_score.stages.posebusters import run_pb
from core_score.stages.comparison import compare_groups
from core_score.stages.scoring import score_groups
from core_score.stages.evaluation import evaluate_groups, summarize_reference
from core_score.results import write_json, summarize_run


def execute(
    groups,
    output,
    *,
    workers=1,
    timeout=None,
    ost="ost",
    chemlib=None,
    radius_fallback=True,
    skip_posebusters=False,
    pb_results=None,
    evaluate_reference=False,
    thresholds=(2.0,),
    plots=False,
):
    """Run the full analysis, save its datasets, and return completion counts/status."""
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a new or empty output directory")
    if evaluate_reference and any(not group.get("reference") for group in groups):
        raise ValueError("Reference evaluation requires a reference for every ensemble")

    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "run.json", dict(
        package_version=__version__, ost=ost, chemlib=chemlib,
        timeout=timeout, radius_fallback=radius_fallback, workers=workers,
        posebusters_config="dock", skip_posebusters=skip_posebusters,
        filter_pb_invalid=not skip_posebusters, thresholds=thresholds,
        evaluate_reference=evaluate_reference, plots=plots,
    ))

    print("Preparing structures and ligand chemistry", flush=True)
    prepared = prepare_groups(groups, output)

    print("Applying PoseBusters checks and filtering settings", flush=True)
    filtered = run_pb(
        prepared, output, skip=skip_posebusters,
        filter_invalid=not skip_posebusters, supplied=pb_results,
    )
    site_groups = prepare_site_groups(filtered, output, ost=ost, timeout=timeout)

    print(f"Comparing poses with {workers} simultaneous OST worker(s)", flush=True)
    rows = compare_groups(
        site_groups, output, workers=workers, ost=ost, chemlib=chemlib,
        timeout=timeout, radius_fallback=radius_fallback,
        provenance=dict(package_version=__version__),
    )
    scores = score_groups(site_groups, rows, output)

    labels = []
    if evaluate_reference:
        labels = evaluate_groups(
            filtered, output, thresholds, workers, ost=ost, chemlib=chemlib,
            timeout=timeout, radius_fallback=radius_fallback,
        )
        summarize_reference(scores, labels, output)

    summary = summarize_run(scores, rows, site_groups, labels)
    write_json(output / "status.json", summary["status"])
    if plots:
        from core_score.stages.plotting import plot_results
        plot_results(output, output / "figures")
    return summary
