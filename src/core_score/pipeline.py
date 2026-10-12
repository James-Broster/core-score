"""Run the analysis stages in order and save every dataset."""

import os
from importlib.metadata import version
from pathlib import Path
from core_score import __version__
from core_score.inputs import load_inputs
from core_score.results import write_json
from core_score.stages.comparison import compare_groups
from core_score.stages.evaluation import evaluate_groups
from core_score.stages.posebusters import run_posebusters
from core_score.stages.preparation import prepare_groups
from core_score.stages.report import write_report
from core_score.stages.scoring import score_groups
from core_score.stages.site_assignment import prepare_site_groups


def run(inputs, out, *, smiles=None, ligand_chains=None, reference=None, evaluate_reference=False, **options):
    """Run CORE Score from Python, exactly as the core-score command does.

    `inputs` is a JSON manifest or a list of structure files. The other settings match
    the command-line options, for example workers=4, thresholds=(1.0, 2.0) or cutoffs=(0.7,).
    Results are saved under `out`, and the completion counts are returned.

        import core_score
        core_score.run("inputs.json", "results/", workers=8)
    """
    inputs = [inputs] if isinstance(inputs, (str, os.PathLike)) else inputs
    groups = load_inputs(inputs, smiles=smiles, ligand_chains=ligand_chains,
                         reference=reference, evaluate_reference=evaluate_reference)
    return execute(groups, out, evaluate_reference=evaluate_reference, **options)


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
    cutoffs=(0.5, 0.6, 0.7, 0.8, 0.9),
):
    """Run the full analysis on loaded ensembles and return completion counts.

    The files a user reads go in `output`. Everything behind them goes in `output/details`.
    """
    output = Path(output).resolve()
    thresholds, cutoffs = tuple(map(float, thresholds)), tuple(map(float, cutoffs))
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a new or empty output directory")
    details = output / "details"
    details.mkdir(parents=True, exist_ok=True)
    write_json(details / "run.json", dict(
        package_version=__version__, posebusters_version=version("posebusters"), timeout=timeout,
        radius_fallback=radius_fallback, workers=workers, skip_posebusters=skip_posebusters,
        thresholds=thresholds, core_cutoffs=cutoffs, evaluate_reference=evaluate_reference,
    ))

    print("Preparing structures and ligand chemistry", flush=True)
    prepared = prepare_groups(groups, details, timeout=timeout)
    pb_records = None
    if not skip_posebusters:
        print("Running PoseBusters dock checks", flush=True)
        pb_records = run_posebusters(prepared, details)
    sites = prepare_site_groups(prepared, details, timeout=timeout)
    print(f"Comparing poses with {workers} simultaneous OST worker(s)", flush=True)
    rows = compare_groups(sites, details, workers=workers, timeout=timeout, radius_fallback=radius_fallback)
    scores = score_groups(sites, rows, thresholds)
    labels = []
    if evaluate_reference:
        print("Comparing poses with their references", flush=True)
        labels = evaluate_groups(prepared, details, thresholds, workers=workers, timeout=timeout,
                                 radius_fallback=radius_fallback)
    write_report(prepared, scores, pb_records, labels, output, cutoffs)

    summary = summarize_run(scores, rows, sites, labels)
    write_json(details / "status.json", summary["status"])
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
            "unknown_reference_labels": sum(label["rmsd"] is None for label in labels),
        },
    }
