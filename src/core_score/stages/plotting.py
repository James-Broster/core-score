"""Optional static plots; never recompute scientific measurements."""

import csv
import json
from pathlib import Path
import numpy as np


def plot_results(results, output):
    try:
        import matplotlib
    except ImportError as exc:
        raise RuntimeError(
            "Install plotting support with: python -m pip install 'core-score[plots]'"
        ) from exc

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from core_score.results import fingerprint, write_json

    results, output = Path(results), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((results / "analysis_manifest.json").read_text())
    with (results / "pairs.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    scores = (
        json.loads((results / "scores.json").read_text())
        if (results / "scores.json").exists()
        else []
    )
    images = []
    for group in manifest["groups"]:
        ids = sorted(p["pose_id"] for p in group["poses"])
        if not ids:
            continue
        index = {p: i for i, p in enumerate(ids)}
        score = next(
            (
                s
                for s in scores
                if (s["ensemble"], s["site"]) == (group["ensemble"], group["site"])
            ),
            {},
        )
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")
        for ax, metric, title in zip(
            axes,
            ["rmsd", "bb_rmsd_mean", "lddt_pli"],
            ["Ligand RMSD (Å)", "Pocket backbone RMSD (Å)", "lDDT-PLI"],
        ):
            matrix = np.full((len(ids), len(ids)), np.nan)
            for r in rows:
                if r["ensemble"] != group["ensemble"] or str(r["site"]) != str(
                    group["site"]
                ):
                    continue
                try:
                    value = float(r[metric])
                except (ValueError, KeyError, TypeError):
                    continue
                if r["pose_i"] in index and r["pose_j"] in index:
                    i, j = index[r["pose_i"]], index[r["pose_j"]]
                    matrix[i, j] = matrix[j, i] = value
            cmap = plt.get_cmap("viridis").copy()
            cmap.set_bad("#dddddd")
            image = ax.imshow(
                matrix, cmap=cmap, vmin=0, vmax=1 if metric == "lddt_pli" else None
            )
            fig.colorbar(image, ax=ax, shrink=0.75)
            labels = [
                p + (" *" if p in score.get("group_pose_ids", []) else "") for p in ids
            ]
            ax.set_xticks(range(len(ids)), labels, rotation=90, fontsize=7)
            ax.set_yticks(range(len(ids)), labels, fontsize=7)
            ax.set_title(title)
            if metric == "rmsd":
                for r in rows:
                    if (
                        r["ensemble"] == group["ensemble"]
                        and str(r["site"]) == str(group["site"])
                        and r.get("radius_fallback_used", "").lower() == "true"
                    ):
                        if r["pose_i"] in index and r["pose_j"] in index:
                            ax.text(
                                index[r["pose_i"]],
                                index[r["pose_j"]],
                                "R",
                                ha="center",
                                va="center",
                                color="red",
                                fontsize=8,
                            )
        label = f'{score["S"]:.3f}' if score.get("S") is not None else "unavailable"
        fig.suptitle(
            f'{group["ensemble"]} · site {group["site"]} · CORE {label}\n* winning group; R = radius recovery; grey = missing/self comparison'
        )
        filename = fingerprint((group["ensemble"], group["site"]))[:16] + ".png"
        fig.savefig(output / filename, dpi=160)
        plt.close(fig)
        images.append(
            dict(ensemble=group["ensemble"], site=group["site"], path=filename)
        )
    write_json(output / "index.json", images)
    summary = results / "reference_summary.csv"
    if summary.exists():
        with summary.open(newline="") as handle:
            evaluations = list(csv.DictReader(handle))
        fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
        for threshold in sorted(
            {r["threshold_angstrom"] for r in evaluations}, key=float
        ):
            values = [
                r
                for r in evaluations
                if r["threshold_angstrom"] == threshold
                and r["winning_group_success_fraction"] != ""
            ]
            if values:
                ax.scatter(
                    [float(r["core_score"]) for r in values],
                    [float(r["winning_group_success_fraction"]) for r in values],
                    label=f"< {threshold} Å",
                    alpha=0.7,
                )
        ax.set(
            xlabel="CORE score",
            ylabel="Winning group fraction successful",
            xlim=(0, 1),
            ylim=(-0.03, 1.03),
            title="Reference evaluation (descriptive; not calibrated)",
        )
        if ax.collections:
            ax.legend()
        fig.savefig(output / "reference_evaluation.png", dpi=160)
        plt.close(fig)
    return images
