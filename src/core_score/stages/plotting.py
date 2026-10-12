"""Optional figures drawn from the saved results; nothing is recalculated."""

import csv
import json
from pathlib import Path
import numpy as np
from core_score.results import fingerprint, write_json

METRICS = [("rmsd", "Ligand RMSD (Å)"), ("bb_rmsd", "Pocket backbone RMSD (Å)"), ("lddt_pli", "lDDT-PLI")]


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def plot_results(results, output):
    try:
        import matplotlib
    except ImportError as exc:
        raise RuntimeError("Install plotting support with: python -m pip install 'core-score[plots]'") from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results, output = Path(results), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pairs = read_csv(results / "pairs.csv")
    scores = json.loads((results / "scores.json").read_text())
    images = []
    for key in dict.fromkeys((row["ensemble"], row["site"]) for row in pairs):
        rows = [row for row in pairs if (row["ensemble"], row["site"]) == key]
        ids = sorted({row["pose_i"] for row in rows} | {row["pose_j"] for row in rows})
        index = {pose_id: i for i, pose_id in enumerate(ids)}
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")
        for ax, (metric, title) in zip(axes, METRICS):
            matrix = np.full((len(ids), len(ids)), np.nan)
            for row in rows:
                if row.get(metric):
                    i, j = index[row["pose_i"]], index[row["pose_j"]]
                    matrix[i, j] = matrix[j, i] = float(row[metric])
            cmap = plt.get_cmap("viridis").copy()
            cmap.set_bad("#dddddd")
            image = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1 if metric == "lddt_pli" else None)
            fig.colorbar(image, ax=ax, shrink=0.75)
            ax.set_xticks(range(len(ids)), ids, rotation=90, fontsize=7)
            ax.set_yticks(range(len(ids)), ids, fontsize=7)
            ax.set_title(title)
        for row in rows:
            if row.get("radius_fallback_used") == "True":
                axes[0].text(index[row["pose_i"]], index[row["pose_j"]], "R",
                             ha="center", va="center", color="red", fontsize=8)
        core = ", ".join(
            f"{s['core_score']:.3f} (T = {s['threshold_angstrom']:g} Å)"
            for s in scores if (s["ensemble"], s["site"]) == key and s["core_score"] is not None
        ) or "unavailable"
        fig.suptitle(f"{key[0]} · site {key[1]} · CORE {core}\nR = radius recovery; grey = missing or self comparison")
        filename = fingerprint(key)[:16] + ".png"
        fig.savefig(output / filename, dpi=160)
        plt.close(fig)
        images.append(dict(ensemble=key[0], site=key[1], path=filename))
    write_json(output / "index.json", images)

    summary = results / "reference_summary.csv"
    if summary.exists():
        rows = [row for row in read_csv(summary) if row["correct_fraction"]]
        fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
        for threshold in sorted({row["threshold_angstrom"] for row in rows}, key=float):
            chosen = [row for row in rows if row["threshold_angstrom"] == threshold]
            ax.scatter([float(row["core_score"]) for row in chosen], [float(row["correct_fraction"]) for row in chosen],
                       label=f"T = {float(threshold):g} Å", alpha=0.7)
        ax.set(xlabel="CORE Score", ylabel="Fraction of poses within T", xlim=(0, 1), ylim=(-0.03, 1.03),
               title="Reference evaluation (descriptive, not calibrated)")
        if rows:
            ax.legend()
        fig.savefig(output / "reference_evaluation.png", dpi=160)
        plt.close(fig)
    return images
