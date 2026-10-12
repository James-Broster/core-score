"""One plot per threshold: every complex as a dot, ordered by its CORE Score."""

from pathlib import Path


def plot_core_scores(rows, output, cutoffs=()):
    """Write the CORE Scores of every complex, lowest first, named for the threshold when
    there is more than one. A dashed line marks each cutoff, labelled with how many
    complexes reach it."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    thresholds = list(dict.fromkeys(row["threshold_angstrom"] for row in rows))
    written = []
    for threshold in thresholds:
        scores = sorted(row["core_score"] for row in rows
                        if row["threshold_angstrom"] == threshold and row["core_score"] is not None)
        figure, axes = plt.subplots(figsize=(7.2, 4.6))
        axes.scatter(range(1, len(scores) + 1), scores, s=dot_size(len(scores)),
                     color="#3f6f8f", alpha=0.6 if len(scores) > 500 else 1.0, linewidths=0, zorder=3)
        for cutoff in cutoffs:
            axes.axhline(cutoff, color="#9a9a9a", linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
            above = sum(score >= cutoff for score in scores)
            axes.text(len(scores) * 0.005, cutoff + 0.018, f"{above} at or above {cutoff:g}",
                      fontsize=9.5, color="#4d4d4d")
        axes.grid(axis="y", color="#e8e8e8", linewidth=0.7, zorder=0)
        axes.set_axisbelow(True)
        for side in ("top", "right"):
            axes.spines[side].set_visible(False)
        axes.set_xlabel("Complexes, ordered by CORE Score")
        axes.set_ylabel("CORE Score")
        axes.set_xlim(0, max(len(scores), 1))
        axes.set_ylim(0, 1.02)
        axes.set_title(f"n = {len(scores)} complexes · threshold {threshold:g} Å", loc="left")
        figure.tight_layout()
        name = "core_scores.png" if len(thresholds) == 1 else f"core_scores_{threshold:g}A.png"
        path = Path(output) / name
        figure.savefig(path, dpi=200)
        plt.close(figure)
        written.append(path)
    return written


def dot_size(count):
    """Smaller dots the more complexes there are, so a big benchmark still reads as a curve."""
    return 24 if count <= 200 else (12 if count <= 1000 else 4)
