"""One figure per threshold: how many complexes have each CORE Score."""

from pathlib import Path


def plot_core_scores(rows, output, cutoffs=(), bin_width=0.05):
    """Write a histogram of the CORE Scores, named for the threshold when there is more
    than one. Dashed lines mark the cutoffs, and the counts above them are in the title."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    written = []
    for threshold in dict.fromkeys(row["threshold_angstrom"] for row in rows):
        scores = [row["core_score"] for row in rows
                  if row["threshold_angstrom"] == threshold and row["core_score"] is not None]
        figure, axes = plt.subplots(figsize=(7.2, 4.6))
        axes.hist(scores, bins=np.arange(0, 1 + bin_width / 2, bin_width),
                  color="#5b8ea6", edgecolor="white", linewidth=0.8, zorder=3)
        for cutoff in cutoffs:
            axes.axvline(cutoff, color="#9a9a9a", linewidth=0.9, linestyle=(0, (4, 3)), zorder=4)
        axes.grid(axis="y", color="#e8e8e8", linewidth=0.7, zorder=0)
        axes.set_axisbelow(True)
        for side in ("top", "right"):
            axes.spines[side].set_visible(False)
        axes.set_xlabel("CORE Score")
        axes.set_ylabel("Complexes")
        axes.set_xlim(0, 1)
        axes.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        axes.set_title(title(scores, threshold, cutoffs), loc="left", fontsize=11.5, linespacing=1.6)
        figure.tight_layout()
        name = "core_scores.png" if len(set(r["threshold_angstrom"] for r in rows)) == 1 \
            else f"core_scores_{threshold:g}A.png"
        path = Path(output) / name
        figure.savefig(path, dpi=200)
        plt.close(figure)
        written.append(path)
    return written


def title(scores, threshold, cutoffs):
    """The run in one line, and how many complexes are at or above each cutoff in the next."""
    heading = f"n = {len(scores)} complexes · threshold {threshold:g} Å"
    counts = "  ·  ".join(f"{sum(score >= cutoff for score in scores)} at or above {cutoff:g}"
                              for cutoff in cutoffs)
    return f"{heading}\n{counts}" if counts else heading
