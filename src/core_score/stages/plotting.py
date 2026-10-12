"""One plot of the CORE Scores: a dot per complex, stacked where scores coincide."""

from pathlib import Path


def plot_core_scores(rows, output, bin_width=0.02):
    """Write one PNG per threshold, named for the threshold when there is more than one."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    thresholds = list(dict.fromkeys(row["threshold_angstrom"] for row in rows))
    written = []
    for threshold in thresholds:
        scores = sorted(row["core_score"] for row in rows
                        if row["threshold_angstrom"] == threshold and row["core_score"] is not None)
        name = "core_scores.png" if len(thresholds) == 1 else f"core_scores_{threshold:g}A.png"
        path = Path(output) / name
        figure, axes = plt.subplots(figsize=(7, 4.5))
        x, y = stack(scores, bin_width)
        axes.plot(x, y, "o", markersize=5)
        axes.set_xlabel("CORE Score")
        axes.set_ylabel("Complexes")
        axes.set_title(f"{len(scores)} complexes, threshold {threshold:g} Å")
        axes.set_xlim(-0.02, 1.02)
        axes.set_ylim(0, max(y, default=1) + 1)
        figure.tight_layout()
        figure.savefig(path, dpi=200)
        plt.close(figure)
        written.append(path)
    return written


def stack(scores, bin_width):
    """Dot positions for a dot plot: scores within one `bin_width` share a column, drawn at
    the middle of that bin and stacked one above the other."""
    x, y, counts = [], [], {}
    for score in scores:
        column = int(score // bin_width)
        counts[column] = counts.get(column, 0) + 1
        x.append((column + 0.5) * bin_width)
        y.append(counts[column])
    return x, y
