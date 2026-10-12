"""One plot per threshold: every complex as a dot, ordered by its CORE Score."""

from pathlib import Path

#: Tried in turn for the figure's type. Optima is the paper's face; the rest are the
#: nearest humanist sans faces on other machines, then whatever matplotlib has.
FONTS = ("Optima", "Optima nova LT Pro", "Gill Sans", "Gill Sans MT", "Candara",
         "Trebuchet MS", "Nimbus Sans", "DejaVu Sans")
INK, MUTED, LINE, DOT, FILL = "#16232b", "#5c6b73", "#c8ccce", "#2f6f8f", "#dce9f0"


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
        with plt.rc_context(style(matplotlib)):
            figure, axes = plt.subplots(figsize=(8, 5))
            draw(axes, scores, threshold, cutoffs)
            figure.tight_layout()
            name = "core_scores.png" if len(thresholds) == 1 else f"core_scores_{threshold:g}A.png"
            path = Path(output) / name
            figure.savefig(path, dpi=220)
            plt.close(figure)
        written.append(path)
    return written


def style(matplotlib):
    """Plain, light and airy: hairline axes, horizontal rules only, and the first of FONTS
    this machine actually has."""
    available = {font.name for font in matplotlib.font_manager.fontManager.ttflist}
    family = [name for name in FONTS if name in available] or list(FONTS)
    return {
        "font.family": "sans-serif", "font.sans-serif": family, "font.size": 12,
        "text.color": INK, "axes.labelcolor": MUTED, "axes.labelsize": 12.5,
        "axes.edgecolor": LINE, "axes.linewidth": 1.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 11, "ytick.labelsize": 11,
        "xtick.major.size": 4, "ytick.major.size": 4, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    }


def draw(axes, scores, threshold, cutoffs):
    count = len(scores)
    if count >= 50:  # a tint under the curve; with a handful of complexes the dots stand alone
        axes.fill_between(range(1, count + 1), scores, color=FILL, linewidth=0, zorder=2)
    axes.scatter(range(1, count + 1), scores, s=dot_size(count), color=DOT,
                 alpha=0.75 if count > 500 else 1.0, linewidths=0, zorder=4)
    for cutoff in cutoffs:
        axes.axhline(cutoff, color=MUTED, linewidth=0.9, linestyle=(0, (5, 4)), alpha=0.7, zorder=3)
        above = sum(score >= cutoff for score in scores)
        share = f" · {above / count:.0%}" if count else ""
        axes.annotate(f"{above:,} at or above {cutoff:g}{share}", xy=(0.008, cutoff),
                      xycoords=("axes fraction", "data"), xytext=(0, 6), textcoords="offset points",
                      ha="left", va="bottom", fontsize=10.5, color=MUTED, zorder=5)
    axes.set_ylim(0, 1.04)
    axes.set_xlim(0, max(count, 1))
    axes.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axes.grid(axis="y", color=LINE, linewidth=0.7, alpha=0.6, zorder=0)
    axes.set_axisbelow(True)
    axes.set_xlabel("Complexes, ordered by CORE Score", labelpad=10)
    axes.set_ylabel("CORE Score", labelpad=10)
    axes.xaxis.set_major_formatter(lambda value, _: f"{int(value):,}")
    axes.set_title("CORE Score", loc="left", fontsize=17, pad=26, color=INK)
    axes.annotate(f"{count:,} complexes · threshold {threshold:g} Å",
                  xy=(0, 1), xycoords="axes fraction", xytext=(0, 10), textcoords="offset points",
                  fontsize=11.5, color=MUTED)


def dot_size(count):
    """Smaller dots the more complexes there are, so a big benchmark still reads as a curve."""
    return 26 if count <= 200 else (12 if count <= 1000 else 5)
