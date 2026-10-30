"""Draw the two README figures that are not made by another script.

Usage: python scripts/make_readme_figures.py
Reads results/fit_summary.json, results/loso_neutralised.json, results/cluster_summary.json,
data/reference/sapt0_scan.jsonl and data/reference/sapt0_clusters.jsonl.
Writes results/fig_loso_neutralised.png (held-out total RMSE as fitted and after neutralising the molecule,
against the SD of the SAPT0 total) and results/fig_cluster_components.png (force field minus SAPT0 per energy
term against cluster size).
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.dataset import load_jsonl  # noqa: E402

RES = ROOT / "results"
SCAN = ROOT / "data" / "reference" / "sapt0_scan.jsonl"
CLUSTERS = ROOT / "data" / "reference" / "sapt0_clusters.jsonl"
LABELS = ("term_resolved", "total_only")
COLORS = {"term_resolved": "#1b6ca8", "total_only": "#c0562b"}
MARKERS = {"term_resolved": "o", "total_only": "s"}
NAMES = {"term_resolved": "term-resolved fit", "total_only": "total-only fit"}
SOLVENT_NAMES = {"acetone": "acetone", "dmc": "dimethyl carbonate", "dme": "dimethyl ether", "ec": "ethylene carbonate"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
DPI = 200

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 8,
})


def loso_figure():
    """Dumbbell per held-out solvent: open marker as fitted, filled marker after the neutralising shift."""
    raw = {r["held_out"]: r for r in json.loads((RES / "fit_summary.json").read_text())["loso"] if "skipped" not in r}
    neutral = [r for r in json.loads((RES / "loso_neutralised.json").read_text())["loso"] if "skipped" not in r]
    ds = load_jsonl(SCAN)
    fig, ax = plt.subplots(figsize=(6.6, 3.5))
    offsets = {"term_resolved": 0.17, "total_only": -0.17}
    for i, row in enumerate(neutral):
        name = row["held_out"]
        y0 = len(neutral) - 1 - i
        sd = float(np.std(ds.energy["total"][ds.solvent == name]))
        ax.plot([sd, sd], [y0 - 0.38, y0 + 0.38], color=INK, lw=1.6, solid_capstyle="butt", zorder=2)
        for label in LABELS:
            y = y0 + offsets[label]
            before, after = raw[name][label]["total"]["rmse"], row[label]["total"]["rmse"]
            ax.plot([after, before], [y, y], color=COLORS[label], lw=1.2, alpha=0.55, zorder=3)
            ax.plot(before, y, MARKERS[label], ms=7, mfc="#fcfcfb", mec=COLORS[label], mew=1.5, zorder=4)
            ax.plot(after, y, MARKERS[label], ms=7, mfc=COLORS[label], mec="#fcfcfb", mew=1.0, zorder=5)
    ax.set_xscale("log")
    ax.set_xlim(8, 500)
    ax.set_xticks([10, 20, 50, 100, 200, 400], ["10", "20", "50", "100", "200", "400"])
    ax.minorticks_off()
    ax.set_yticks(range(len(neutral)), [SOLVENT_NAMES[r["held_out"]] for r in neutral][::-1])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.set_ylim(-0.6, len(neutral) - 0.4)
    ax.set_xlabel("total-energy RMSE on the held-out solvent (kJ/mol, log scale)")
    handles = [
        Line2D([], [], color=COLORS[l], marker=MARKERS[l], ms=7, lw=0, label=NAMES[l]) for l in LABELS
    ] + [
        Line2D([], [], color=INK2, marker="o", ms=7, mfc="#fcfcfb", mew=1.5, lw=0, label="as fitted"),
        Line2D([], [], color=INK2, marker="o", ms=7, lw=0, label="after neutralising the molecule"),
        Line2D([], [], color=INK, lw=1.6, label="SD of SAPT0 total"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, columnspacing=1.4, handletextpad=0.4)
    fig.tight_layout()
    fig.savefig(RES / "fig_loso_neutralised.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def cluster_figure():
    """Force field minus SAPT0 for each energy term against n, three structures per n as faint points."""
    summary = json.loads((RES / "cluster_summary.json").read_text())
    records = [json.loads(line) for line in CLUSTERS.read_text().splitlines() if line.strip()]
    unconverged = sorted({r["n"] for r in records if r.get("converged") is False})
    panels = [
        ("elst", "electrostatics", 1.0), ("elst", "electrostatics per water", None), ("exch", "exchange", 1.0),
        ("ind", "induction", 1.0), ("disp", "dispersion", 1.0), ("total", "total", 1.0),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.6), sharex=True)
    for ax, (term, title, scale) in zip(axes.flat, panels):
        for n in unconverged:
            ax.axvspan(n - 0.4, n + 0.4, color="#f0efec", lw=0, zorder=0)
        ax.axhline(0, color="#898781", lw=0.8, zorder=1)
        for label in LABELS:
            rows = summary[label]
            ns = sorted({r["n"] for r in rows})
            err = {n: np.array([(r[term]["ff"] - r[term]["sapt0"]) / (n if scale is None else 1) for r in rows if r["n"] == n]) for n in ns}
            ax.plot([n for n in ns for _ in err[n]], np.concatenate([err[n] for n in ns]), MARKERS[label], ms=3, color=COLORS[label], alpha=0.35, zorder=2)
            ax.plot(ns, [err[n].mean() for n in ns], MARKERS[label] + "-", ms=5, lw=1.6, color=COLORS[label], mec="#fcfcfb", mew=0.8, zorder=3)
        ax.set_title(title, loc="left", fontsize=9, color=INK)
        ax.set_xticks(range(1, 7))
    for ax in axes[1]:
        ax.set_xlabel("water molecules n")
    for ax in (axes[0, 0], axes[1, 0]):
        ax.set_ylabel("force field minus SAPT0 (kJ/mol)")
    handles = [Line2D([], [], color=COLORS[l], marker=MARKERS[l], ms=5, lw=1.6, label=NAMES[l]) for l in LABELS]
    handles.append(Line2D([], [], color="#f0efec", lw=8, label="not converged (n = " + ", ".join(map(str, unconverged)) + ")"))
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(RES / "fig_cluster_components.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    loso_figure()
    cluster_figure()
    print("wrote", RES / "fig_loso_neutralised.png", RES / "fig_cluster_components.png")
