"""Score the fitted force fields against the Na+(H2O)n SAPT0 references and plot the result."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.evaluate import cluster_scores  # noqa: E402
from ionff.model import Params  # noqa: E402

REF = ROOT / "data" / "reference" / "sapt0_clusters.jsonl"
OUT = ROOT / "results"
LABELS = ("term_resolved", "total_only")

if __name__ == "__main__":
    records = [json.loads(line) for line in REF.read_text().splitlines() if line.strip()]
    summary = {}
    for label in LABELS:
        params = Params.from_dict(json.loads((OUT / f"params_{label}.json").read_text()))
        rows = cluster_scores(params, records)
        summary[label] = [
            {"n": r["n"], **{k: {"ff": r[k]["ff"], "sapt0": r[k]["sapt0"]} for k in ("elst", "exch", "ind", "disp", "total")}}
            for r in rows
        ]
    (OUT / "cluster_summary.json").write_text(json.dumps(summary, indent=1))
    for label, rows in summary.items():
        print(label)
        for n in sorted({r["n"] for r in rows}):
            sel = [r for r in rows if r["n"] == n]
            err = np.mean([abs(r["total"]["ff"] - r["total"]["sapt0"]) for r in sel])
            ind = np.mean([r["ind"]["ff"] - r["ind"]["sapt0"] for r in sel])
            print(f"  n={n}: |total error| {err:6.1f} kJ/mol, mean induction error {ind:6.1f} kJ/mol")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharex=True)
    colors = {"term_resolved": "#1b6ca8", "total_only": "#c0562b"}
    for label, rows in summary.items():
        ns = sorted({r["n"] for r in rows})
        for ax, key in zip(axes, ("total", "ind")):
            err = [np.mean([r[key]["ff"] - r[key]["sapt0"] for r in rows if r["n"] == n]) for n in ns]
            ax.plot(ns, err, "o-", color=colors[label], label=label)
    for ax, key in zip(axes, ("total", "induction")):
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlabel("n in Na$^+$(H$_2$O)$_n$")
        ax.set_title(f"{key}: force field minus SAPT0")
    axes[0].set_ylabel("kJ/mol")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_clusters.png", dpi=200)
