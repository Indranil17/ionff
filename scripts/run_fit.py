"""Fit the force field both ways and test transfer to held-out solvents.

Reads data/reference/sapt0_scan.jsonl, writes fitted parameters, a JSON summary
and two figures to results/. Settings are the constants below; pass a different
reference file as the first argument if needed.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.dataset import load_jsonl  # noqa: E402
from ionff.evaluate import KEYS, loso, predict, score, score_by_solvent  # noqa: E402
from ionff.fit import fit_terms, fit_total  # noqa: E402

REF = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "reference" / "sapt0_scan.jsonl"
OUT = ROOT / "results"
FULL_KW = dict(popsize=15, maxiter=600)
LOSO_KW = dict(popsize=10, maxiter=300)
SEED = 0
FITTERS = {"term_resolved": fit_terms, "total_only": fit_total}


def make_figures(ds, params, summary):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"term_resolved": "#1b6ca8", "total_only": "#c0562b"}
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))
    for ax, k in zip(axes, ("elst", "exch", "ind", "disp")):
        ref = ds.energy[k]
        for label, p in params.items():
            ax.scatter(ref, predict(p, ds.r, ds.t)[k], s=6, alpha=0.5, color=colors[label], label=label)
        lim = [min(ref.min(), -1), ref.max()]
        ax.plot(lim, lim, color="0.4", lw=0.8)
        ax.set_title(k)
        ax.set_xlabel("SAPT0 (kJ/mol)")
    axes[0].set_ylabel("force field (kJ/mol)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_term_parity.png", dpi=200)
    plt.close(fig)

    rows = [r for r in summary["loso"] if "skipped" not in r]
    if rows:
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        x = np.arange(len(rows))
        for j, label in enumerate(FITTERS):
            ax.bar(x + 0.2 * (2 * j - 1), [r[label]["total"]["rmse"] for r in rows], 0.38, color=colors[label], label=label)
        ax.set_xticks(x, [r["held_out"] for r in rows])
        ax.set_ylabel("total-energy RMSE on held-out solvent (kJ/mol)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "fig_loso.png", dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    ds = load_jsonl(REF)
    print(f"{len(ds)} points, solvents: {ds.solvents()}", flush=True)
    params, summary = {}, {"n_points": len(ds), "full_kw": FULL_KW, "loso_kw": LOSO_KW, "seed": SEED}
    for label, fit in FITTERS.items():
        t0 = time.time()
        params[label] = fit(ds, seed=SEED, **FULL_KW)
        (OUT / f"params_{label}.json").write_text(json.dumps(params[label].to_dict(), indent=1))
        summary[label] = {"all_data": score(params[label], ds), "by_solvent": score_by_solvent(params[label], ds)}
        print(f"{label}: fitted in {time.time() - t0:.0f} s, total RMSE {summary[label]['all_data']['total']['rmse']:.2f}", flush=True)
    t0 = time.time()
    summary["loso"] = loso(ds, FITTERS, seed=SEED, **LOSO_KW)
    print(f"leave-one-solvent-out done in {time.time() - t0:.0f} s", flush=True)
    for r in summary["loso"]:
        if "skipped" in r:
            print(f"  {r['held_out']}: skipped ({r['skipped']})")
        else:
            print(f"  {r['held_out']}: " + ", ".join(f"{k} {r[k]['total']['rmse']:.2f}" for k in FITTERS))
    (OUT / "fit_summary.json").write_text(json.dumps(summary, indent=1))
    make_figures(ds, params, summary)
