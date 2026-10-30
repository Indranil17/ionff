"""Collect the JSON outputs in results/ into one Markdown table file, results/REPORT.md.

Run after run_fit.py, run_cluster_test.py, run_md.py and analyze_md.py. Missing parts are skipped.
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
LABELS = ("term_resolved", "total_only")
TERMS = ("elst", "exch", "ind", "disp", "total")


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


def fit_section(fit):
    out = [f"## Fit to the SAPT0 scan ({fit['n_points']} points)", "",
           "RMSE in kJ/mol over all points. Term-resolved fits each SAPT0 component on its own; "
           "total-only fits the sum.", "",
           "| fit | " + " | ".join(TERMS) + " |", "|---|" + "---|" * len(TERMS)]
    for label in LABELS:
        out.append(f"| {label} | " + " | ".join(f"{fit[label]['all_data'][k]['rmse']:.2f}" for k in TERMS) + " |")
    out += ["", "## Leave one solvent out", "", "Total-energy RMSE on the held-out solvent (kJ/mol).", "",
            "| held out | term_resolved | total_only |", "|---|---|---|"]
    for r in fit["loso"]:
        if "skipped" in r:
            out.append(f"| {r['held_out']} | skipped | skipped |")
        else:
            out.append(f"| {r['held_out']} | " + " | ".join(f"{r[k]['total']['rmse']:.2f}" for k in LABELS) + " |")
    out += ["", "Skipped means the atom types of that solvent do not occur in the remaining training data.", ""]
    return out


def cluster_section(cl):
    out = ["## Na+(H2O)n many-body test", "",
           "Mean of force field minus SAPT0 over relaxed and perturbed clusters (kJ/mol). "
           "The force field is summed pairwise over the waters.", "",
           "| n | " + " | ".join(f"{k} {m}" for k in LABELS for m in ("total", "ind")) + " |",
           "|---|" + "---|" * (2 * len(LABELS))]
    sizes = sorted({r["n"] for r in cl[LABELS[0]]})
    for n in sizes:
        cells = []
        for k in LABELS:
            rows = [r for r in cl[k] if r["n"] == n]
            for m in ("total", "ind"):
                cells.append(f"{np.mean([r[m]['ff'] - r[m]['sapt0'] for r in rows]):+.1f}")
        out.append(f"| {n} | " + " | ".join(cells) + " |")
    out.append("")
    return out


def md_pooled_section(pooled):
    lo, hi = pooled[LABELS[0]]["fit_window_ps"]

    def pm(v, scale=1.0, nd=3):
        return f"{v[0] * scale:.{nd}f} +/- {v[1] * scale:.{nd}f}"

    out = ["## Na+ in SPC/E water, pooled over seeds", "",
           "Mean over seeds with the standard error over seeds (each seed has 4 Na+). "
           f"D is the MSD slope between {lo:g} and {hi:g} ps, without and with the Yeh-Hummer correction. "
           "Experiment: Na-O 0.24 nm, coordination 4 to 8, D = 1.33e-9 m2/s at 298 K.", "",
           "| parameters | seeds | Na-O peak (nm) | first minimum (nm) | coordination | D raw (1e-9 m2/s) | D corrected (1e-9 m2/s) |",
           "|---|---|---|---|---|---|---|"]
    for label in LABELS:
        s = pooled[label]
        out.append(f"| {label} | {s['n_seeds']} | {pm(s['Na_O_peak_nm'])} | {pm(s['first_minimum_nm'])} | "
                   f"{pm(s['coordination_number'], nd=2)} | {pm(s['D_raw_m2_per_s'], 1e9)} | "
                   f"{pm(s['D_yeh_hummer_m2_per_s'], 1e9)} |")
    out.append("")
    return out


def md_section(summaries):
    out = ["## Na+ in SPC/E water, per seed", "",
           "Experiment: Na-O 0.24 nm, coordination 4 to 8, D = 1.33e-9 m2/s at 298 K. "
           "D is shown without and with the Yeh-Hummer finite-size correction.", "",
           "| parameters | seed | Na-O peak (nm) | coordination | D raw (1e-9 m2/s) | D corrected (1e-9 m2/s) | block error (1e-9) |",
           "|---|---|---|---|---|---|---|"]
    for s in summaries:
        out.append(f"| {s['params']} | {s['seed']} | {s['Na_O_peak_nm']:.3f} | {s['coordination_number']:.1f} | "
                   f"{s['D_raw_m2_per_s'] * 1e9:.2f} | {s['D_yeh_hummer_m2_per_s'] * 1e9:.2f} | "
                   f"{s['D_block_error_m2_per_s'] * 1e9:.2f} |")
    out.append("")
    return out


def main():
    lines = ["# Results", ""]
    fit, cl = load(RES / "fit_summary.json"), load(RES / "cluster_summary.json")
    md = [json.loads(p.read_text()) for p in sorted((RES / "md").glob("summary_*.json"))]
    pooled = load(RES / "md_pooled.json")
    if fit:
        lines += fit_section(fit)
    if cl:
        lines += cluster_section(cl)
    if pooled:
        lines += md_pooled_section(pooled)
    if md:
        lines += md_section(md)
    (RES / "REPORT.md").write_text("\n".join(lines))
    print(f"wrote {RES / 'REPORT.md'}")


if __name__ == "__main__":
    main()
