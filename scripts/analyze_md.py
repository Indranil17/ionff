"""Pool the seeds of the bulk Na+ runs: g(r), coordination, MSD and D per parameter set.

Usage: python scripts/analyze_md.py
Reads results/md/<label>_seed<N>.npz, writes results/md_pooled.json and results/fig_md.png.
D is the mean over seeds (4 ions each) with the standard error over seeds, from the MSD slope between FIT_FROM_PS and FIT_TO_PS.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.md import KB, SPCE_VISCOSITY, YEH_HUMMER, first_shell, msd_fft, rdf_and_cn, unwrap  # noqa: E402

MD = ROOT / "results" / "md"
LABELS = ["term_resolved", "total_only"]
COLORS = {"term_resolved": "#1b6ca8", "total_only": "#c0562b"}
TEMPERATURE = 298.15
FIT_FROM_PS, FIT_TO_PS = 10.0, 80.0
D_EXP = 1.33e-9


def slope_d(msd, stride_ps, lo, hi):
    """D in m^2/s from a straight-line fit of the MSD (nm^2) between lo and hi ps."""
    t = np.arange(len(msd)) * stride_ps
    sel = (t >= lo) & (t <= hi)
    return np.polyfit(t[sel], msd[sel], 1)[0] / 6.0 * 1e-6


def sem(x):
    x = np.asarray(x, float)
    return float(x.std(ddof=1) / np.sqrt(len(x)))


def analyze(label):
    files = sorted(MD.glob(f"{label}_seed*.npz"))
    g_all, cn_all, msd_all, d_raw, d_corr, shells, occupancy = [], [], [], [], [], [], []
    for f in files:
        data = np.load(f)
        box, stride = float(data["box"]), float(data["stride_ps"])
        r, g, cn = rdf_and_cn(data["dist"], box)
        r_peak, r_min, cn_shell = first_shell(r, g, cn)
        msd = msd_fft(unwrap(data["na"], box))
        msd = msd[: len(msd) // 2]
        d = slope_d(msd, stride, FIT_FROM_PS, FIT_TO_PS)
        d_raw.append(d)
        d_corr.append(d + YEH_HUMMER * KB * TEMPERATURE / (6 * np.pi * SPCE_VISCOSITY * box * 1e-9))
        g_all.append(g)
        cn_all.append(cn)
        msd_all.append(msd)
        shells.append((r_peak, r_min, cn_shell))
        occupancy.append((data["dist"] < r_min).sum(axis=-1).ravel())
    shells = np.array(shells)
    occ = np.concatenate(occupancy)
    return {
        "r": r, "g": np.array(g_all), "cn": np.array(cn_all), "msd": np.array(msd_all), "stride": stride,
        "occupancy": np.bincount(occ, minlength=10) / len(occ),
        "summary": {
            "n_seeds": len(files),
            "Na_O_peak_nm": [float(shells[:, 0].mean()), sem(shells[:, 0])],
            "first_minimum_nm": [float(shells[:, 1].mean()), sem(shells[:, 1])],
            "coordination_number": [float(shells[:, 2].mean()), sem(shells[:, 2])],
            "D_raw_m2_per_s": [float(np.mean(d_raw)), sem(d_raw)],
            "D_yeh_hummer_m2_per_s": [float(np.mean(d_corr)), sem(d_corr)],
            "fit_window_ps": [FIT_FROM_PS, FIT_TO_PS],
        },
    }


def main():
    res = {label: analyze(label) for label in LABELS}
    (ROOT / "results" / "md_pooled.json").write_text(
        json.dumps({k: v["summary"] for k, v in res.items()}, indent=2))

    fig, axes = plt.subplots(1, 4, figsize=(14, 3.4))
    for label, v in res.items():
        c = COLORS[label]
        for ax, key in ((axes[0], "g"), (axes[1], "cn")):
            m, e = v[key].mean(axis=0), v[key].std(axis=0, ddof=1) / np.sqrt(len(v[key]))
            ax.plot(v["r"], m, color=c, label=label)
            ax.fill_between(v["r"], m - e, m + e, color=c, alpha=0.25, lw=0)
        axes[2].bar(np.arange(len(v["occupancy"])) + (0.2 if label == "total_only" else -0.2), v["occupancy"], 0.4,
                    color=c, label=label)
        t = np.arange(v["msd"].shape[1]) * v["stride"]
        m = v["msd"].mean(axis=0)
        axes[3].loglog(t[1:], m[1:], color=c, label=label)
    axes[0].axvline(0.24, color="0.4", lw=0.8, ls="--", label="experiment")
    axes[0].set(xlim=(0.15, 0.6), xlabel="r (nm)", ylabel="g(Na-O)")
    axes[1].axhspan(4, 8, color="0.8", alpha=0.6, lw=0, label="experiment 4 to 8")
    axes[1].set(xlim=(0.15, 0.6), ylim=(0, 12), xlabel="r (nm)", ylabel="running coordination number")
    axes[2].set(xlim=(2.5, 9.5), xlabel="waters inside first minimum", ylabel="fraction of frames")
    t_ref = np.array([5.0, 50.0])
    axes[3].loglog(t_ref, 6 * D_EXP * 1e6 * t_ref, color="0.4", lw=0.8, ls="--", label="experiment D")
    axes[3].set(xlabel="lag (ps)", ylabel="MSD (nm$^2$)")
    for ax in axes:
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "results" / "fig_md.png", dpi=200)

    for label, v in res.items():
        s = v["summary"]
        print(label, {k: (np.round(x, 3 if "nm" in k or "coord" in k else 12).tolist() if isinstance(x, list) else x)
                      for k, x in s.items()})
        print("  occupancy 3..9:", np.round(v["occupancy"][3:10], 3).tolist())


if __name__ == "__main__":
    main()
