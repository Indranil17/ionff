"""Bulk Na+ in SPC/E water with the fitted pair terms, then g(r), coordination number and diffusion.

Usage: python scripts/run_md.py [params.json] [prod_ps]
Defaults: results/params_term_resolved.json and 1000 ps of NVE after 200 ps of NPT.
Set MD_PLATFORM (CPU, CUDA, OpenCL) to pick an OpenMM platform and MD_SEED to change the seed.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.md import diffusion, first_shell, rdf_and_cn, run  # noqa: E402
from ionff.model import Params  # noqa: E402

PARAMS = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "params_term_resolved.json"
PROD_PS = float(sys.argv[2]) if len(sys.argv) > 2 else 1000.0
EQUIL_PS = float(os.environ.get("MD_EQUIL_PS", 200.0))
BOX_NM = float(os.environ.get("MD_BOX_NM", 3.0))
N_NA = int(os.environ.get("MD_N_NA", 4))
SEED = int(os.environ.get("MD_SEED", 1))
PLATFORM = os.environ.get("MD_PLATFORM") or None
TEMPERATURE = 298.15
OUT = ROOT / "results" / "md"

# experiment: D from limiting conductivity 50.1 S cm2/mol, first shell from X-ray and neutron diffraction
EXPERIMENT = {"D_m2_per_s": 1.33e-9, "Na_O_nm": 0.24, "coordination": "4 to 8 (Ohtaki and Radnai 1993)"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    params = Params.from_dict(json.loads(PARAMS.read_text()))
    tag = f"{PARAMS.stem.replace('params_', '')}_seed{SEED}"
    npz = run(params, OUT / tag, box_nm=BOX_NM, n_na=N_NA, temperature=TEMPERATURE, equil_ps=EQUIL_PS,
              prod_ps=PROD_PS, platform=PLATFORM, seed=SEED)
    data = np.load(npz)
    box = float(data["box"])
    r, g, cn = rdf_and_cn(data["dist"], box)
    r_peak, r_min, cn_shell = first_shell(r, g, cn)
    d_raw, d_corr, d_err = diffusion(data["na"], box, float(data["stride_ps"]), temperature=TEMPERATURE)
    summary = {
        "params": PARAMS.name, "seed": SEED, "box_nm": box, "n_na": N_NA, "equil_ps": EQUIL_PS, "prod_ps": PROD_PS,
        "Na_O_peak_nm": float(r_peak), "first_minimum_nm": float(r_min), "coordination_number": float(cn_shell),
        "D_raw_m2_per_s": d_raw, "D_yeh_hummer_m2_per_s": d_corr, "D_block_error_m2_per_s": d_err,
        "experiment": EXPERIMENT,
    }
    (OUT / f"summary_{tag}.json").write_text(json.dumps(summary, indent=2))
    np.savetxt(OUT / f"rdf_{tag}.dat", np.c_[r, g, cn], header="r_nm g_NaO CN")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
