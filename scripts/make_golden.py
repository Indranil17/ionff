"""Write reference values from the Python terms for the C++ test (cxx/tests/golden.csv).

Random parameters inside the fitting bounds and distances from 1.8 to 8 Angstrom,
fixed seed, so the file is reproducible.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.model import BOUNDS  # noqa: E402
from ionff.terms import e_disp, e_elst, e_exch, e_ind  # noqa: E402

OUT = ROOT / "cxx" / "tests" / "golden.csv"
N_ROWS = 400
SEED = 0
ORDER = ("q", "beta", "logA", "B", "alpha", "c", "C6", "d")


def main():
    rng = np.random.default_rng(SEED)
    r = rng.uniform(1.8, 8.0, N_ROWS)
    p = {k: rng.uniform(BOUNDS[k][0], BOUNDS[k][1], N_ROWS) for k in ORDER}
    cols = [
        r,
        *[p[k] for k in ORDER],
        e_elst(r, p["q"], p["beta"]),
        e_exch(r, np.exp(p["logA"]), p["B"]),
        e_ind(r, p["alpha"], p["c"]),
        e_disp(r, p["C6"], p["d"]),
    ]
    header = "# r,q,beta,logA,B,alpha,c,C6,d,elst,exch,ind,disp (Angstrom, e, kJ/mol)"
    np.savetxt(OUT, np.column_stack(cols), delimiter=",", fmt="%.17g", header=header, comments="")
    print(f"wrote {N_ROWS} rows to {OUT}")


if __name__ == "__main__":
    main()
