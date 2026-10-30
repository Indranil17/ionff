"""Leave one solvent out again, with the held-out molecule made neutral before scoring.

The charge neutrality penalty in the fits acts on training molecules only, so a
held-out molecule can come out charged and the electrostatic error then contains a
net charge error. This script repeats the leave-one-solvent-out test of run_fit.py with
the same fit settings and seed, shifts the charges of the held-out molecule to zero net
charge (ionff.evaluate.neutralise) and writes results/loso_neutralised.json. Reads
data/reference/sapt0_scan.jsonl, or the file given as the first argument.
"""

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.dataset import load_jsonl  # noqa: E402
from ionff.evaluate import loso  # noqa: E402
from run_fit import FITTERS, LOSO_KW, OUT, REF, SEED  # noqa: E402

if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    ds = load_jsonl(REF)
    t0 = time.time()
    rows = loso(ds, FITTERS, neutral=True, seed=SEED, **LOSO_KW)
    for r in rows:
        if "skipped" in r:
            print(r["held_out"], "skipped", flush=True)
            continue
        print(
            r["held_out"],
            " ".join(
                f"{label}: net charge {r[label + '_net_charge']:+.3f} e, elst {r[label]['elst']['rmse']:.1f}, total {r[label]['total']['rmse']:.1f}"
                for label in FITTERS
            ),
            flush=True,
        )
    (OUT / "loso_neutralised.json").write_text(json.dumps({"seed": SEED, "loso_kw": LOSO_KW, "loso": rows}, indent=2))
    print(f"done in {time.time() - t0:.0f} s")
