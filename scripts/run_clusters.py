"""Relax Na+(H2O)n for n = 1 to 6 and compute SAPT0 of the cluster against Na+ (resumable).

Pass cluster sizes as arguments to run a subset, for example: python scripts/run_clusters.py 4 5 6
Set CLUSTER_OUT to write to another file, for example one file per size in an array job.
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.clusters import run_clusters  # noqa: E402

OUT = Path(os.environ.get("CLUSTER_OUT", ROOT / "data" / "reference" / "sapt0_clusters.jsonl"))
SIZES = tuple(int(a) for a in sys.argv[1:]) or (1, 2, 3, 4, 5, 6)
THREADS = int(os.environ.get("NCPUS", 8))
MEMORY = os.environ.get("PSI4_MEMORY", "8 GB")

if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    run_clusters(OUT, sizes=SIZES, threads=THREADS, memory=MEMORY)
