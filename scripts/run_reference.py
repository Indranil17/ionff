"""Generate the SAPT0 reference scan (resumable).

Edit the settings below or pass solvent names as arguments:
    python scripts/run_reference.py water methanol
Set REF_OUT to write to another file, for example one file per solvent in an array job.
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ionff.reference import run_reference  # noqa: E402

OUT = Path(os.environ.get("REF_OUT", ROOT / "data" / "reference" / "sapt0_scan.jsonl"))
N_DIR = 24
BASIS = "aug-cc-pVDZ"
THREADS = int(os.environ.get("NCPUS", 8))
MEMORY = os.environ.get("PSI4_MEMORY", "8 GB")
ORDER = ["water", "methanol", "dme", "acetone", "ec", "dmc"]

if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    run_reference(OUT, solvents=sys.argv[1:] or ORDER, n_dir=N_DIR, basis=BASIS, threads=THREADS, memory=MEMORY)
