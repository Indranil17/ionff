"""SAPT0 reference energies for Na+ at fixed positions around a solvent molecule."""

import json
import time
from pathlib import Path

from .geometry import SOLVENTS, load_monomer, sample_na_positions

HARTREE_TO_KJ = 2625.4996394799
COMPONENTS = {
    "elst": "SAPT ELST ENERGY",
    "exch": "SAPT EXCH ENERGY",
    "ind": "SAPT IND ENERGY",
    "disp": "SAPT DISP ENERGY",
    "total": "SAPT0 TOTAL ENERGY",
}


def _geometry_string(symbols, pos, na_xyz, charge_b=1):
    lines = [f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, (x, y, z) in zip(symbols, pos)]
    lines = ["0 1", *lines, "--", f"{charge_b} 1", "Na {:.8f} {:.8f} {:.8f}".format(*na_xyz)]
    lines += ["units angstrom", "no_com", "no_reorient", "symmetry c1"]
    return "\n".join(lines)


def sapt0_point(psi4, symbols, pos, na_xyz, basis="aug-cc-pVDZ"):
    """One SAPT0 calculation of solvent (fragment A) and Na+ (fragment B), kJ/mol."""
    mol = psi4.geometry(_geometry_string(symbols, pos, na_xyz))
    psi4.set_options({"basis": basis, "scf_type": "df", "guess": "sad", "freeze_core": True})
    psi4.energy("sapt0", molecule=mol)
    v = psi4.core.variables()
    return {k: v[name] * HARTREE_TO_KJ for k, name in COMPONENTS.items()}


def run_reference(out_path, solvents=SOLVENTS, n_dir=24, basis="aug-cc-pVDZ", threads=8, memory="8 GB"):
    """Append SAPT0 records to a JSON lines file, skipping points already there.

    The file is written point by point, so a job that hits its wall time can be
    resubmitted and continues where it stopped.
    """
    import psi4

    psi4.set_memory(memory)
    psi4.set_num_threads(threads)
    psi4.core.set_output_file("/dev/null", False)
    out_path = Path(out_path)
    done = set()
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            r = json.loads(line)
            done.add((r["solvent"], r["idx"]))
    with open(out_path, "a") as fh:
        for name in solvents:
            symbols, pos, _ = load_monomer(name)
            points = sample_na_positions(name, n_dir=n_dir)
            t0 = time.time()
            for idx, na in enumerate(points):
                if (name, idx) in done:
                    continue
                try:
                    e = sapt0_point(psi4, symbols, pos, na, basis)
                except Exception as exc:  # a failed point is logged, not fatal
                    print(f"skip {name} {idx}: {exc!r}"[:200], flush=True)
                    psi4.core.clean()
                    continue
                rec = {"solvent": name, "idx": idx, "na": [float(x) for x in na], "basis": basis, **e}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                psi4.core.clean()
            print(f"{name}: {len(points)} points, {time.time() - t0:.0f} s", flush=True)
