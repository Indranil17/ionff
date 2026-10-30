import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ionff.atomtypes import TYPE_INDEX
from ionff.dataset import Dataset
from ionff.geometry import load_monomer, sample_na_positions
from ionff.model import Params
from ionff.terms import TERMS


def truth_params():
    """A neutral, hand-built parameter set used to generate synthetic references."""
    p = Params.default()
    v = p.values
    for name, (q, beta, logA, B, alpha, c, C6, d) in {
        "OW": (-0.8, 3.0, 11.0, 3.6, 0.6, 2.0, 150.0, 2.0),
        "HW": (0.4, 2.0, 7.0, 3.2, 0.1, 2.5, 20.0, 2.5),
        "OC": (-0.5, 3.0, 10.5, 3.4, 0.5, 2.0, 120.0, 2.0),
        "OE": (-0.4, 3.0, 10.8, 3.5, 0.4, 2.0, 100.0, 2.0),
        "CC": (0.7, 2.5, 9.0, 3.3, 0.3, 2.2, 80.0, 2.2),
        "CA": (0.0, 2.5, 9.5, 3.3, 0.4, 2.2, 90.0, 2.2),
        "H": (0.0, 2.0, 6.5, 3.0, 0.1, 2.5, 10.0, 2.5),
    }.items():
        i = TYPE_INDEX[name]
        for key, val in zip(("q", "beta", "logA", "B", "alpha", "c", "C6", "d"), (q, beta, logA, B, alpha, c, C6, d)):
            v[key][i] = val
    return p


def synthetic_dataset(solvents=("water", "dme", "ec", "acetone"), n_dir=6):
    """Dataset whose references are the truth_params energies (no psi4 needed).

    Charges of CC/OC/OE do not cancel inside every molecule, so the neutrality
    penalty is switched off by giving truth charges to water only and zero
    elsewhere.
    """
    from ionff.model import all_terms

    p = truth_params()
    p.values["q"][[TYPE_INDEX[x] for x in ("OC", "OE", "CC", "CA", "H")]] = 0.0
    pts, rs, ts, names = [], [], [], []
    width = 12
    for s in solvents:
        _, pos, types = load_monomer(s)
        for na in sample_na_positions(s, n_dir=n_dir):
            d = np.linalg.norm(pos - na, axis=1)
            r = np.full(width, 60.0)
            t = np.full(width, -1)
            r[: len(d)] = d
            t[: len(d)] = [TYPE_INDEX[x] for x in types]
            rs.append(r), ts.append(t), names.append(s)
    r, t = np.array(rs), np.array(ts)
    terms = all_terms(p, r, t)
    energy = {**terms, "total": sum(terms.values())}
    keep = energy["total"] < 60.0
    counts = {}
    for s in solvents:
        _, _, types = load_monomer(s)
        counts[s] = np.bincount([TYPE_INDEX[x] for x in types], minlength=len(TYPE_INDEX))
    ds = Dataset(np.array(names)[keep], r[keep], t[keep], {k: v[keep] for k, v in energy.items()}, counts)
    return ds, p


@pytest.fixture(scope="session")
def synth():
    return synthetic_dataset()
