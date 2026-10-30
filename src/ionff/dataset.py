"""Reference scans as padded arrays that the force-field terms can evaluate in one pass."""

import json
from dataclasses import dataclass

import numpy as np

from .atomtypes import TYPE_INDEX, TYPES
from .geometry import load_monomer
from .terms import PAD_R, TERMS

E_MAX = 60.0  # kJ/mol, drop deeper into the repulsive wall than this


@dataclass
class Dataset:
    solvent: np.ndarray  # (N,) solvent name per point
    r: np.ndarray  # (N, S) Na+ to site distance, padded with PAD_R
    t: np.ndarray  # (N, S) site type index, -1 for padding
    energy: dict  # term name or "total" to (N,) reference energies, kJ/mol
    counts: dict  # solvent name to (n_types,) type counts of one molecule

    def __len__(self):
        return len(self.solvent)

    def subset(self, mask):
        mask = np.asarray(mask, bool)
        names = set(self.solvent[mask])
        return Dataset(
            self.solvent[mask],
            self.r[mask],
            self.t[mask],
            {k: v[mask] for k, v in self.energy.items()},
            {k: v for k, v in self.counts.items() if k in names},
        )

    def solvents(self):
        return sorted(set(self.solvent))

    def types_present(self):
        return sorted(set(self.t[self.t >= 0].tolist()))


def _padded(distances, types, width):
    r = np.full(width, PAD_R)
    t = np.full(width, -1, dtype=int)
    r[: len(distances)] = distances
    t[: len(types)] = [TYPE_INDEX[x] for x in types]
    return r, t


def from_records(records, e_max=E_MAX):
    """Build a Dataset from reference records (dicts as written by reference.run_reference)."""
    monomers = {n: load_monomer(n) for n in {r["solvent"] for r in records}}
    width = max(len(m[1]) for m in monomers.values())
    rs, ts, names = [], [], []
    en = {k: [] for k in (*TERMS, "total")}
    for rec in records:
        if rec["total"] > e_max:
            continue
        _, pos, types = monomers[rec["solvent"]]
        d = np.linalg.norm(pos - np.asarray(rec["na"]), axis=1)
        r, t = _padded(d, types, width)
        rs.append(r)
        ts.append(t)
        names.append(rec["solvent"])
        for k in en:
            en[k].append(rec[k])
    counts = {n: np.bincount([TYPE_INDEX[x] for x in m[2]], minlength=len(TYPES)) for n, m in monomers.items()}
    return Dataset(
        np.array(names),
        np.array(rs),
        np.array(ts),
        {k: np.array(v) for k, v in en.items()},
        counts,
    )


def load_jsonl(path, e_max=E_MAX):
    with open(path) as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    return from_records(records, e_max)
