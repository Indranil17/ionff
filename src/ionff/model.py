"""Parameter container and vectorised energy evaluation for the four Na+ / site terms."""

from dataclasses import dataclass

import numpy as np

from .atomtypes import TYPES
from .terms import PAIR_FUNCS, TERMS

N_T = len(TYPES)

# (lower, upper, start) per parameter. A is stored as log(A) so that the search
# works on a sensible scale.
BOUNDS = {
    "q": (-1.5, 1.5, 0.0),
    "beta": (0.5, 6.0, 2.5),
    "logA": (2.0, 15.0, 10.0),
    "B": (1.5, 5.0, 3.5),
    "alpha": (0.0, 3.0, 0.5),
    "c": (0.5, 5.0, 2.0),
    "C6": (0.0, 5000.0, 200.0),
    "d": (0.5, 5.0, 2.0),
}
TERM_PARAMS = {"elst": ("q", "beta"), "exch": ("logA", "B"), "ind": ("alpha", "c"), "disp": ("C6", "d")}


@dataclass
class Params:
    """One value per atom type for every parameter name."""

    values: dict

    @classmethod
    def default(cls):
        return cls({k: np.full(N_T, v[2]) for k, v in BOUNDS.items()})

    def copy(self):
        return Params({k: v.copy() for k, v in self.values.items()})

    def pack(self, names, active):
        """Flatten the chosen parameter names over the active types."""
        return np.concatenate([self.values[n][active] for n in names])

    def unpack(self, x, names, active):
        out = self.copy()
        k = len(active)
        for i, n in enumerate(names):
            out.values[n][active] = x[i * k : (i + 1) * k]
        return out

    def to_dict(self):
        return {n: {t: float(v[i]) for i, t in enumerate(TYPES)} for n, v in self.values.items()}

    @classmethod
    def from_dict(cls, d):
        return cls({n: np.array([d[n][t] for t in TYPES]) for n in BOUNDS})


def bounds_for(names, active):
    lo = np.concatenate([np.full(len(active), BOUNDS[n][0]) for n in names])
    hi = np.concatenate([np.full(len(active), BOUNDS[n][1]) for n in names])
    return lo, hi


def term_energy(term, params, r, t):
    """Energy of one term for N configurations; r and t are the padded (N, S) arrays."""
    func = PAIR_FUNCS[term][0]
    names = TERM_PARAMS[term]
    tt = np.where(t >= 0, t, 0)
    a = params.values[names[0]][tt]
    b = params.values[names[1]][tt]
    if names[0] == "logA":
        a = np.exp(a)
    e = func(r, a, b)
    return np.where(t >= 0, e, 0.0).sum(axis=1)


def all_terms(params, r, t):
    return {k: term_energy(k, params, r, t) for k in TERMS}


def total_energy(params, r, t):
    return sum(all_terms(params, r, t).values())


def molecule_charges(params, counts):
    """Net charge of one molecule of every solvent, from its type counts."""
    return {n: float(c @ params.values["q"]) for n, c in counts.items()}
