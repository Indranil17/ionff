"""Force-field fitting: one SAPT term at a time, or the total energy only.

Both routes use the same functional forms, so any difference in how well they
generalise comes from what they are trained against. The search is a global
differential evolution run followed by a bounded least-squares polish.
"""

import numpy as np
from scipy.optimize import differential_evolution, least_squares

from .model import TERM_PARAMS, Params, bounds_for, molecule_charges, term_energy, total_energy
from .terms import TERMS

NEUTRALITY_WEIGHT = 20.0  # kJ/mol per elementary charge of net molecular charge


def _neutrality(params, counts):
    return NEUTRALITY_WEIGHT * np.array(list(molecule_charges(params, counts).values()))


def minimise(resid, lo, hi, x0, seed=0, popsize=12, maxiter=300):
    """Differential evolution on the sum of squares, then least squares from its best point."""

    def cost(x):
        r = resid(x)
        return float(r @ r)

    de = differential_evolution(
        cost, list(zip(lo, hi)), seed=seed, popsize=popsize, maxiter=maxiter, tol=1e-8,
        polish=False, init="sobol", x0=np.clip(x0, lo, hi),
    )
    ls = least_squares(resid, de.x, bounds=(lo, hi), x_scale="jac", max_nfev=400)
    return ls.x if ls.cost * 2 <= de.fun else de.x


def fit_term(ds, term, seed=0, **kw):
    """Fit the two parameters per active type of one term against its SAPT component."""
    active = ds.types_present()
    names = TERM_PARAMS[term]
    base = Params.default()
    ref = ds.energy[term]

    def resid(x):
        p = base.unpack(x, names, active)
        r = term_energy(term, p, ds.r, ds.t) - ref
        return np.concatenate([r, _neutrality(p, ds.counts)]) if term == "elst" else r

    lo, hi = bounds_for(names, active)
    x = minimise(resid, lo, hi, base.pack(names, active), seed=seed, **kw)
    return base.unpack(x, names, active)


def fit_terms(ds, seed=0, **kw):
    """Term-resolved fit: elst, exch, ind and disp each against their own reference."""
    out = Params.default()
    for term in TERMS:
        fitted = fit_term(ds, term, seed=seed, **kw)
        for n in TERM_PARAMS[term]:
            out.values[n] = fitted.values[n]
    return out


def fit_total(ds, seed=0, **kw):
    """Total-only fit: all eight parameters per active type against the SAPT0 total."""
    active = ds.types_present()
    names = [n for term in TERMS for n in TERM_PARAMS[term]]
    base = Params.default()
    ref = ds.energy["total"]

    def resid(x):
        p = base.unpack(x, names, active)
        return np.concatenate([total_energy(p, ds.r, ds.t) - ref, _neutrality(p, ds.counts)])

    lo, hi = bounds_for(names, active)
    x = minimise(resid, lo, hi, base.pack(names, active), seed=seed, **kw)
    return base.unpack(x, names, active)
