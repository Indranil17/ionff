"""Error metrics, leave-one-solvent-out testing and the Na+(H2O)n cluster test."""

import numpy as np

from .model import all_terms
from .terms import TERMS

KEYS = (*TERMS, "total")


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def mae(a, b):
    return float(np.mean(np.abs(np.asarray(a) - np.asarray(b))))


def predict(params, r, t):
    """Per-term energies plus their sum under the key 'total'."""
    terms = all_terms(params, r, t)
    terms["total"] = sum(terms.values())
    return terms


def score(params, ds):
    """RMSE and MAE (kJ/mol) of every term and the total on a dataset."""
    pred = predict(params, ds.r, ds.t)
    return {k: {"rmse": rmse(pred[k], ds.energy[k]), "mae": mae(pred[k], ds.energy[k])} for k in KEYS}


def score_by_solvent(params, ds):
    return {s: score(params, ds.subset(ds.solvent == s)) for s in ds.solvents()}


def neutralise(params, counts):
    """Copy of `params` with one molecule made neutral by a common shift of its type charges.

    `counts` is the (n_types,) type count of the molecule. Every type the molecule contains
    gets the same shift, -Q / n_atoms with Q the net charge, so the molecule ends at zero
    net charge. The copy is meant for scoring that one molecule only, since the other
    molecules sharing those types would see the shifted charges too.
    """
    out = params.copy()
    q = out.values["q"]
    out.values["q"] = np.where(counts > 0, q - float(counts @ q) / counts.sum(), q)
    return out


def loso(ds, fitters, holdable=None, neutral=False, **kw):
    """Leave one solvent out: fit on the rest, score on the held-out solvent.

    `fitters` maps a label to a function fit(ds, **kw) returning Params. A solvent
    is held out only if every atom type it contains also appears in the remaining
    training data; otherwise its parameters could not be learned and the test
    would say nothing about transfer.

    The charge neutrality penalty acts on training molecules only, so a held-out
    molecule can come out charged. With `neutral=True` the held-out molecule is made
    neutral with `neutralise` before scoring, and its net charge before the shift is
    stored under `<label>_net_charge`.
    """
    results = []
    for s in holdable or ds.solvents():
        train = ds.subset(ds.solvent != s)
        test = ds.subset(ds.solvent == s)
        if not set(test.types_present()) <= set(train.types_present()):
            results.append({"held_out": s, "skipped": "atom types not covered by the remaining solvents"})
            continue
        row = {"held_out": s, "n_train": len(train), "n_test": len(test)}
        for label, fit in fitters.items():
            fitted = fit(train, **kw)
            if neutral:
                row[f"{label}_net_charge"] = float(test.counts[s] @ fitted.values["q"])
                fitted = neutralise(fitted, test.counts[s])
            row[label] = score(fitted, test)
        results.append(row)
    return results


def cluster_scores(params, records):
    """Compare the additive force field with SAPT0 for Na+ against a water cluster.

    Each record carries `symbols`, `xyz` (cluster atoms), `na`, the SAPT0 terms and
    the cluster size `n`. The force field sums the Na+ to site terms over every
    water in the cluster, so the error combines non-additivity with any misfit of
    the single-water terms. The two are not separated here.
    """
    from .atomtypes import TYPE_INDEX, assign_types

    rows = []
    for rec in records:
        xyz = np.asarray(rec["xyz"])
        types = []
        for i in range(0, len(xyz), 3):
            types += assign_types(rec["symbols"][i : i + 3], xyz[i : i + 3])
        r = np.linalg.norm(xyz - np.asarray(rec["na"]), axis=1)[None, :]
        t = np.array([[TYPE_INDEX[x] for x in types]])
        pred = predict(params, r, t)
        row = {"n": rec["n"]}
        for k in KEYS:
            row[k] = {"ff": float(pred[k][0]), "sapt0": rec[k]}
        rows.append(row)
    return rows
