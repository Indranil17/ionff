import numpy as np
import pytest

from ionff.evaluate import loso, neutralise, rmse, score
from ionff.fit import fit_term, fit_terms, fit_total
from ionff.model import term_energy

FAST = dict(popsize=6, maxiter=60)


def test_exchange_term_recovers_reference(synth):
    ds, _ = synth
    p = fit_term(ds, "exch", **FAST)
    assert rmse(term_energy("exch", p, ds.r, ds.t), ds.energy["exch"]) < 1.0


def test_term_resolved_fit_recovers_every_term(synth):
    ds, _ = synth
    s = score(fit_terms(ds, **FAST), ds)
    # induction and electrostatics are the slow terms to pin down with a short search
    assert s["exch"]["rmse"] < 1.0 and s["disp"]["rmse"] < 0.5
    assert s["total"]["rmse"] < 4.0


def test_total_fit_reaches_low_total_error(synth):
    ds, _ = synth
    s = score(fit_total(ds, **FAST), ds)
    assert s["total"]["rmse"] < 4.0


def test_loso_skips_uncovered_solvents(synth):
    ds, _ = synth
    res = loso(ds, {"terms": fit_terms}, **FAST)
    by = {r["held_out"]: r for r in res}
    assert "skipped" in by["water"]  # OW and HW occur only in water
    assert "terms" in by["acetone"]  # OC, CC, CA, H all covered by EC and DME


def test_neutralise_zeroes_net_charge_of_one_molecule(synth):
    ds, _ = synth
    p = fit_term(ds, "elst", **FAST)
    counts = ds.counts["acetone"]
    before = p.values["q"].copy()
    q = neutralise(p, counts).values["q"]
    assert abs(counts @ q) < 1e-9
    assert np.array_equal(q[counts == 0], before[counts == 0])
    assert np.array_equal(p.values["q"], before)  # the input is not modified


def test_loso_neutral_records_net_charge_before_the_shift(synth):
    ds, _ = synth
    res = loso(ds, {"terms": fit_terms}, neutral=True, **FAST)
    row = {r["held_out"]: r for r in res}["acetone"]
    assert "terms_net_charge" in row and np.isfinite(row["terms_net_charge"])
    plain = {r["held_out"]: r for r in loso(ds, {"terms": fit_terms}, **FAST)}["acetone"]
    assert "terms_net_charge" not in plain


def test_cluster_n1_matches_scan_path(synth):
    """A one-water cluster record must give the same energies as the scan evaluation."""
    import numpy as np

    from ionff.evaluate import cluster_scores, predict
    from ionff.geometry import load_monomer
    from conftest import truth_params

    _, pos, _ = load_monomer("water")
    sym = ["O", "H", "H"]
    na = np.array([0.1, 0.2, 2.5])
    p = truth_params()
    rec = {"n": 1, "symbols": sym, "xyz": pos.tolist(), "na": na.tolist(),
           **{k: 0.0 for k in ("elst", "exch", "ind", "disp", "total")}}
    row = cluster_scores(p, [rec])[0]
    from ionff.atomtypes import TYPE_INDEX

    r = np.linalg.norm(pos - na, axis=1)[None]
    t = np.array([[TYPE_INDEX["OW"], TYPE_INDEX["HW"], TYPE_INDEX["HW"]]])
    ref = predict(p, r, t)
    assert np.isclose(row["total"]["ff"], ref["total"][0])
