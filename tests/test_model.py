import numpy as np

from ionff.atomtypes import TYPE_INDEX
from ionff.model import Params, term_energy, total_energy
from conftest import truth_params


def test_pack_unpack_roundtrip():
    p = Params.default()
    active = [0, 2, 4]
    x = np.arange(6, dtype=float)
    q = p.unpack(x, ("q", "beta"), active)
    assert np.allclose(q.pack(("q", "beta"), active), x)
    assert q.values["q"][1] == p.values["q"][1]  # inactive type untouched


def test_params_dict_roundtrip():
    p = truth_params()
    assert np.allclose(Params.from_dict(p.to_dict()).values["logA"], p.values["logA"])


def test_padding_adds_nothing():
    p = truth_params()
    r = np.array([[2.3, 2.9, 3.1, 60.0, 60.0]])
    t = np.array([[TYPE_INDEX["OW"], TYPE_INDEX["HW"], TYPE_INDEX["HW"], -1, -1]])
    short = total_energy(p, r[:, :3], t[:, :3])
    assert np.allclose(short, total_energy(p, r, t))


def test_single_site_matches_pair_function():
    from ionff.terms import e_exch

    p = truth_params()
    i = TYPE_INDEX["OW"]
    e = term_energy("exch", p, np.array([[2.4]]), np.array([[i]]))
    assert np.isclose(e[0], e_exch(2.4, np.exp(p.values["logA"][i]), p.values["B"][i]))
