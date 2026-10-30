import numpy as np

from ionff.terms import K_COUL, e_disp, e_elst, e_exch, e_ind, tt_damp


def test_damping_limits():
    assert tt_damp(np.array(0.0), 6) == 0.0
    assert abs(tt_damp(np.array(60.0), 6) - 1.0) < 1e-12


def test_elst_long_range_is_coulomb():
    r = np.array([30.0])
    assert abs(e_elst(r, 0.5, 2.0)[0] - K_COUL * 0.5 / 30.0) < 1e-6


def test_elst_is_finite_at_contact():
    assert np.isfinite(e_elst(np.array([1e-3]), -0.8, 3.0)[0])


def test_signs_and_decay():
    r = np.array([2.0, 3.0, 4.0])
    assert np.all(np.diff(e_exch(r, 1e4, 3.5)) < 0) and np.all(e_exch(r, 1e4, 3.5) > 0)
    assert np.all(e_ind(r, 0.5, 2.0) < 0) and np.all(e_disp(r, 100.0, 2.0) < 0)
    assert abs(e_ind(r, 0.5, 2.0)[2]) < abs(e_ind(r, 0.5, 2.0)[0])


def test_ind_long_range_scaling():
    r = np.array([20.0, 40.0])
    e = e_ind(r, 0.5, 2.0)
    assert abs(e[0] / e[1] - 16.0) < 1e-3  # 1 / r^4
