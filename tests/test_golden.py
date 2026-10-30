"""cxx/tests/golden.csv must still match the Python terms, so the C++ test checks the current code."""

from pathlib import Path

import numpy as np

from ionff.terms import e_disp, e_elst, e_exch, e_ind

GOLDEN = Path(__file__).resolve().parents[1] / "cxx" / "tests" / "golden.csv"


def test_golden_file_matches_python_terms():
    g = np.loadtxt(GOLDEN, delimiter=",", comments="#")
    r, q, beta, logA, B, alpha, c, C6, d = g[:, :9].T
    got = np.column_stack([e_elst(r, q, beta), e_exch(r, np.exp(logA), B), e_ind(r, alpha, c), e_disp(r, C6, d)])
    assert np.allclose(got, g[:, 9:], rtol=1e-12, atol=1e-15)
