import numpy as np
import pytest

from ionff.clusters import build_cluster


@pytest.mark.parametrize("n", [1, 2, 4, 6])
def test_cluster_geometry(n):
    sym, xyz = build_cluster(n)
    assert len(sym) == 3 * n and xyz.shape == (3 * n, 3)
    for i in range(n):
        o = xyz[3 * i]
        assert abs(np.linalg.norm(o) - 2.4) < 1e-6  # O to Na+ at the origin
        assert np.linalg.norm(xyz[3 * i + 1] - o) < 1.05  # intact O-H bond


def test_waters_do_not_overlap():
    _, xyz = build_cluster(6)
    ox = xyz[::3]
    d = np.linalg.norm(ox[:, None] - ox[None], axis=-1) + np.eye(6) * 9
    assert d.min() > 2.5
