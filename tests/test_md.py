"""The OpenMM correction force must reproduce the analytic pair terms."""

import numpy as np
import pytest

pytest.importorskip("openmm")

from ionff.atomtypes import TYPE_INDEX  # noqa: E402
from ionff.md import SPCE_CHARGE, build_system, diffusion, msd_fft, rdf_and_cn, unwrap  # noqa: E402
from ionff.model import Params  # noqa: E402
from ionff.terms import K_COUL, e_disp, e_elst, e_exch, e_ind  # noqa: E402


def switch(r_nm, rs=0.9, rc=1.0):
    x = np.clip((r_nm - rs) / (rc - rs), 0.0, 1.0)
    return 1 - 6 * x**5 + 15 * x**4 - 10 * x**3


@pytest.mark.slow
def test_correction_force_matches_analytic_terms():
    import openmm as mm

    params = Params.default()
    for name, (o, h) in {"q": (-0.7, 0.35), "alpha": (1.2, 0.3), "C6": (300.0, 20.0), "logA": (11.0, 8.0)}.items():
        params.values[name][TYPE_INDEX["OW"]] = o
        params.values[name][TYPE_INDEX["HW"]] = h

    system, top, pos, na_idx = build_system(params, box_nm=2.0, n_na=1, seed=3)
    force = system.getForce(system.getNumForces() - 1)
    force.setForceGroup(1)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(pos)
    got = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)

    x = np.array(pos.value_in_unit(mm.unit.nanometer))
    box = 2.0
    expected = 0.0
    for atom in top.atoms():
        if atom.index in na_idx:
            continue
        sym = atom.element.symbol
        t = TYPE_INDEX["OW" if sym == "O" else "HW"]
        d = x[atom.index] - x[na_idx[0]]
        d -= box * np.round(d / box)
        r_nm = np.linalg.norm(d)
        if r_nm >= 1.0:
            continue
        r = 10.0 * r_nm
        v = params.values
        e = (e_elst(r, v["q"][t], v["beta"][t]) + e_exch(r, np.exp(v["logA"][t]), v["B"][t])
             + e_ind(r, v["alpha"][t], v["c"][t]) + e_disp(r, v["C6"][t], v["d"][t])
             - K_COUL * SPCE_CHARGE[sym] / r)
        expected += e * switch(r_nm)
    assert got == pytest.approx(expected, rel=1e-4, abs=1e-3)


def test_msd_fft_matches_brute_force():
    x = np.cumsum(np.random.default_rng(2).normal(size=(50, 3, 3)), axis=0)
    brute = [np.mean(np.sum((x[l:] - x[: len(x) - l]) ** 2, axis=-1)) for l in range(len(x))]
    assert np.allclose(msd_fft(x), brute)


def test_unwrap_and_diffusion_of_a_random_walk():
    rng = np.random.default_rng(0)
    box = 3.0
    d_true = 1.0e-9  # m^2/s
    stride = 0.5  # ps
    sigma = np.sqrt(2 * d_true * stride * 1e-12) * 1e9  # nm per axis per step
    steps = rng.normal(0.0, sigma, size=(4000, 40, 3))
    path = np.cumsum(steps, axis=0)
    wrapped = np.mod(path, box)
    assert np.allclose(unwrap(wrapped, box) - unwrap(wrapped, box)[0], path - path[0], atol=1e-9)
    d_raw, d_corr, err = diffusion(wrapped, box, stride)
    assert d_raw == pytest.approx(d_true, rel=0.25)
    assert d_corr > d_raw and err > 0


def test_rdf_of_uniform_distances_is_flat():
    rng = np.random.default_rng(1)
    box, n_w = 3.0, 900
    pts = rng.uniform(0, box, size=(200, 1, n_w, 3))
    ion = rng.uniform(0, box, size=(200, 1, 1, 3))
    d = pts - ion
    d -= box * np.round(d / box)
    dist = np.linalg.norm(d, axis=-1)[:, :, :]
    r, g, cn = rdf_and_cn(dist, box, r_max=1.0, n_bins=20)
    assert np.mean(g[5:]) == pytest.approx(1.0, abs=0.05)
    assert cn[-1] > cn[0]
