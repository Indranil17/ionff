"""Na+(H2O)n clusters: starting geometries, relaxation and SAPT0 against the water cluster."""

import json
from pathlib import Path

import numpy as np

from .geometry import fibonacci_sphere, load_monomer
from .reference import COMPONENTS, HARTREE_TO_KJ

ANGSTROM_PER_BOHR = 0.52917721067


def _align(a, b):
    """Rotation matrix taking unit vector a onto unit vector b (Rodrigues)."""
    v = np.cross(a, b)
    c = float(a @ b)
    if np.linalg.norm(v) < 1e-10:
        return np.eye(3) if c > 0 else -np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k / (1.0 + c)


def build_cluster(n, r_na_o=2.4):
    """Na+ at the origin with n waters on a spherical spread, oxygen towards the ion.

    Returns (symbols, positions) of the waters only, three atoms each (O, H, H).
    """
    sym, pos, _ = load_monomer("water")
    axis = pos[1:].mean(axis=0) - pos[0]
    axis /= np.linalg.norm(axis)
    atoms, symbols = [], []
    for u in fibonacci_sphere(n):
        rot = _align(axis, u)
        atoms.append((pos - pos[0]) @ rot.T + r_na_o * u)
        symbols += sym
    return symbols, np.vstack(atoms)


def _psi4_geometry(symbols, xyz, na):
    lines = [f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, (x, y, z) in zip(symbols, xyz)]
    return lines, "Na {:.8f} {:.8f} {:.8f}".format(*na)


def relax_cluster(psi4, symbols, xyz, maxiter=200):
    """B3LYP/def2-SVP optimisation of Na+(H2O)n; returns (water xyz, Na xyz in Angstrom, converged).

    Floppy hydrogen bond networks sometimes need more steps than psi4 allows by default. If the limit is
    reached the last geometry is kept and flagged, since SAPT0 against the force field is still a valid
    comparison at any structure.
    """
    lines, na_line = _psi4_geometry(symbols, xyz, np.zeros(3))
    mol = psi4.geometry("\n".join(["1 1", *lines, na_line, "units angstrom", "no_com", "no_reorient", "symmetry c1"]))
    psi4.set_options({"basis": "def2-svp", "scf_type": "df", "g_convergence": "gau", "geom_maxiter": maxiter})
    converged = True
    try:
        psi4.optimize("b3lyp", molecule=mol)
    except psi4.OptimizationConvergenceError as err:
        converged = False
        mol = err.wfn.molecule()
    g = mol.geometry().np * ANGSTROM_PER_BOHR
    return g[:-1], g[-1], converged


def sapt0_cluster(psi4, symbols, xyz, na, basis="aug-cc-pVDZ"):
    """SAPT0 of the water cluster (fragment A) against Na+ (fragment B), kJ/mol."""
    lines, na_line = _psi4_geometry(symbols, xyz, na)
    geom = "\n".join(["0 1", *lines, "--", "1 1", na_line, "units angstrom", "no_com", "no_reorient", "symmetry c1"])
    mol = psi4.geometry(geom)
    psi4.set_options({"basis": basis, "scf_type": "df", "guess": "sad", "freeze_core": True})
    psi4.energy("sapt0", molecule=mol)
    v = psi4.core.variables()
    return {k: v[name] * HARTREE_TO_KJ for k, name in COMPONENTS.items()}


def run_clusters(out_path, sizes=(1, 2, 3, 4, 5, 6), n_perturbed=2, basis="aug-cc-pVDZ", threads=8, memory="8 GB", seed=0):
    """Write relaxed and perturbed Na+(H2O)n SAPT0 records to a JSON lines file."""
    import psi4

    psi4.set_memory(memory)
    psi4.set_num_threads(threads)
    psi4.core.set_output_file("/dev/null", False)
    rng = np.random.default_rng(seed)
    out_path = Path(out_path)
    done = set()
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            r = json.loads(line)
            done.add((r["n"], r["kind"], r["idx"]))
    with open(out_path, "a") as fh:
        for n in sizes:
            symbols, start = build_cluster(n)
            xyz, na, converged = None, None, True
            todo = [("relaxed", 0)] + [("perturbed", i) for i in range(n_perturbed)]
            for kind, idx in todo:
                if (n, kind, idx) in done:
                    continue
                if xyz is None:
                    xyz, na, converged = relax_cluster(psi4, symbols, start)
                use = xyz if kind == "relaxed" else xyz + rng.normal(scale=0.08, size=xyz.shape)
                e = sapt0_cluster(psi4, symbols, use, na, basis)
                rec = {"n": n, "kind": kind, "idx": idx, "symbols": symbols, "xyz": use.tolist(), "na": na.tolist(), "converged": converged, **e}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                psi4.core.clean()
            print(f"n={n} done", flush=True)
