"""Solvent monomers and Na+ placement for the reference scans."""

from pathlib import Path

import numpy as np
from ase.io import read

from .atomtypes import assign_types, bond_matrix

MONOMER_DIR = Path(__file__).resolve().parents[2] / "data" / "monomers"
SOLVENTS = ("water", "methanol", "dme", "ec", "dmc", "acetone")

# Nearest-atom distances (Angstrom) of the radial scan. The short ones sample the
# wall and the minimum, the long ones the tail where induction and dispersion matter.
SCAN_DISTANCES = (2.0, 2.2, 2.5, 3.0, 3.6, 4.5)
MIN_CLEARANCE = 1.7  # no Na+ closer than this to any atom


def load_monomer(name):
    """Return (symbols, positions, types) of a solvent monomer."""
    atoms = read(MONOMER_DIR / f"{name}.xyz")
    sym = atoms.get_chemical_symbols()
    return sym, atoms.positions.copy(), assign_types(sym, atoms.positions)


def fibonacci_sphere(n):
    """n roughly uniform unit vectors."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5**0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def _ray_point(pos, origin, direction, target):
    """Point on the ray origin + t*direction whose nearest atom is at `target`.

    Walks outward in 0.02 Angstrom steps from inside the molecule and returns the
    first point at or beyond `target`; None if the ray never gets there.
    """
    for t in np.arange(0.0, 12.0, 0.02):
        p = origin + t * direction
        if np.linalg.norm(pos - p, axis=1).min() >= target:
            return p
    return None


def lone_pair_directions(symbols, pos):
    """Unit vectors pointing away from each oxygen, opposite its neighbours."""
    bonded = bond_matrix(symbols, pos)
    out = []
    for i, s in enumerate(symbols):
        if s != "O":
            continue
        v = pos[i] - pos[np.flatnonzero(bonded[i])].mean(axis=0)
        out.append((i, v / np.linalg.norm(v)))
    return out


def sample_na_positions(name, n_dir=24, n_cone=2, cone_deg=25.0, seed=0):
    """Na+ positions around a fixed monomer.

    n_dir roughly uniform directions from the centre of mass, each scanned over
    SCAN_DISTANCES, plus n_cone jittered lone-pair directions per oxygen scanned
    the same way. Returns an (N, 3) array.
    """
    rng = np.random.default_rng(seed)
    sym, pos, _ = load_monomer(name)
    com = pos.mean(axis=0)
    rays = [(com, d) for d in fibonacci_sphere(n_dir)]
    for i, u in lone_pair_directions(sym, pos):
        for _ in range(n_cone):
            v = u + np.tan(np.radians(cone_deg)) * rng.normal(size=3) / 3**0.5
            rays.append((pos[i], v / np.linalg.norm(v)))
    pts = []
    for origin, d in rays:
        for target in SCAN_DISTANCES:
            p = _ray_point(pos, origin, d, target)
            if p is not None and np.linalg.norm(pos - p, axis=1).min() >= MIN_CLEARANCE:
                pts.append(p)
    return np.array(pts)
