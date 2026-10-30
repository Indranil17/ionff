"""Analytic Na+ / solvent-site interaction terms.

Units: distances in Angstrom, energies in kJ/mol, charges in e, polarisabilities
in Angstrom^3. Every term is a sum over solvent sites j of a pair function of
the Na+ to site distance r_j. The four terms mirror the four SAPT0 components,
so each can be fitted against its own reference.
"""

import numpy as np

K_COUL = 1389.35457520  # kJ mol^-1 Angstrom e^-2, 1 / (4 pi eps0)
TERMS = ("elst", "exch", "ind", "disp")
Q_NA = 1.0


def tt_damp(x, n):
    """Tang-Toennies damping, 1 - exp(-x) * sum_{k<=n} x^k / k!."""
    x = np.asarray(x, float)
    total = np.ones_like(x)
    term = np.ones_like(x)
    for k in range(1, n + 1):
        term = term * x / k
        total = total + term
    return 1.0 - np.exp(-x) * total


def e_elst(r, q, beta):
    """Coulomb energy of Na+ with a Slater-smeared site charge q."""
    return K_COUL * Q_NA * q / r * (1.0 - np.exp(-beta * r) * (1.0 + beta * r))


def e_exch(r, A, B):
    """Born-Mayer exchange repulsion."""
    return A * np.exp(-B * r)


def e_ind(r, alpha, c):
    """Charge to induced-dipole energy, -1/2 alpha q^2 / r^4, damped."""
    return -0.5 * K_COUL * Q_NA**2 * alpha / r**4 * tt_damp(c * r, 4)


def e_disp(r, C6, d):
    """Damped London dispersion, -C6 / r^6."""
    return -C6 / r**6 * tt_damp(d * r, 6)


PAIR_FUNCS = {
    "elst": (e_elst, ("q", "beta")),
    "exch": (e_exch, ("A", "B")),
    "ind": (e_ind, ("alpha", "c")),
    "disp": (e_disp, ("C6", "d")),
}

PAD_R = 60.0  # padded sites sit far away and carry zero parameters
