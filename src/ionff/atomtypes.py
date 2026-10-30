"""Rule-based atom typing for the small set of solvents used here."""

import numpy as np
from ase.data import atomic_numbers, covalent_radii

TYPES = ("OW", "HW", "OH", "HO", "OC", "OE", "CC", "CA", "H")
TYPE_INDEX = {t: i for i, t in enumerate(TYPES)}

BOND_SCALE = 1.25


def bond_matrix(symbols, positions):
    """Boolean connectivity from covalent radii (BOND_SCALE times the radius sum)."""
    pos = np.asarray(positions, float)
    rc = covalent_radii[[atomic_numbers[s] for s in symbols]]
    d = np.linalg.norm(pos[:, None] - pos[None], axis=-1)
    bonded = d < BOND_SCALE * (rc[:, None] + rc[None])
    np.fill_diagonal(bonded, False)
    return bonded


def assign_types(symbols, positions):
    """Return one type label per atom of a neutral closed-shell solvent molecule.

    OW, HW: water oxygen and hydrogen. OH, HO: hydroxyl oxygen and its hydrogen.
    OC: carbonyl oxygen (one neighbour). OE: ether or ester oxygen (two heavy
    neighbours). CC: carbon with a carbonyl oxygen or three oxygen neighbours.
    CA: any other carbon. H: hydrogen not bonded to oxygen.
    """
    symbols = list(symbols)
    bonded = bond_matrix(symbols, positions)
    nb = [np.flatnonzero(bonded[i]) for i in range(len(symbols))]
    out = []
    for i, s in enumerate(symbols):
        if s == "O":
            n_h = sum(symbols[j] == "H" for j in nb[i])
            if n_h == 2:
                out.append("OW")
            elif n_h == 1:
                out.append("OH")
            elif len(nb[i]) == 1:
                out.append("OC")
            elif len(nb[i]) == 2:
                out.append("OE")
            else:
                raise ValueError(f"cannot type oxygen {i} with {len(nb[i])} neighbours")
        elif s == "H":
            host = [symbols[j] for j in nb[i]]
            if "O" in host:
                o = nb[i][host.index("O")]
                n_h = sum(symbols[k] == "H" for k in nb[o])
                out.append("HW" if n_h == 2 else "HO")
            else:
                out.append("H")
        elif s == "C":
            carbonyl = any(symbols[j] == "O" and len(nb[j]) == 1 for j in nb[i])
            n_o = sum(symbols[j] == "O" for j in nb[i])
            out.append("CC" if carbonyl or n_o >= 3 else "CA")
        else:
            raise ValueError(f"no type for element {s}")
    return out
