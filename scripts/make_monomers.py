"""Build the solvent monomer geometries.

RDKit embeds each molecule from SMILES (ETKDG plus MMFF), then Psi4 relaxes it at
B3LYP/def2-SVP. The relaxed geometries are written to data/monomers/*.xyz.
Run once: python scripts/make_monomers.py
"""

import sys
from pathlib import Path

import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

SOLVENTS = {
    "water": "O",
    "methanol": "CO",
    "dme": "COC",
    "ec": "O=C1OCCO1",
    "dmc": "COC(=O)OC",
    "acetone": "CC(C)=O",
}
OUT = Path(__file__).resolve().parents[1] / "data" / "monomers"
LEVEL = "b3lyp/def2-svp"


def embed(smiles, seed=7):
    """Lowest-energy MMFF conformer out of 30 ETKDG embeddings, as (symbols, xyz)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    cids = AllChem.EmbedMultipleConfs(mol, 30, randomSeed=seed)
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    best = min(range(len(cids)), key=lambda i: res[i][1])
    conf = mol.GetConformer(cids[best])
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return symbols, xyz


def relax(symbols, xyz):
    """Psi4 geometry optimisation; returns relaxed Angstrom coordinates."""
    lines = "\n".join(f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, (x, y, z) in zip(symbols, xyz))
    mol = psi4.geometry(f"0 1\n{lines}\nunits angstrom\nsymmetry c1\nno_com\nno_reorient")
    psi4.set_options({"basis": "def2-svp", "scf_type": "df", "g_convergence": "gau_tight"})
    psi4.optimize("b3lyp", molecule=mol)
    return mol.geometry().np * 0.52917721067


def write_xyz(path, symbols, xyz, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(symbols)}\n{comment}\n")
        for s, (x, y, z) in zip(symbols, xyz):
            fh.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


if __name__ == "__main__":
    psi4.set_memory("4 GB")
    psi4.set_num_threads(8)
    psi4.core.set_output_file("/dev/null", False)
    OUT.mkdir(parents=True, exist_ok=True)
    names = sys.argv[1:] or list(SOLVENTS)
    for name in names:
        symbols, xyz = embed(SOLVENTS[name])
        xyz = relax(symbols, xyz)
        write_xyz(OUT / f"{name}.xyz", symbols, xyz, f"{name} {LEVEL}")
        print("wrote", OUT / f"{name}.xyz", flush=True)
