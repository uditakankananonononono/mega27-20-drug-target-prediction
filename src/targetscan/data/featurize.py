"""Molecule and protein featurization. Pure functions over RDKit / strings."""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem

ATOM_SYMBOLS = ["C", "N", "O", "S", "F", "P", "Cl", "Br", "I", "H", "OTHER"]
HYBRIDIZATIONS = [rdchem.HybridizationType.SP,
                  rdchem.HybridizationType.SP2,
                  rdchem.HybridizationType.SP3,
                  rdchem.HybridizationType.SP3D,
                  rdchem.HybridizationType.SP3D2]
N_ATOM_FEATS = len(ATOM_SYMBOLS) + len(HYBRIDIZATIONS) + 5

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard AAs; others -> 0 (UNK)
AA_TO_INT = {aa: i + 1 for i, aa in enumerate(AA_ALPHABET)}


def atom_features(atom: rdchem.Atom) -> np.ndarray:
    sym = atom.GetSymbol()
    v = np.zeros(N_ATOM_FEATS, dtype=np.float32)
    idx = ATOM_SYMBOLS.index(sym) if sym in ATOM_SYMBOLS else len(ATOM_SYMBOLS) - 1
    v[idx] = 1.0
    hyb = atom.GetHybridization()
    if hyb in HYBRIDIZATIONS:
        v[len(ATOM_SYMBOLS) + HYBRIDIZATIONS.index(hyb)] = 1.0
    off = len(ATOM_SYMBOLS) + len(HYBRIDIZATIONS)
    v[off + 0] = float(atom.GetDegree())
    v[off + 1] = float(atom.GetTotalNumHs())
    v[off + 2] = float(atom.GetFormalCharge())
    v[off + 3] = float(atom.GetTotalValence())
    v[off + 4] = 1.0 if atom.GetIsAromatic() else 0.0
    return v


def smiles_to_graph(smiles: str) -> tuple:
    """SMILES -> (X: (n, N_ATOM_FEATS) float32, A: (n, n) float32 with self loops)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"unparsable SMILES: {smiles!r}")
    n = mol.GetNumAtoms()
    X = np.stack([atom_features(mol.GetAtomWithIdx(i)) for i in range(n)])
    A = np.zeros((n, n), dtype=np.float32)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        A[i, j] = A[j, i] = 1.0
    A += np.eye(n, dtype=np.float32)
    return X, A


def encode_sequence(seq: str, max_len: int = 1200) -> np.ndarray:
    """Amino-acid sequence -> int64 (max_len,) padded with 0 (UNK/pad)."""
    out = np.zeros(max_len, dtype=np.int64)
    for i, ch in enumerate(seq[:max_len]):
        out[i] = AA_TO_INT.get(ch.upper(), 0)
    return out
