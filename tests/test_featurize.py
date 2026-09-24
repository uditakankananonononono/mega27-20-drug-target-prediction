import numpy as np
import pytest

from targetscan.data.featurize import (ATOM_SYMBOLS, N_ATOM_FEATS,
                                       encode_sequence, smiles_to_graph)


def test_aspirin_graph_handchecked():
    # aspirin: CC(=O)Oc1ccccc1C(=O)O - 13 heavy atoms, 13 bonds
    X, A = smiles_to_graph("CC(=O)Oc1ccccc1C(=O)O")
    assert X.shape == (13, N_ATOM_FEATS)
    assert A.shape == (13, 13)
    # adjacency symmetric, self loops present
    assert np.allclose(A, A.T)
    assert np.allclose(np.diag(A), 1.0)
    # bond count: (sum - self loops) / 2
    assert (A.sum() - 13) / 2 == 13
    # atom symbols: 9 carbons, 4 oxygens
    sym_idx = X[:, :len(ATOM_SYMBOLS)].argmax(1)
    assert (sym_idx == ATOM_SYMBOLS.index("C")).sum() == 9
    assert (sym_idx == ATOM_SYMBOLS.index("O")).sum() == 4
    # aromatic ring carbons flagged
    assert X[:, -1].sum() == 6


def test_methane_single_atom():
    X, A = smiles_to_graph("C")
    assert X.shape == (1, N_ATOM_FEATS)
    assert A.shape == (1, 1) and A[0, 0] == 1.0


def test_invalid_smiles_raises():
    with pytest.raises(ValueError):
        smiles_to_graph("not_a_molecule")


def test_encode_sequence_mapping_and_padding():
    enc = encode_sequence("ACD", max_len=8)
    assert enc.tolist()[:3] == [1, 2, 3]  # A, C, D in alphabet order
    assert (enc[3:] == 0).all()
    enc2 = encode_sequence("A" * 2000, max_len=1200)
    assert len(enc2) == 1200


def test_encode_unknown_residue_is_pad():
    enc = encode_sequence("AXA", max_len=4)
    assert enc[1] == 0
