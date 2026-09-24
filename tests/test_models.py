import numpy as np
import torch

from targetscan.data.featurize import encode_sequence, smiles_to_graph
from targetscan.models.dti import DTINet
from targetscan.models.mol_gnn import MolGCN, collate_graphs
from targetscan.models.prot_cnn import ProtCNN
from targetscan.metrics import concordance_index, concordance_index_fast, mse


def _batch(smiles_list):
    return collate_graphs([smiles_to_graph(s) for s in smiles_list])


def test_collate_padding_and_norm():
    x, a_norm, mask = _batch(["C", "CCO"])
    assert x.shape[0] == 2 and x.shape[1] == 3  # padded to max atoms
    assert mask.tolist() == [[1, 0, 0], [1, 1, 1]]
    # padded rows of a_norm are zero (no NaN)
    assert torch.isfinite(a_norm).all()
    assert torch.allclose(a_norm[0, 1:, :], torch.zeros(2, 3))


def test_molgcn_forward_and_grad():
    x, a_norm, mask = _batch(["CCO", "c1ccccc1"])
    net = MolGCN(x.shape[-1], hidden=16, out_dim=8)
    z = net(x, a_norm, mask)
    assert z.shape == (2, 8)
    z.sum().backward()
    assert all(p.grad is not None for p in net.parameters())


def test_molgcn_isomorphism_invariance():
    # same molecule, different atom order in SMILES -> same embedding
    x1, a1, m1 = _batch(["CCO"])
    x2, a2, m2 = _batch(["OCC"])
    net = MolGCN(x1.shape[-1], hidden=16, out_dim=8)
    net.eval()
    with torch.no_grad():
        z1 = net(x1, a1, m1)
        z2 = net(x2, a2, m2)
    assert torch.allclose(z1, z2, atol=1e-4)


def test_protcnn_forward():
    net = ProtCNN(out_dim=8)
    seq = torch.tensor([encode_sequence("MKTAYIAK"), encode_sequence("ACDEFGH")])
    z = net(seq)
    assert z.shape == (2, 8)


def test_dtinet_end_to_end():
    store_graphs = [smiles_to_graph(s) for s in ["CCO", "c1ccccc1", "CC(=O)O"]]
    x, a_norm, mask = collate_graphs(store_graphs)
    seq = torch.tensor([encode_sequence("MKT"), encode_sequence("ACD"),
                        encode_sequence("GGG")])
    net = DTINet(x.shape[-1], out_dim=8)
    out = net(x, a_norm, mask, seq)
    assert out.shape == (3,)


def test_metrics_hand_computed():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    p = np.array([1.0, 2.0, 3.0, 4.0])
    assert concordance_index(y, p) == 1.0
    assert concordance_index(y, p[::-1]) == 0.0
    assert mse(y, p) == 0.0
    rng = np.random.default_rng(0)
    yr = rng.normal(size=300)
    pr = yr + 0.3 * rng.normal(size=300)
    assert abs(concordance_index(yr, pr) - concordance_index_fast(yr, pr)) < 1e-9


def test_dtinet_can_overfit_tiny_batch():
    torch.manual_seed(0)
    graphs = [smiles_to_graph(s) for s in ["CCO", "c1ccccc1"] * 4]
    x, a_norm, mask = collate_graphs(graphs)
    seq = torch.tensor([encode_sequence(s) for s in
                        ["MKT", "ACD", "GGG", "WWW"] * 2])
    y = torch.tensor([5.0, 7.0, 6.0, 8.0] * 2)
    net = DTINet(x.shape[-1], out_dim=16)
    opt = torch.optim.Adam(net.parameters(), lr=5e-3)
    loss0 = None
    for _ in range(300):
        pred = net(x, a_norm, mask, seq)
        loss = torch.nn.functional.mse_loss(pred, y)
        if loss0 is None:
            loss0 = loss.item()
        opt.zero_grad(); loss.backward(); opt.step()
    assert loss.item() < 0.05 * loss0
