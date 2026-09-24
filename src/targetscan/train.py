"""Pair-minibatch training for the DTI model."""
from __future__ import annotations

import numpy as np
import torch

from .data.featurize import encode_sequence, smiles_to_graph
from .metrics import concordance_index_fast, mse
from .models.mol_gnn import collate_graphs


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)


class FeatureStore:
    """Precomputed drug graphs and protein encodings for a benchmark."""

    def __init__(self, smiles, sequences, max_len: int = 1200):
        self.graphs = [smiles_to_graph(s) for s in smiles]
        self.seqs = [torch.tensor(encode_sequence(s, max_len))
                     for s in sequences]

    def batch(self, pairs: np.ndarray):
        d_idx, t_idx = pairs[:, 0], pairs[:, 1]
        x, a_norm, mask = collate_graphs([self.graphs[i] for i in d_idx])
        seq = torch.stack([self.seqs[i] for i in t_idx])
        return x, a_norm, mask, seq


def train_dti(model, store: FeatureStore, pairs, y, val_pairs=None, val_y=None,
              epochs: int = 20, bs: int = 128, lr: float = 5e-4, wd: float = 1e-5,
              patience: int = 4, seed: int = 0, log=lambda *a: None,
              opt_state: dict = None):
    set_seed(seed)
    yt = torch.tensor(np.asarray(y), dtype=torch.float32)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    if opt_state is not None:
        opt.load_state_dict(opt_state)
        for g in opt.param_groups:  # allow lr schedule changes across chunks
            g["lr"] = lr
    best = {"loss": np.inf, "state": None, "stall": 0}
    n = len(pairs)
    for ep in range(epochs):
        model.train()
        perm = np.random.permutation(n)
        tot = 0.0
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            x, a_norm, mask, seq = store.batch(pairs[idx])
            pred = model(x, a_norm, mask, seq)
            # affinity-weighted MSE: strong binders (pKd>=7) are ~5% of pairs
            # but carry the discovery signal - upweight them
            w = torch.where(yt[idx] >= 7.0, torch.full_like(yt[idx], 6.0),
                            torch.ones_like(yt[idx]))
            loss = (w * (pred - yt[idx]) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * len(idx)
        train_mse = tot / n
        if val_pairs is not None:
            vm = evaluate(model, store, val_pairs, val_y, bs=bs)
            log(f"epoch {ep}: train_mse {train_mse:.4f} val_mse {vm['mse']:.4f} "
                f"val_ci {vm['ci']:.4f}")
            crit = vm["mse"]
        else:
            log(f"epoch {ep}: train_mse {train_mse:.4f}")
            crit = train_mse
        if crit < best["loss"] - 1e-4:
            best.update(loss=crit, stall=0,
                        state={k: v.detach().clone()
                               for k, v in model.state_dict().items()})
        else:
            best["stall"] += 1
            if best["stall"] >= patience:
                break
    if best["state"] is not None:
        model.load_state_dict(best["state"])
    return model, opt


@torch.no_grad()
def predict(model, store: FeatureStore, pairs, bs: int = 256) -> np.ndarray:
    model.eval()
    out = []
    for i in range(0, len(pairs), bs):
        x, a_norm, mask, seq = store.batch(pairs[i:i + bs])
        out.append(model(x, a_norm, mask, seq).numpy())
    return np.concatenate(out)


def evaluate(model, store: FeatureStore, pairs, y, bs: int = 256) -> dict:
    p = predict(model, store, pairs, bs)
    return {"mse": mse(y, p), "ci": concordance_index_fast(np.asarray(y), p),
            "n": int(len(y))}
