"""GCN over molecular graphs (batched dense, padded)."""
from __future__ import annotations

import torch
import torch.nn as nn


class GCNLayer(nn.Module):
    def __init__(self, d_in: int, d_out: int):
        super().__init__()
        self.lin = nn.Linear(d_in, d_out)

    def forward(self, x: torch.Tensor, a_norm: torch.Tensor) -> torch.Tensor:
        return a_norm @ self.lin(x)


class MolGCN(nn.Module):
    """GraphDTA-style GCN encoder -> molecule embedding via masked mean pool."""

    def __init__(self, d_in: int, hidden: int = 64, out_dim: int = 128,
                 dropout: float = 0.2):
        super().__init__()
        self.gc1 = GCNLayer(d_in, hidden)
        self.gc2 = GCNLayer(hidden, hidden)
        self.proj = nn.Linear(hidden, out_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, a_norm: torch.Tensor,
                mask: torch.Tensor) -> torch.Tensor:
        # x: (b, n, f), a_norm: (b, n, n), mask: (b, n) 1=real atom
        h = torch.relu(self.gc1(x, a_norm))
        h = self.dropout(h)
        h = torch.relu(self.gc2(h, a_norm))
        h = h * mask.unsqueeze(-1)
        pooled = h.sum(1) / mask.sum(1, keepdim=True).clamp(min=1.0)
        return torch.relu(self.proj(pooled))


def collate_graphs(graphs: list):
    """Pad a list of (X, A) to batch tensors. Returns x, a_norm, mask."""
    b = len(graphs)
    n_max = max(X.shape[0] for X, _ in graphs)
    f = graphs[0][0].shape[1]
    x = torch.zeros(b, n_max, f)
    a = torch.zeros(b, n_max, n_max)
    mask = torch.zeros(b, n_max)
    for i, (X, A) in enumerate(graphs):
        n = X.shape[0]
        x[i, :n] = torch.tensor(X)
        a[i, :n, :n] = torch.tensor(A)
        mask[i, :n] = 1.0
    # symmetric norm per graph: D^-1/2 (A+I already inside A) D^-1/2
    deg = a.sum(-1)
    dinv = torch.where(deg > 0, deg, torch.ones_like(deg)).pow(-0.5)
    a_norm = a * dinv.unsqueeze(-1) * dinv.unsqueeze(-2)
    return x, a_norm, mask
