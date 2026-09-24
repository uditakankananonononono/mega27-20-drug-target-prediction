"""DTI model: MolGCN + ProtCNN -> affinity regression."""
from __future__ import annotations

import torch
import torch.nn as nn

from .mol_gnn import MolGCN
from .prot_cnn import ProtCNN


class DTINet(nn.Module):
    def __init__(self, n_atom_feats: int, out_dim: int = 128,
                 dropout: float = 0.2):
        super().__init__()
        self.mol = MolGCN(n_atom_feats, out_dim=out_dim, dropout=dropout)
        self.prot = ProtCNN(out_dim=out_dim, dropout=dropout)
        self.head = nn.Sequential(
            nn.Linear(2 * out_dim, 256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 1),
        )

    def forward(self, x, a_norm, mask, seq) -> torch.Tensor:
        z = torch.cat([self.mol(x, a_norm, mask), self.prot(seq)], dim=1)
        return self.head(z).squeeze(-1)
