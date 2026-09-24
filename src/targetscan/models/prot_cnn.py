"""1D-CNN protein-sequence encoder (DeepDTA-style)."""
from __future__ import annotations

import torch
import torch.nn as nn


class ProtCNN(nn.Module):
    def __init__(self, vocab: int = 21, embed: int = 64,
                 channels: tuple = (64, 96), kernel: int = 8,
                 out_dim: int = 128, dropout: float = 0.2):
        super().__init__()
        self.embed = nn.Embedding(vocab, embed, padding_idx=0)
        blocks, prev = [], embed
        for c in channels:
            blocks += [nn.Conv1d(prev, c, kernel, padding=kernel // 2),
                       nn.BatchNorm1d(c), nn.ReLU(),
                       nn.MaxPool1d(2)]
            prev = c
        self.features = nn.Sequential(*blocks)
        self.proj = nn.Linear(prev, out_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, seq: torch.Tensor) -> torch.Tensor:
        # seq: (b, L) int -> (b, embed, L)
        h = self.embed(seq).transpose(1, 2)
        h = self.features(h).mean(dim=-1)
        return torch.relu(self.proj(self.dropout(h)))
