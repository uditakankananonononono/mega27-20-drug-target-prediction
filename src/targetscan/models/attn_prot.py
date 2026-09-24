"""Attention-pooling protein encoder: a learned query attends over sequence
positions so a single mutated residue can dominate the representation -
the fix for the mean-pooling mutant-blindness found in v1."""
from __future__ import annotations

import torch
import torch.nn as nn


class AttnProtCNN(nn.Module):
    def __init__(self, vocab: int = 21, embed: int = 64,
                 channels: tuple = (64, 96), kernel: int = 8,
                 out_dim: int = 128, dropout: float = 0.2):
        super().__init__()
        self.embed = nn.Embedding(vocab, embed, padding_idx=0)
        blocks, prev = [], embed
        for c in channels:
            blocks += [nn.Conv1d(prev, c, kernel, padding=kernel // 2),
                       nn.BatchNorm1d(c), nn.ReLU()]
            prev = c
        self.features = nn.Sequential(*blocks)
        self.attn_score = nn.Linear(prev, 1)
        self.proj = nn.Linear(prev, out_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, seq: torch.Tensor, return_attn: bool = False):
        # seq: (b, L) int
        pad_mask = (seq > 0).float()                    # (b, L)
        h = self.embed(seq).transpose(1, 2)             # (b, e, L)
        h = self.features(h).transpose(1, 2)[:, :seq.shape[1], :]  # (b, L, c)
        scores = self.attn_score(h).squeeze(-1)         # (b, L)
        scores = scores.masked_fill(pad_mask == 0, float("-inf"))
        alpha = torch.softmax(scores, dim=-1)
        alpha = torch.nan_to_num(alpha, nan=0.0)
        pooled = (h * alpha.unsqueeze(-1)).sum(1)       # (b, c)
        z = torch.relu(self.proj(self.dropout(pooled)))
        if return_attn:
            return z, alpha
        return z
