"""DAVIS and KIBA DTI benchmarks in DeepDTA layout (hkmztrk/DeepDTA mirror).

DAVIS: 68 kinase inhibitors x 442 kinases, Kd (nM); we model pKd =
-log10(Kd/1e9) exactly as DeepDTA/GraphDTA do, so numbers are comparable.
KIBA: 2,068 drugs x 229 targets, KIBA score (higher = tighter).
Network access lives only in fetch_text; parsing is pure and testable.
"""
from __future__ import annotations

import ast
import io
import json
import os
import pickle
import urllib.request
from dataclasses import dataclass

import numpy as np

BASE = "https://raw.githubusercontent.com/hkmztrk/DeepDTA/master/data"
CACHE_DIR = os.environ.get(
    "TARGETSCAN_CACHE",
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data_cache"),
)


class DataUnavailable(RuntimeError):
    pass


def _fetch(rel: str, binary: bool = False):
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, rel.replace("/", "_"))
    if not os.path.exists(path):
        try:
            with urllib.request.urlopen(f"{BASE}/{rel}", timeout=90) as resp:
                data = resp.read()
        except Exception as exc:
            raise DataUnavailable(f"{rel}: {exc}") from exc
        with open(path, "wb") as fh:
            fh.write(data)
    with open(path, "rb") as fh:
        return fh.read() if binary else fh.read().decode("utf-8", "replace")


@dataclass
class DTIBenchmark:
    name: str
    smiles: list          # (n_drugs,) canonical SMILES
    sequences: list       # (n_targets,) amino-acid sequences
    affinity: np.ndarray  # (n_drugs, n_targets) float32, modeled scale
    train_pairs: np.ndarray  # (m, 2) int64 (drug_idx, target_idx)
    test_pairs: np.ndarray
    y_train: np.ndarray
    y_test: np.ndarray
    transform: str


def _load_matrix(raw: bytes) -> np.ndarray:
    return pickle.loads(raw, encoding="latin1")


def _load_fold(rel: str) -> list:
    """DeepDTA fold files are Python literals; train files hold 5 CV folds
    which we flatten into one training index set (GraphDTA protocol)."""
    parsed = ast.literal_eval(_fetch(rel))
    if isinstance(parsed[0], list):
        return [int(i) for fold in parsed for i in fold]
    return [int(i) for i in parsed]


def load_davis() -> DTIBenchmark:
    smiles = list(json.loads(_fetch("davis/ligands_can.txt")).values())
    sequences = list(json.loads(_fetch("davis/proteins.txt")).values())
    Y = _load_matrix(_fetch("davis/Y", binary=True)).astype(np.float64)
    train_idx = _load_fold("davis/folds/train_fold_setting1.txt")
    test_idx = _load_fold("davis/folds/test_fold_setting1.txt")
    kd = Y.copy()
    kd[kd == 0] = 100000.0  # DeepDTA convention: missing -> weakest bin
    pkd = -np.log10(kd / 1e9)
    n_d, n_t = pkd.shape
    flat = pkd.ravel()
    pairs = np.array([[i // n_t, i % n_t] for i in range(n_d * n_t)],
                     dtype=np.int64)
    return DTIBenchmark(
        name="davis", smiles=smiles, sequences=sequences,
        affinity=pkd.astype(np.float32),
        train_pairs=pairs[train_idx], test_pairs=pairs[test_idx],
        y_train=flat[train_idx].astype(np.float32),
        y_test=flat[test_idx].astype(np.float32),
        transform="pKd = -log10(Kd/1e9), missing Kd=0 -> 1e5 nM (DeepDTA)",
    )


def load_kiba() -> DTIBenchmark:
    smiles = list(json.loads(_fetch("kiba/ligands_can.txt")).values())
    sequences = list(json.loads(_fetch("kiba/proteins.txt")).values())
    Y = _load_matrix(_fetch("kiba/Y", binary=True)).astype(np.float64)
    train_idx = _load_fold("kiba/folds/train_fold_setting1.txt")
    test_idx = _load_fold("kiba/folds/test_fold_setting1.txt")
    n_d, n_t = Y.shape
    flat = Y.ravel()
    pairs = np.array([[i // n_t, i % n_t] for i in range(n_d * n_t)],
                     dtype=np.int64)
    return DTIBenchmark(
        name="kiba", smiles=smiles, sequences=sequences,
        affinity=Y.astype(np.float32),
        train_pairs=pairs[train_idx], test_pairs=pairs[test_idx],
        y_train=flat[train_idx].astype(np.float32),
        y_test=flat[test_idx].astype(np.float32),
        transform="KIBA score as published",
    )


LOADERS = {"davis": load_davis, "kiba": load_kiba}
