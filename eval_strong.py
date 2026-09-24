"""Strong-binder (pKd>=7) eval for a trained DAVIS checkpoint."""
import sys
sys.path.insert(0, "src")
import numpy as np
import torch
from targetscan.data.davis_kiba import load_davis
from targetscan.data.featurize import N_ATOM_FEATS
from targetscan.models.dti import DTINet
from targetscan.train import FeatureStore, predict

ckpt = sys.argv[1] if len(sys.argv) > 1 else "results/davis_ckpt.pt"
prot = sys.argv[2] if len(sys.argv) > 2 else "cnn"
ds = load_davis()
store = FeatureStore(ds.smiles, ds.sequences, max_len=600)
net = DTINet(N_ATOM_FEATS, out_dim=32, prot_encoder=prot)
state = torch.load(ckpt, weights_only=False)
net.load_state_dict(state["model"])
net.eval()
pred = predict(net, store, ds.test_pairs)
true = ds.y_test
m = true >= 7.0
import json
rec = {"ckpt": ckpt, "epoch": state["epoch"], "n_strong": int(m.sum()),
                  "strong_pred_mean": round(float(pred[m].mean()), 3),
                  "strong_true_mean": round(float(true[m].mean()), 3),
                  "strong_mse": round(float(((pred[m]-true[m])**2).mean()), 3),
                  "all_pred_mean": round(float(pred.mean()), 3),
                  "all_pred_std": round(float(pred.std()), 3)}
print(json.dumps(rec))
with open("results/strong_eval_log.jsonl", "a") as fh:
    fh.write(json.dumps(rec) + "\n")
