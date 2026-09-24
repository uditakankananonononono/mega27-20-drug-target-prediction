"""Chunked DAVIS training. Each call: resume ckpt, train EPOCHS passes over a
random 8k-pair subset (SGD over the full 25k across calls), eval test, save."""
import json, os, sys, time
sys.path.insert(0, "src")
import numpy as np
import torch
from targetscan.data.davis_kiba import load_davis
from targetscan.data.featurize import N_ATOM_FEATS
from targetscan.models.dti import DTINet
from targetscan.train import FeatureStore, evaluate, train_dti

EPOCHS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CKPT = sys.argv[2] if len(sys.argv) > 2 else "results/davis_ckpt.pt"
PROT = sys.argv[3] if len(sys.argv) > 3 else "cnn"
LOG = sys.argv[4] if len(sys.argv) > 4 else "results/davis_log.jsonl"
NSUB = int(sys.argv[5]) if len(sys.argv) > 5 else 8000  # 0 = full training set
ds = load_davis()
store = FeatureStore(ds.smiles, ds.sequences, max_len=600)
rng = np.random.default_rng(0)
n = len(ds.train_pairs)
val_idx = rng.permutation(n)[:1200]
tr_mask = np.ones(n, bool); tr_mask[val_idx] = False
tr_all, val_pairs = ds.train_pairs[tr_mask], ds.train_pairs[val_idx]
y = ds.y_train
y_all, y_val = y[tr_mask], y[val_idx]

ckpt = CKPT
net = DTINet(N_ATOM_FEATS, out_dim=32, prot_encoder=PROT)
start_ep = 0
opt_state = None
if os.path.exists(ckpt):
    state = torch.load(ckpt, weights_only=False)
    net.load_state_dict(state["model"])
    opt_state = state.get("opt")
    start_ep = state["epoch"]
    print("resumed at epoch", start_ep, flush=True)

if NSUB > 0:
    sub = np.random.default_rng(start_ep + 1).permutation(len(tr_all))[:NSUB]
elif NSUB < 0:  # fixed subset, same pairs every chunk (no sampling noise)
    sub = np.random.default_rng(0).permutation(len(tr_all))[:-NSUB]
else:
    sub = np.arange(len(tr_all))
t0 = time.time()
net, opt = train_dti(net, store, tr_all[sub], y_all[sub], val_pairs, y_val,
                     epochs=EPOCHS, bs=512, lr=3e-4, patience=99,
                     seed=start_ep + 1, log=lambda *a: print(*a, flush=True),
                     opt_state=opt_state)
abs_ep = start_ep + EPOCHS
torch.save({"model": net.state_dict(), "epoch": abs_ep, "opt": opt.state_dict()}, ckpt)
m = evaluate(net, store, ds.test_pairs, ds.y_test)
rec = {"epoch": abs_ep, "test_mse": m["mse"], "test_ci": m["ci"],
       "wall_s": round(time.time() - t0, 1)}
with open(LOG, "a") as fh:
    fh.write(json.dumps(rec) + "\n")
print("CHUNK DONE", rec, flush=True)
