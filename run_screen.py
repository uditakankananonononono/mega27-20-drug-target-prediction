"""Kinome repurposing screen: score all DAVIS drug-kinase pairs, flag pairs
where the model predicts strong binding (pKd >= THRESH) but DAVIS itself
records weak/no binding (Kd >= 10 uM or missing). Candidates go to ChEMBL
corroboration - the independent evidence layer."""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import torch
from targetscan.data.davis_kiba import load_davis, _fetch, _load_matrix
from targetscan.data.featurize import N_ATOM_FEATS
from targetscan.models.dti import DTINet
from targetscan.train import FeatureStore, predict

THRESH = 7.0  # pKd >= 7 ~ Kd <= 100 nM predicted
ds = load_davis()
kd_raw = _load_matrix(_fetch("davis/Y", binary=True)).astype(np.float64)
store = FeatureStore(ds.smiles, ds.sequences, max_len=600)
net = DTINet(N_ATOM_FEATS, out_dim=32)
state = torch.load("results/davis_ckpt.pt", weights_only=False)
net.load_state_dict(state["model"])
print("model epoch:", state["epoch"], flush=True)

n_d, n_t = ds.affinity.shape
pairs = np.array([[i, j] for i in range(n_d) for j in range(n_t)], dtype=np.int64)
pred = predict(net, store, pairs, bs=512).reshape(n_d, n_t)
np.save("results/davis_full_pred.npy", pred)

cands = []
for i in range(n_d):
    for j in range(n_t):
        if pred[i, j] >= THRESH and (kd_raw[i, j] == 0 or kd_raw[i, j] >= 10000):
            cands.append({"drug_idx": i, "target_idx": j,
                          "pred_pkd": float(pred[i, j]),
                          "davis_kd_nm": float(kd_raw[i, j])})
cands.sort(key=lambda c: -c["pred_pkd"])
json.dump(cands, open("results/repurposing_candidates.json", "w"), indent=2)
print(f"candidates (pred pKd>={THRESH}, DAVIS weak/untested): {len(cands)}", flush=True)
for c in cands[:15]:
    print(f"  drug {c['drug_idx']} x target {c['target_idx']}: pred pKd {c['pred_pkd']:.2f}, DAVIS Kd {c['davis_kd_nm']:.0f} nM", flush=True)
