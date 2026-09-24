"""Mutant-rescue experiment: does the model separate WT from TRUE mutant
sequences once the dataset erasure is fixed?

Protocol (fixed before running):
  - Rebuild KIT(V559D) and KIT(L576P) sequences by editing the WT string at
    the ClinVar-validated positions (both < 600, inside the encoder window).
  - Predict all 68 DAVIS drugs against WT and each rebuilt mutant with the
    v1 CNN checkpoint.
  - Criteria: R1 mean |pred delta| > 0 (inputs now differ);
    R2 Spearman rho between predicted and true DAVIS deltas (reported; the
    true deltas here are small - max 1.33 pKd - so rho is indicative only).
  - Negative control: rebuild KIT(D816V) at position 816 (> 600 window);
    predictions must stay bit-identical - quantifies the residual
    truncation limit AFTER the dataset fix.
"""
import json, sys
sys.path.insert(0, "src")
import numpy as np
import torch
from scipy.stats import spearmanr
from targetscan.data.davis_kiba import load_davis
from targetscan.data.featurize import N_ATOM_FEATS
from targetscan.models.dti import DTINet
from targetscan.train import FeatureStore, predict

ds = load_davis()
raw = json.load(open("data_cache/davis_proteins.txt"))
names = list(raw.keys())
Y = ds.affinity
wt_seq = raw["KIT"]

def edit(seq, pos1, wt_res, mut_res):
    assert seq[pos1 - 1] == wt_res, f"expected {wt_res} at {pos1}, got {seq[pos1-1]}"
    return seq[:pos1 - 1] + mut_res + seq[pos1:]

variants = {
    "KIT(V559D)_rescued": edit(wt_seq, 559, "V", "D"),
    "KIT(L576P)_rescued": edit(wt_seq, 576, "L", "P"),
    "KIT(D816V)_rescued": edit(wt_seq, 816, "D", "V"),  # negative control: outside window
}

seqs = list(ds.sequences) + list(variants.values())
store = FeatureStore(ds.smiles, seqs, max_len=600)
net = DTINet(N_ATOM_FEATS, out_dim=32, prot_encoder="cnn")
state = torch.load("results/davis_ckpt.pt", weights_only=False)
net.load_state_dict(state["model"])
net.eval()

wt_idx = names.index("KIT")
n_new = len(variants)
t_idx = [wt_idx] + [len(ds.sequences) + i for i in range(n_new)]
pairs = np.array([[d, t] for d in range(len(ds.smiles)) for t in t_idx])
preds = predict(net, store, pairs).reshape(len(ds.smiles), len(t_idx))

res = {"protocol": {"ckpt": "results/davis_ckpt.pt", "ckpt_epoch": state["epoch"],
                    "encoder": "cnn", "max_len": 600}}
for k, (vname, _) in enumerate(variants.items(), start=1):
    d_pred = preds[:, k] - preds[:, 0]
    rec = {"mean_abs_pred_delta": round(float(np.abs(d_pred).mean()), 5),
           "bit_identical_preds": bool(np.all(d_pred == 0.0))}
    davis_name = vname.replace("_rescued", "")
    if davis_name in names:
        d_true = Y[:, names.index(davis_name)] - Y[:, wt_idx]
        rec["true_mean_abs_delta"] = round(float(np.abs(d_true).mean()), 4)
        rec["spearman_pred_vs_true"] = round(float(spearmanr(d_pred, d_true).statistic), 4)
    res[vname] = rec

res["R1_visible_inputs_give_nonzero_deltas"] = all(
    not res[v]["bit_identical_preds"] for v in res if v.endswith("_rescued") and "D816V" not in v)
res["negative_control_D816V_still_invisible"] = res["KIT(D816V)_rescued"]["bit_identical_preds"]
json.dump(res, open("results/mutant_rescue.json", "w"), indent=1)
print(json.dumps(res, indent=1))
