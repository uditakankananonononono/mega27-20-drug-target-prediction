"""Mutant-sensitivity experiment (the deciding test for the tool claim).

DAVIS carries 7 KIT mutants and 10 EGFR mutants alongside WT, all assayed
against the same 68 drugs. The protein encoder truncates at max_len=600, so
mutations at positions >= 600 are structurally invisible to the model
(predicted affinities must be bit-identical to WT). Mutations at positions
< 600 are visible.

Falsifiable criteria (fixed before running):
  C1 (truncation hypothesis): mutants with first-diff position >= 600 give
      EXACTLY bit-identical predictions to WT for all 68 drugs.
  C2 (sensitivity): mutants with first-diff position < 600 give nonzero
      prediction deltas (mean |delta| > 0 over drugs with measured WT and
      mutant affinities).
  C3 (direction): for visible mutants, Spearman rho between predicted and
      true WT->mutant affinity deltas is reported; rho > 0 = directionally
      useful, rho <= 0 = mutation-blind where it matters.
Verdict: tool mutation-aware iff C2 passes AND median visible rho > 0.
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

MAX_LEN = 600

ds = load_davis()
raw = json.load(open("data_cache/davis_proteins.txt"))
names = list(raw.keys())
assert list(ds.sequences) == list(raw.values())

def variants(wt):
    return [n for n in names if n == wt or (n.startswith(wt + "(") and n.endswith(")"))]

def first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i + 1  # 1-based position
    return None

store = FeatureStore(ds.smiles, ds.sequences, max_len=MAX_LEN)
net = DTINet(N_ATOM_FEATS, out_dim=32, prot_encoder="cnn")
state = torch.load("results/davis_ckpt.pt", weights_only=False)
net.load_state_dict(state["model"])
net.eval()
epoch = state["epoch"]

# predict every drug against every WT/mutant target of interest
fams = {}
for wt in ["KIT", "EGFR"]:
    vs = variants(wt)
    t_idx = [names.index(v) for v in vs]
    pairs = np.array([[d, t] for d in range(len(ds.smiles)) for t in t_idx])
    preds = predict(net, store, pairs).reshape(len(ds.smiles), len(vs))
    Y = ds.affinity[:, t_idx]  # (n_drugs, n_variants), NaN where unmeasured
    wt_col = 0
    muts = []
    for j, v in enumerate(vs[1:], start=1):
        pos = first_diff(raw[wt], raw[v])
        visible = pos is not None and pos < MAX_LEN
        d_pred = preds[:, j] - preds[:, wt_col]
        d_true = Y[:, j] - Y[:, wt_col]
        measured = ~np.isnan(d_true)
        bit_identical = bool(np.all(d_pred == 0.0))
        rho = None
        if visible and measured.sum() >= 10:
            rho = round(float(spearmanr(d_pred[measured], d_true[measured]).statistic), 4)
        muts.append({
            "mutant": v, "first_diff_pos": pos, "visible": visible,
            "bit_identical_preds": bit_identical,
            "mean_abs_pred_delta": round(float(np.abs(d_pred).mean()), 5),
            "n_measured_pairs": int(measured.sum()),
            "spearman_pred_vs_true_delta": rho,
        })
    fams[wt] = {"ckpt_epoch": epoch, "mutants": muts}

c1 = all(m["bit_identical_preds"] for f in fams.values() for m in f["mutants"] if not m["visible"])
c2 = all(m["mean_abs_pred_delta"] > 0 for f in fams.values() for m in f["mutants"] if m["visible"])
rhos = [m["spearman_pred_vs_true_delta"] for f in fams.values() for m in f["mutants"] if m["spearman_pred_vs_true_delta"] is not None]
median_rho = round(float(np.median(rhos)), 4) if rhos else None
verdict = "mutation-aware" if (c2 and median_rho is not None and median_rho > 0) else "mutation-blind"

out = {"protocol": {"max_len": MAX_LEN, "ckpt": "results/davis_ckpt.pt",
                    "ckpt_epoch": epoch, "encoder": "cnn",
                    "criteria": "C1 truncation bit-identity, C2 nonzero deltas for visible, C3 median visible Spearman>0"},
       "families": fams,
       "C1_truncation_bit_identity": c1, "C2_visible_nonzero_deltas": c2,
       "median_visible_spearman": median_rho, "n_visible_with_rho": len(rhos),
       "verdict": verdict}
json.dump(out, open("results/mutant_sensitivity.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ["C1_truncation_bit_identity", "C2_visible_nonzero_deltas", "median_visible_spearman", "verdict"]}, indent=1))
for wt, f in fams.items():
    for m in f["mutants"]:
        print(wt, m["mutant"], "pos", m["first_diff_pos"], "visible" if m["visible"] else "INVISIBLE",
              "bitident" if m["bit_identical_preds"] else "diff", m["spearman_pred_vs_true_delta"])
