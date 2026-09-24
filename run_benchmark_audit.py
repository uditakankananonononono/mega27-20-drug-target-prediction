"""DAVIS benchmark audit: mutation-erasure defect in the DeepDTA-mirror
sequence distribution. Falsifiable claim + evidence pack.

CLAIM: in the DAVIS files shipped by the standard DeepDTA GitHub mirror
(hkmztrk/DeepDTA), every mutant kinase entry's protein sequence is
byte-identical to its wild type (verified by MD5 over the JSON values),
while the affinity matrix carries WT-vs-mutant differences up to 4.49 pKd
for the same drug. Any sequence-based model trained on this distribution
is structurally incapable of representing the mutation, so mutant-vs-WT
affinity differences are unlearnable label noise from its perspective.
"""
import hashlib, json, sys
sys.path.insert(0, "src")
import numpy as np
from targetscan.data.davis_kiba import load_davis

ds = load_davis()
raw = json.load(open("data_cache/davis_proteins.txt"))
names = list(raw.keys())
Y = ds.affinity

muts = [n for n in names if "(" in n and n.split("(")[0] in raw]
def md5(s):
    return hashlib.md5(s.encode()).hexdigest()

per_target, tot_pairs, big_pairs = [], 0, 0
for n in muts:
    wt = n.split("(")[0]
    i, j = names.index(wt), names.index(n)
    identical = raw[n] == raw[wt]
    d = np.abs(Y[:, i] - Y[:, j])
    tot_pairs += Y.shape[0]
    big = int((d >= 1.0).sum())
    big_pairs += big
    per_target.append({"mutant": n, "wt": wt, "seq_byte_identical": identical,
                       "seq_md5_mutant": md5(raw[n]), "seq_md5_wt": md5(raw[wt]),
                       "n_drugs": Y.shape[0], "mean_abs_delta": round(float(d.mean()), 4),
                       "n_abs_delta_ge_1": big, "max_abs_delta": round(float(d.max()), 3)})

lig = json.load(open("data_cache/davis_ligands_can.txt"))
drugs = list(lig.keys())
k = drugs.index("5291")  # imatinib, CID verified via PubChem pull (MW 493.6)
headline = []
for wt, mu in [("ABL1", "ABL1(T315I)"), ("KIT", "KIT(D816V)")]:
    i, j = names.index(wt), names.index(mu)
    headline.append({"drug": "imatinib (CID 5291)", "wt": wt, "wt_pkd": round(float(Y[k, i]), 2),
                     "mutant": mu, "mutant_pkd": round(float(Y[k, j]), 2),
                     "delta_pkd": round(float(Y[k, j] - Y[k, i]), 2),
                     "seq_inputs_identical": raw[wt] == raw[mu]})

out = {
  "claim": "DAVIS-as-distributed (DeepDTA mirror): all mutant entries share WT sequences",
  "evidence": {
    "n_targets": len(names), "n_mutant_entries": len(muts),
    "n_mutants_byte_identical_to_wt": sum(1 for t in per_target if t["seq_byte_identical"]),
    "mutant_drug_pairs": tot_pairs,
    "pairs_abs_pkd_delta_ge_1": big_pairs,
    "pct_pairs_ge_1_pkd": round(100.0 * big_pairs / tot_pairs, 2),
    "max_abs_pkd_delta": max(t["max_abs_delta"] for t in per_target),
    "headline_examples": headline},
  "consequences": [
    "WT-vs-mutant affinity shifts (up to 4.49 pKd, e.g. imatinib vs ABL1(T315I) resistance) are unlearnable from the shipped inputs",
    "published DAVIS CI scores (incl. DeepDTA 0.878, ours 0.708) are computed on a matrix where 54/442 targets carry mutation labels with no mutation signal",
    "model mutation-blindness previously attributed to encoder truncation is actually guaranteed by the dataset itself"],
  "fix_direction": "re-source mutant sequences from UniProt canonical + variant features before any mutation-sensitivity claim",
  "per_target": per_target}
json.dump(out, open("results/benchmark_audit_davis_mutants.json", "w"), indent=1)
print(json.dumps(out["evidence"], indent=1)[:600])
