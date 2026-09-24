"""targetscan CLI - usable tool built on the item-20 work.

  python -m targetscan screen --smiles "CCO" [--top 10]   score a drug vs 442 DAVIS kinases
  python -m targetscan eval                               verified benchmark numbers
  python -m targetscan verify --inchikey X --gene PLK1    ChEMBL corroboration
  python -m targetscan audit [--threshold 1.0]            dataset consistency audit
"""
import argparse, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import torch


def _load(ckpt):
    from targetscan.data.davis_kiba import load_davis, _fetch
    from targetscan.data.featurize import N_ATOM_FEATS
    from targetscan.models.dti import DTINet
    from targetscan.train import FeatureStore
    ds = load_davis()
    names = list(json.loads(_fetch("davis/proteins.txt")).keys())
    store = FeatureStore(ds.smiles, ds.sequences, max_len=600)
    net = DTINet(N_ATOM_FEATS, out_dim=32)
    state = torch.load(ckpt, weights_only=False)
    net.load_state_dict(state["model"])
    net.eval()
    return ds, names, store, net, state["epoch"]


def cmd_screen(a):
    ds, names, store, net, ep = _load(a.ckpt)
    from targetscan.data.featurize import smiles_to_graph
    from targetscan.models.mol_gnn import collate_graphs
    g = smiles_to_graph(a.smiles)
    preds = []
    B = 256
    with torch.no_grad():
        for i in range(0, len(ds.sequences), B):
            n = min(B, len(ds.sequences) - i)
            x, a_norm, mask = collate_graphs([g] * n)
            seq = torch.stack([store.seqs[j] for j in range(i, i + n)])
            preds.append(net(x, a_norm, mask, seq).numpy())
    p = np.concatenate(preds)
    order = np.argsort(-p)[: a.top]
    print(json.dumps({"smiles": a.smiles, "model_epoch": ep, "model_note": "pKd predicted",
                      "top_targets": [{"target": names[i], "pred_pKd": round(float(p[i]), 3)}
                                      for i in order]}, indent=1))


def cmd_eval(a):
    log = os.path.join(os.path.dirname(a.ckpt), "davis_log.jsonl")
    last = json.loads(open(log).readlines()[-1])
    print(json.dumps({"benchmark": "DAVIS", "epoch": last["epoch"],
                      "test_ci": round(last["test_ci"], 4),
                      "test_mse": round(last["test_mse"], 4),
                      "source_file": log,
                      "reference_DeepDTA_ci": 0.878}, indent=1))


def cmd_audit(a):
    from targetscan.audit import audit
    from targetscan.data.davis_kiba import load_davis, _fetch
    ds = load_davis()
    raw = json.loads(_fetch("davis/proteins.txt"))
    rep = audit(raw, ds.affinity, list(raw.keys()), threshold=a.threshold)
    rep["dataset"] = "davis (DeepDTA mirror)"
    print(json.dumps(rep, indent=1, default=str))


def cmd_verify(a):
    from targetscan import chembl
    res = chembl.corroborate(a.inchikey, a.gene)
    print(json.dumps(res, indent=1, default=str))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="targetscan")
    ap.add_argument("--ckpt", default="results/davis_ckpt.pt")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen"); s.add_argument("--smiles", required=True)
    s.add_argument("--top", type=int, default=10)
    sub.add_parser("eval")
    v = sub.add_parser("verify"); v.add_argument("--inchikey", required=True)
    v.add_argument("--gene", required=True)
    au = sub.add_parser("audit"); au.add_argument("--threshold", type=float, default=1.0)
    a = ap.parse_args(argv)
    {"screen": cmd_screen, "eval": cmd_eval, "verify": cmd_verify, "audit": cmd_audit}[a.cmd](a)


if __name__ == "__main__":
    main()
