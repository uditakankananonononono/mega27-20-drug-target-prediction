"""R2 (PREREG_DAVIS_R2.md): frozen ESM-2 (t6_8M, 320-d) per-residue embeddings of the
CORRECTED DAVIS protein sequences (data_cache/esm2/*.npy, provenance
results/uniprot_variant_pull_provenance.json), learned attention-pool + R1 drug tower
(MolGCN, out_dim=32) + R1 head. Same chunk-resume protocol as run_davis_r1.py
(each call: EPOCHS passes over a random 8k-pair train subset, 1200-pair val split,
affinity-weighted MSE, eval test, checkpoint).

Usage: python3 run_davis_r2.py [EPOCHS] [CKPT] [LOG] [NSUB]
Falsification gate (locked): val_ci < 0.80 at ep40 -> declare, escalate to R3.
Primary gate (locked): test CI > 0.878 AND test MSE < 0.261 (DeepDTA)."""
import json, os, sys, time
sys.path.insert(0, "src")
import numpy as np
import torch
import torch.nn as nn
from targetscan.data.davis_kiba import load_davis
from targetscan.data.featurize import N_ATOM_FEATS, smiles_to_graph
from targetscan.models.mol_gnn import MolGCN, collate_graphs
from targetscan.metrics import concordance_index_fast, mse

EPOCHS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CKPT = sys.argv[2] if len(sys.argv) > 2 else "results/davis_r2_ckpt.pt"
LOG = sys.argv[3] if len(sys.argv) > 3 else "results/davis_r2_log.jsonl"
NSUB = int(sys.argv[4]) if len(sys.argv) > 4 else 8000
ESM_DIR = "data_cache/esm2"


def safe(name):
    return name.replace("(", "_").replace(")", "").replace("/", "_")


class AttnPool(nn.Module):
    """Learned attention pooling over frozen per-residue embeddings."""

    def __init__(self, d_in=320, d_attn=64):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(d_in, d_attn), nn.Tanh())
        self.score = nn.Linear(d_attn, 1)

    def forward(self, h, mask):  # h [B,L,320], mask [B,L] bool (True=pad)
        s = self.score(self.proj(h)).squeeze(-1)
        s = s.masked_fill(mask, -1e9)
        w = torch.softmax(s, dim=1)
        return (h * w.unsqueeze(-1)).sum(1)


class ESM2DTINet(nn.Module):
    def __init__(self, n_atom_feats, out_dim=32, dropout=0.2):
        super().__init__()
        self.mol = MolGCN(n_atom_feats, out_dim=out_dim, dropout=dropout)
        self.pool = AttnPool(320)
        self.prot_head = nn.Sequential(
            nn.Linear(320, out_dim), nn.ReLU(), nn.Dropout(dropout))
        self.head = nn.Sequential(
            nn.Linear(2 * out_dim, 256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 1))

    def forward(self, x, a_norm, mask, h, pad):
        z = torch.cat([self.mol(x, a_norm, mask),
                       self.prot_head(self.pool(h, pad))], dim=1)
        return self.head(z).squeeze(-1)


class ESM2Store:
    """Drug graphs + memmapped frozen ESM-2 embeddings, batched on demand."""

    def __init__(self, smiles, prot_names):
        self.graphs = [smiles_to_graph(s) for s in smiles]
        self.paths = [os.path.join(ESM_DIR, safe(n) + ".npy") for n in prot_names]
        for p in self.paths:
            if not os.path.exists(p):
                raise FileNotFoundError(p)
        self._mm = [None] * len(self.paths)

    def emb(self, i):
        if self._mm[i] is None:
            self._mm[i] = np.load(self.paths[i], mmap_mode="r")
        return self._mm[i]

    def batch(self, pairs):
        d_idx, t_idx = pairs[:, 0], pairs[:, 1]
        x, a_norm, mask = collate_graphs([self.graphs[i] for i in d_idx])
        embs = [np.asarray(self.emb(i), dtype=np.float32) for i in t_idx]
        L = max(e.shape[0] for e in embs)
        h = torch.zeros(len(embs), L, 320)
        pad = torch.ones(len(embs), L, dtype=torch.bool)
        for j, e in enumerate(embs):
            h[j, : e.shape[0]] = torch.from_numpy(e)
            pad[j, : e.shape[0]] = False
        return x, a_norm, mask, h, pad


def train(model, store, pairs, y, val_pairs, val_y, epochs, bs, lr, seed,
          opt_state=None, log=print):
    np.random.seed(seed); torch.manual_seed(seed)
    yt = torch.tensor(np.asarray(y), dtype=torch.float32)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    if opt_state is not None:
        opt.load_state_dict(opt_state)
        for g in opt.param_groups:
            g["lr"] = lr
    n = len(pairs)
    for ep in range(epochs):
        model.train()
        perm = np.random.permutation(n)
        tot = 0.0
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            x, a_norm, mask, h, pad = store.batch(pairs[idx])
            pred = model(x, a_norm, mask, h, pad)
            w = torch.where(yt[idx] >= 7.0, torch.full_like(yt[idx], 6.0),
                            torch.ones_like(yt[idx]))
            loss = (w * (pred - yt[idx]) ** 2).mean()
            opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.item() * len(idx)
        vm = evaluate(model, store, val_pairs, val_y, bs)
        log(f"epoch {ep}: train_mse {tot/n:.4f} val_mse {vm['mse']:.4f} "
            f"val_ci {vm['ci']:.4f}", flush=True)
    return model, opt


@torch.no_grad()
def predict(model, store, pairs, bs=256):
    model.eval()
    out = []
    for i in range(0, len(pairs), bs):
        x, a_norm, mask, h, pad = store.batch(pairs[i:i + bs])
        out.append(model(x, a_norm, mask, h, pad).numpy())
    return np.concatenate(out)


def evaluate(model, store, pairs, y, bs=256):
    p = predict(model, store, pairs, bs)
    return {"mse": mse(y, p), "ci": concordance_index_fast(np.asarray(y), p),
            "n": int(len(y))}


def main():
    ds = load_davis()
    prot_names = list(json.loads(
        open("data_cache/davis_proteins_corrected.json").read()).keys())
    assert len(prot_names) == len(ds.sequences), (
        len(prot_names), len(ds.sequences))
    store = ESM2Store(ds.smiles, prot_names)
    rng = np.random.default_rng(0)
    n = len(ds.train_pairs)
    val_idx = rng.permutation(n)[:1200]
    tr_mask = np.ones(n, bool); tr_mask[val_idx] = False
    tr_all, val_pairs = ds.train_pairs[tr_mask], ds.train_pairs[val_idx]
    y_all, y_val = ds.y_train[tr_mask], ds.y_train[val_idx]

    net = ESM2DTINet(N_ATOM_FEATS, out_dim=32)
    start_ep, opt_state = 0, None
    if os.path.exists(CKPT):
        state = torch.load(CKPT, weights_only=False)
        net.load_state_dict(state["model"])
        opt_state = state.get("opt")
        start_ep = state["epoch"]
        print("resumed at epoch", start_ep, flush=True)

    if NSUB > 0:
        sub = np.random.default_rng(start_ep + 1).permutation(len(tr_all))[:NSUB]
    else:
        sub = np.arange(len(tr_all))
    t0 = time.time()
    net, opt = train(net, store, tr_all[sub], y_all[sub], val_pairs, y_val,
                     epochs=EPOCHS, bs=256, lr=1e-3, seed=start_ep + 1,
                     opt_state=opt_state)
    abs_ep = start_ep + EPOCHS
    torch.save({"model": net.state_dict(), "epoch": abs_ep,
                "opt": opt.state_dict()}, CKPT)
    m = evaluate(net, store, ds.test_pairs, ds.y_test)
    rec = {"epoch": abs_ep, "test_mse": m["mse"], "test_ci": m["ci"],
           "val_ci": None, "wall_s": round(time.time() - t0, 1)}
    with open(LOG, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print("CHUNK DONE", rec, flush=True)


if __name__ == "__main__":
    main()
