"""R3 (PREREG_DAVIS_R3.md, locked 2026-09-27): D1+D2 hybrid.
Protein = Pfam PF00069 domain crop (data_cache/davis_pfam_crop.json, frozen;
full-length fallbacks counted + declared) -> ESM-2 t6_8M with TOP 2 transformer
layers unfrozen (lr 1e-5 on those, 1e-3 on heads) -> attention pooling ->
R1 drug tower (MolGCN out_dim=32, unchanged) + R1 head.
Locked protocol: 8k-pair chunks, batch 256, resume-safe ckpt, affinity-weighted
MSE, internal val = same 1200-pair rng-seed-0 split as R2.
Locked early stop: val_ci fails to improve by >= 0.002 over 8 consecutive
chunks. Locked falsification gate: val_ci < 0.80 at the earlier of chunk 40 or
early stop -> FAILED. Locked test evaluation runs ONCE on the stopped ckpt.
Primary gate (locked): test CI > 0.878 AND test MSE < 0.261 (DeepDTA).

Usage: python3 run_davis_r3.py [NCHUNKS] [NSUB]   (defaults 1, 8000)
"""
import json, os, sys, time
sys.path.insert(0, "src")
import numpy as np
import torch
import torch.nn as nn
from targetscan.data.davis_kiba import load_davis
from targetscan.data.featurize import N_ATOM_FEATS, smiles_to_graph
from targetscan.models.mol_gnn import MolGCN, collate_graphs
from targetscan.metrics import concordance_index_fast, mse

NCHUNKS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
NSUB = int(sys.argv[2]) if len(sys.argv) > 2 else 8000
CKPT = "results/davis_r3_ckpt.pt"
LOG = "results/davis_r3_log.jsonl"
VERDICT = "results/davis_r3_verdict.json"
MAX_CHUNKS = 40
PATIENCE, MIN_DELTA = 8, 0.002


class AttnPool(nn.Module):
    def __init__(self, d_in=320, d_attn=64):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(d_in, d_attn), nn.Tanh())
        self.score = nn.Linear(d_attn, 1)

    def forward(self, h, mask):  # h [B,L,320], mask [B,L] bool (True=pad)
        s = self.score(self.proj(h)).squeeze(-1)
        s = s.masked_fill(mask, -1e9)
        w = torch.softmax(s, dim=1)
        return (h * w.unsqueeze(-1)).sum(1)


class R3Net(nn.Module):
    def __init__(self, n_atom_feats, esm, alphabet, out_dim=32, dropout=0.2):
        super().__init__()
        self.esm = esm
        self.alphabet = alphabet
        for p in self.esm.parameters():
            p.requires_grad = False
        for layer in self.esm.layers[-2:]:  # top 2 transformer layers unfrozen
            for p in layer.parameters():
                p.requires_grad = True
        self.mol = MolGCN(n_atom_feats, out_dim=out_dim, dropout=dropout)
        self.pool = AttnPool(320)
        self.prot_head = nn.Sequential(
            nn.Linear(320, out_dim), nn.ReLU(), nn.Dropout(dropout))
        self.head = nn.Sequential(
            nn.Linear(2 * out_dim, 256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(128, 1))

    def esm_forward(self, tokens):
        """Frozen stage (embed + layers[:-2]) under no_grad, grad stage after."""
        padding_mask = tokens == self.alphabet.padding_idx
        with torch.no_grad():
            x = self.esm.embed_scale * self.esm.embed_tokens(tokens)
            x = x.transpose(0, 1)  # [L, B, E] for the transformer layers
            for layer in self.esm.layers[:-2]:
                x = layer(x, self_attn_padding_mask=padding_mask)[0]
        x = x.detach().requires_grad_(True)
        for layer in self.esm.layers[-2:]:
            x = layer(x, self_attn_padding_mask=padding_mask)[0]
        x = self.esm.emb_layer_norm_after(x)
        x = x.transpose(0, 1)  # back to [B, L, E]
        return x[:, 1:-1], padding_mask[:, 1:-1]  # strip BOS/EOS

    def forward(self, gx, a_norm, gmask, tokens):
        h, pad = self.esm_forward(tokens)
        z = torch.cat([self.mol(gx, a_norm, gmask),
                       self.prot_head(self.pool(h, pad))], dim=1)
        return self.head(z).squeeze(-1)


class R3Store:
    def __init__(self, smiles, cropped_seqs):
        self.graphs = [smiles_to_graph(s) for s in smiles]
        self.seqs = cropped_seqs

    def batch(self, pairs, bc):
        d_idx, t_idx = pairs[:, 0], pairs[:, 1]
        x, a_norm, mask = collate_graphs([self.graphs[i] for i in d_idx])
        _, _, toks = bc([(str(i), self.seqs[i]) for i in t_idx])
        return x, a_norm, mask, toks


@torch.no_grad()
def predict(model, store, pairs, bc, bs=256):
    model.eval()
    seqlen = np.array([len(s) for s in store.seqs])
    out = []
    for i in range(0, len(pairs), bs):
        chunk_idx = np.arange(i, min(i + bs, len(pairs)))
        chunk_idx = chunk_idx[np.argsort(seqlen[pairs[chunk_idx][:, 1]])]
        res = {}
        j = 0
        while j < len(chunk_idx):
            L = int(seqlen[pairs[chunk_idx[j], 1]]) + 2
            mb = max(1, min(8, int(1_200_000 / (L * L))))
            sub = chunk_idx[j:j + mb]; j += mb
            x, a_norm, mask, toks = store.batch(pairs[sub], bc)
            for k, p_ in zip(sub, model(x, a_norm, mask, toks).numpy()):
                res[k] = p_
        out.append(np.array([res[k] for k in range(i, i + len(chunk_idx))]))
    return np.concatenate(out)


def evaluate(model, store, pairs, y, bc, bs=256):
    p = predict(model, store, pairs, bc, bs)
    return {"mse": mse(y, p), "ci": concordance_index_fast(np.asarray(y), p),
            "n": int(len(y))}


def main():
    torch.set_num_threads(2)
    ds = load_davis()
    prot_names = list(json.loads(
        open("data_cache/davis_proteins_corrected.json").read()).keys())
    assert len(prot_names) == len(ds.sequences)
    seqs = json.loads(open("data_cache/davis_proteins_corrected.json").read())
    crop = json.loads(open("data_cache/davis_pfam_crop.json").read())
    # ESM-2 hard context = 1024 tokens (1022 residues + BOS/EOS); full-length
    # fallbacks longer than 1022 are capped to 1022 (same recorded convention
    # as src/precompute_esm2.py MAXLEN) - a model hard limit, not an
    # architecture change; count logged for the paper.
    cropped, n_full, n_trunc = [], 0, 0
    for n_ in prot_names:
        c = crop[n_]
        if c["source"].startswith("full_length"):
            n_full += 1
        s = seqs[n_][c["start"]:c["end"]]
        if len(s) > 1022:
            s = s[:1022]; n_trunc += 1
        cropped.append(s)
    print(f"crop: {len(cropped) - n_full} domain-cropped, {n_full} full-length "
          f"fallbacks ({n_trunc} context-capped at 1022) - declared", flush=True)

    esm, alphabet = torch.hub.load("facebookresearch/esm:main",
                                   "esm2_t6_8M_UR50D", verbose=False)
    esm.eval()
    bc = alphabet.get_batch_converter()
    net = R3Net(N_ATOM_FEATS, esm, alphabet, out_dim=32)
    store = R3Store(ds.smiles, cropped)

    rng = np.random.default_rng(0)  # identical val split to R2
    n = len(ds.train_pairs)
    val_idx = rng.permutation(n)[:1200]
    tr_mask = np.ones(n, bool); tr_mask[val_idx] = False
    tr_all, val_pairs = ds.train_pairs[tr_mask], ds.train_pairs[val_idx]
    y_all, y_val = ds.y_train[tr_mask], ds.y_train[val_idx]
    yt_val = None

    unfrozen = [p for p in net.esm.parameters() if p.requires_grad]
    heads = ([p for p in net.mol.parameters()] + [p for p in net.pool.parameters()]
             + [p for p in net.prot_head.parameters()]
             + [p for p in net.head.parameters()])
    opt = torch.optim.Adam([
        {"params": unfrozen, "lr": 1e-5},
        {"params": heads, "lr": 1e-3}], weight_decay=1e-5)

    start_chunk, start_batch, best_val, stall, done = 0, 0, -1.0, 0, False
    torch_rng, np_rng = None, None
    if os.path.exists(CKPT):
        state = torch.load(CKPT, weights_only=False)
        net.load_state_dict(state["model"])
        opt.load_state_dict(state["opt"])
        start_chunk = state["chunk"]
        start_batch = state.get("batch_in_chunk", 0)
        best_val, stall, done = state["best_val"], state["stall"], state["done"]
        torch_rng, np_rng = state.get("torch_rng"), state.get("np_rng")
        print(f"resumed chunk {start_chunk} batch {start_batch} "
              f"best_val {best_val:.4f} stall {stall} done {done}", flush=True)
    if done:
        print("R3 already stopped; nothing to do", flush=True)
        return

    yvt = torch.tensor(np.asarray(y_val), dtype=torch.float32)
    for chunk in range(start_chunk, min(start_chunk + NCHUNKS, MAX_CHUNKS)):
        sub = np.random.default_rng(chunk + 1).permutation(len(tr_all))[:NSUB]
        pairs, ys = tr_all[sub], y_all[sub]
        yt = torch.tensor(np.asarray(ys), dtype=torch.float32)
        net.train()
        np.random.seed(chunk + 1); torch.manual_seed(chunk + 1)
        perm = np.random.permutation(len(pairs))
        if chunk == start_chunk and start_batch and torch_rng is not None:
            torch.set_rng_state(torch_rng)
            np.random.set_state(np_rng)
        t0 = time.time(); tot = 0.0
        seqlen = np.array([len(s) for s in store.seqs])
        b0 = start_batch if chunk == start_chunk else 0
        for i in range(b0, len(pairs), 256):
            idx = perm[i:i + 256]
            # locked 256-pair optimizer step; attention memory forces
            # length-sorted adaptive micro-batches (implementation detail,
            # optimizer semantics identical: weighted MSE mean over the 256)
            idx = idx[np.argsort(seqlen[pairs[idx][:, 1]])]
            opt.zero_grad()
            j = 0
            while j < len(idx):
                L = int(seqlen[pairs[idx[j], 1]]) + 2
                mb = max(1, min(8, int(1_200_000 / (L * L))))
                sub_idx = idx[j:j + mb]; j += mb
                x, a_norm, mask, toks = store.batch(pairs[sub_idx], bc)
                pred = net(x, a_norm, mask, toks)
                w = torch.where(yt[sub_idx] >= 7.0,
                                torch.full_like(yt[sub_idx], 6.0),
                                torch.ones_like(yt[sub_idx]))
                loss = (w * (pred - yt[sub_idx]) ** 2).sum() / len(idx)
                loss.backward()
                tot += loss.item() * len(idx)
            opt.step()
            torch.save({"model": net.state_dict(), "opt": opt.state_dict(),
                        "chunk": chunk, "batch_in_chunk": i + 256,
                        "best_val": best_val, "stall": stall, "done": False,
                        "torch_rng": torch.get_rng_state(),
                        "np_rng": np.random.get_state()}, CKPT)
        vm = evaluate(net, store, val_pairs, y_val, bc)
        improved = vm["ci"] >= best_val + MIN_DELTA
        best_val = max(best_val, vm["ci"])
        stall = 0 if improved else stall + 1
        rec = {"chunk": chunk + 1, "train_mse": round(tot / len(pairs), 4),
               "val_mse": round(vm["mse"], 4), "val_ci": round(vm["ci"], 4),
               "best_val_ci": round(best_val, 4), "stall": stall,
               "wall_s": round(time.time() - t0, 1)}
        with open(LOG, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print("CHUNK", rec, flush=True)
        stop = (stall >= PATIENCE) or (chunk + 1 >= MAX_CHUNKS)
        torch.save({"model": net.state_dict(), "opt": opt.state_dict(),
                    "chunk": chunk + 1, "batch_in_chunk": 0,
                    "best_val": best_val, "stall": stall, "done": stop}, CKPT)
        if stop:
            tm = evaluate(net, store, ds.test_pairs, ds.y_test, bc)  # ONCE
            verdict = {
                "stopped_at_chunk": chunk + 1,
                "stop_reason": "early_stop" if stall >= PATIENCE else "chunk_cap_40",
                "final_val_ci": vm["ci"], "best_val_ci": best_val,
                "falsification_gate_val_ci_ge_0.80": bool(vm["ci"] >= 0.80),
                "test_ci": tm["ci"], "test_mse": tm["mse"],
                "primary_gate": {"test_ci_gt_0.878": bool(tm["ci"] > 0.878),
                                 "test_mse_lt_0.261": bool(tm["mse"] < 0.261),
                                 "MET": bool(tm["ci"] > 0.878 and tm["mse"] < 0.261)},
                "n_full_length_fallbacks": n_full,
                "labels_splits_val_locks": "verbatim R2 carries (pKd, setting1, 1200/rng0)"}
            json.dump(verdict, open(VERDICT, "w"), indent=1)
            print("R3 STOPPED", json.dumps(verdict), flush=True)
            return


if __name__ == "__main__":
    main()
