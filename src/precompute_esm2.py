#!/usr/bin/env python3
"""PREREG_DAVIS_R2 step: precompute FROZEN ESM-2 (esm2_t6_8M_UR50D) per-residue
embeddings for all targets in data_cache/davis_proteins_corrected.json.

Output: data_cache/esm2/<safe_target_name>.npy  (fp16, shape [L, 320])
Skips targets already cached (resume-safe). No gradient, eval mode.
Usage: python3 src/precompute_esm2.py [--limit N]
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data_cache/davis_proteins_corrected.json"
OUT = ROOT / "data_cache/esm2"; OUT.mkdir(exist_ok=True)
MAXLEN = 1022  # ESM-2 t6 context 1024 incl BOS/EOS; longer sequences truncated (recorded)

def safe(name):
    return name.replace("(", "_").replace(")", "").replace("/", "_")

def main():
    import torch
    torch.set_num_threads(2)
    limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 10**9
    prots = json.loads(SRC.read_text())
    t0 = time.time()
    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t6_8M_UR50D", verbose=False)
    model.eval()
    bc = alphabet.get_batch_converter()
    done = 0
    trunc_log = {}
    for name, seq in prots.items():
        out = OUT / f"{safe(name)}.npy"
        if out.exists():
            continue
        s = seq[:MAXLEN]
        if len(seq) > MAXLEN:
            trunc_log[name] = len(seq)
        _, _, toks = bc([(name, s)])
        with torch.no_grad():
            r = model(toks, repr_layers=[6], return_contacts=False)
        emb = r["representations"][6][0, 1:len(s) + 1].to(torch.float16).numpy()
        np.save(out, emb)
        done += 1
        if done % 10 == 0:
            print(f"{done} embedded this call ({time.time()-t0:.0f}s)", flush=True)
        if done >= limit:
            break
    if trunc_log:
        log = OUT / "_truncated.json"
        old = json.loads(log.read_text()) if log.exists() else {}
        old.update(trunc_log); log.write_text(json.dumps(old, indent=1))
    print(f"DONE this call: {done} embedded, {len(list(OUT.glob('*.npy')))} cached total", flush=True)

if __name__ == "__main__":
    main()
