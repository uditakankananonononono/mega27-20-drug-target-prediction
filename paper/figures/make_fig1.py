"""Regenerate figures/fig1_training_curves.pdf from the committed logs."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
def read(name):
    rows = [json.loads(l) for l in open(os.path.join(ROOT, "results", name))]
    return [r["epoch"] for r in rows], [r["test_ci"] for r in rows]

e1, c1 = read("davis_log.jsonl")
e2, c2 = read("davis_log_attn.jsonl")
plt.figure(figsize=(6, 3.6))
plt.plot(e1, c1, marker="o", ms=3, label="CNN encoder (v1)")
plt.plot(e2, c2, marker="s", ms=3, label="attention encoder (v2)")
plt.axhline(0.878, ls="--", lw=1, c="gray", label="DeepDTA reference 0.878")
plt.xlabel("epoch"); plt.ylabel("held-out concordance index")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), "fig1_training_curves.pdf"))
print("wrote fig1:", max(c1), max(c2))
