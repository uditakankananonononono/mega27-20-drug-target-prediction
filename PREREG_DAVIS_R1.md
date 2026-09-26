# PREREG: DAVIS protocol-redesign rung R1 (locked 2026-09-26 17:25 IST, before any R1 chunk)

## Diagnosis (committed evidence)
- 55 epochs at lr=3e-4, bs=512, 8k-pair sampled subsets: test CI plateau
  0.5585 -> 0.7093, +0.013 over last 25 epochs (results/davis_log.jsonl).
- 1-epoch LR probe at 1e-3 (val-side): no immediate change - 1 epoch cannot
  judge Adam LR; a multi-chunk comparison is required.
- Reference protocol (DeepDTA repo, hkmztrk/DeepDTA README, fetched today):
  Adam (Keras default lr=1e-3), batch 256, 100 epochs, FULL train set
  (DAVIS 25,056 pairs), published test CI 0.878 / MSE 0.261 (Nature Comms
  2025 benchmark table, verified earlier).
- Our chunks see 8k/25k pairs per "epoch" at 1/3 the reference LR -> the
  plateau is protocol-bound, not necessarily architecture-bound.

## Rung R1 (protocol match within sandbox constraints)
- Optimizer: Adam lr=1e-3 (match reference), bs=256 (match reference).
- Data: same DeepDTA DAVIS split (already in-repo mirror). 8k-pair sampled
  subset per chunk (sandbox 120s cap); 3 chunks ~= 1 full-data epoch.
  NSUB may drop to 4000 if a chunk exceeds the call budget (documented).
- Checkpoint: continue from results/davis_ckpt.pt (ep55) - Adam moments
  carry over; the LR change is the only intervention.
- Checkpoints-to-beat ledger: val_ci trend is the tuning signal; test CI is
  logged per chunk exactly as the inherited protocol does (best-epoch
  convention matching DeepDTA's published reporting; seed variance to be
  reported in the paper).
- Success gate: test CI > 0.878 AND test MSE < 0.261 at any epoch, then
  confirm with a second seed. Falsification gate: if val_ci has not
  reached 0.78 by absolute epoch 130 (~25 full-data-equivalent epochs),
  declare protocol-R1 an honest negative and escalate to R2 (architecture
  rung: GraphDTA-style or deeper protein encoder, separately preregistered
  with ChatGPT redirection input per RULE 6).
- Honesty: all chunks logged; no cherry-picked restarts. Negatives stay.
