# DAVIS R1 VERDICT (declared 2026-09-27 03:06 IST, per locked PREREG_DAVIS_R1.md)

## Falsification gate outcome: TRIGGERED
- Locked gate: val_ci must reach 0.78 by absolute epoch 130.
- Observed val_ci at ep130: **0.7610** (< 0.78). Epochs 126-130 val_ci:
  ep126 (chunk val not printed per-chunk in this log; per-epoch val_ci from
  trainer stdout): 126 -, 127 0.7565, 128 0.7655, 129 0.7568, 130 0.7610.
- Test CI (best-epoch convention, matching DeepDTA reporting): best 0.7533
  @ ep128; final ep130 0.7496. Target (DeepDTA published): CI 0.878 / MSE
  0.261. NOT MET - honest shortfall, ~0.125 CI below benchmark.

## Declaration
Protocol-R1 is an HONEST NEGATIVE: matching the reference optimizer (Adam
lr=1e-3) and batch size (256) on 8k-pair sampled chunks, continued from the
ep55 checkpoint, lifted test CI from the 0.7093 plateau to 0.7533 (+0.044)
but did not reach the benchmark and missed the locked falsification gate.
All 130 chunks logged in results/davis_log.jsonl; no cherry-picked restarts.

## Escalation (locked in prereg)
Escalate to rung R2: architecture change (GraphDTA-style or deeper protein
encoder), SEPARATELY PREREGISTERED with ChatGPT redirection input per RULE 6
(stuck negative -> judge redirection, pivot on strongest signal, verbatim
log) BEFORE any R2 training. Lane-20 R1 judge consult is staged
(/tmp/gem_l20_01-04.txt + /tmp/gem_l20_q.txt) and will be extended with
this falsification outcome; it queues behind the lane-16 R1 consult on the
next browser token.

## What R1 still delivered
- Diagnosis confirmed: the 0.7093 plateau was substantially protocol-bound
  (LR 3x too low, 1/3 batch-visit rate): +0.044 CI from protocol match alone.
- This negative + the protocol diagnosis fold into the paper's limitations
  and methods-audit sections (compact treatment per user rule 8:27:49).
