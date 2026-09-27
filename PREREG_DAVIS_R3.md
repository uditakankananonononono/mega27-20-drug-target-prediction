# PREREG_DAVIS_R3 — locked R3 rung (DAVIS DTI), 2026-09-27

Locked BEFORE any R3 training. Ordered via parent confirmation 1:36 PM IST,
after R2's locked-gate falsification (@959493a: val_ci 0.7379 < 0.80 at ep40)
and the Gemini R3 consult (archived isef_judge/GEMINI_R3_CONSULT.txt,
supplementary history, key label-transform lead VERIFIED NON-APPLICABLE:
labels already pKd = -log10(Kd/1e9) from the hkmztrk/DeepDTA mirror,
setting1 splits; R2 MSE 0.873 is a true architecture gap in correct units).

## Carried-forward locks (verbatim from R2)
- Labels: pKd = -log10(Kd/1e9), missing Kd=0 -> 1e5 nM (DeepDTA convention).
- Splits: DeepDTA setting1 train/test fold files, untouched.
- Internal val: 1200 train pairs, rng seed 0 (same split as R2, reused for
  early stop so R2/R3 val numbers are directly comparable).
- Primary benchmark: CI > 0.878 AND MSE < 0.261 vs DeepDTA on the benchmark
  test set. Any beat claim requires BOTH on the locked test set.

## R3 design (frozen)
PRIMARY = D1+D2 hybrid:
1. Protein: crop each of the 442 kinase sequences to its Pfam PF00069
   (Pkinase) domain boundaries. Crop table built ONCE from sequence-only
   evidence (HMM/profile alignment or annotated domain map), saved to the
   repo, and frozen before training; kinases without a resolvable domain
   keep full length and are counted + declared.
2. Encoder: ESM-2 t6_8M with the TOP 2 transformer layers unfrozen
   (lr 1e-5 on unfrozen layers, 1e-3 on heads), attention pooling over the
   cropped sequence; drug tower = R1 graph tower unchanged.
3. Training: 8k-pair chunks, batch 256, resume-safe checkpoints, locked
   EARLY STOP: stop when val_ci fails to improve by >= 0.002 over 8
   consecutive chunks; the locked test evaluation runs ONCE on the
   early-stopped checkpoint.
4. Falsification gate: if val_ci < 0.80 at the earlier of ep40 or early stop,
   R3 is declared failed with full numbers; the negative goes in the paper;
   no post-hoc architecture tweaks within R3 (a tweak = R4, new prereg).

ALTERNATIVE (recorded fallback, not executed in R3): DeepDTA-style end-to-end
1D-CNNs on both modalities with the same label/split/val locks. Executed only
by a future owner/parent decision as its own locked rung.

## Compute honesty
Same 2-core box; ESM-2 8M forward+backward on cropped domains (~250-350 aa)
keeps memory inside the 1GB budget (full-length unfrozen would not - the
crop is also a compute necessity, stated plainly). Unfrozen-layer gradients
are checkpointed with the same resume discipline as R2.
