# Result: H4 kinase-family-held-out DAVIS diagnostic (plan frozen and committed first: PLAN_H4_FAMILY_HOLDOUT_2026-10-08.md, commit 6e69b2b). One run, no tuning. Script: src/targetscan/audit_h4.py.
Raw output:
```
mapped 381 families 89 unmapped 61
proteins per fold [77, 76, 76, 76, 76]
M0 global mean pooled MSE 0.8306680271637772 CI 0.4739511760394331 per-fold [0.733, 1.11, 0.928, 0.599, 0.785]
M1 ligand mean pooled MSE 0.6470290956263556 CI 0.7443580083859388 per-fold [0.495, 0.903, 0.757, 0.405, 0.677]
M2 seq-kNN pooled MSE 0.7019509630161755 CI 0.7197606738418375 per-fold [0.627, 0.747, 0.875, 0.517, 0.743]
M2-M1 MSE diff 0.05492186738981997 bootstrap 95% CI [0.02639179635207398, 0.08195891957711526] n scored rows 25908 scored proteins 381
```
Reading (measured only):
- 381 of 442 entries mapped to 89 KinHub families; 61 unmapped stayed in train and were never scored. 25908 rows (381 proteins x 68 ligands) scored across 5 family-grouped folds.
- Leakage-safe rule held: all means and neighbors come from train rows only; held-out proteins get zero protein effect, so the additive model equals M1 (train ligand mean).
- Comparator is not trivial: M1 MSE 0.647 vs global mean 0.831 (gap 0.184, far above the 0.02 stop threshold). M1 CI 0.744.
- M2 (sequence 3-mer kNN, k=5, the one interaction-aware model) MSE 0.702, CI 0.720. It is WORSE than the ligand mean by 0.055 MSE (bootstrap-over-proteins 95% CI 0.026 to 0.082). So under family hold-out this crude sequence-similarity model adds nothing over the ligand baseline.
- Reference only (H3, different split, not recomputed): setting-1 additive 0.558, ligand-only 0.645. M1 is nearly the same under both splits (0.647 vs 0.645), expected because it uses no protein information. The additive protein effect is what family hold-out removes: setting-1 additive 0.558 vs ligand-only 0.645 on that fold. This crosses different test sets and is not a like-for-like family-hold-out gap; no claim of "harder" is made.
- Per-fold MSE varies a lot (M1 0.405 to 0.903), so a single pooled number hides large fold dependence.
Caveats: diagnostic study only, no novelty or benchmark claim, M2 is a deliberately simple fixed model (a negative for it says nothing about published DTA models). Y unverified against historical labels. KinHub family labels are an external annotation (pinned sha256); mapping of aliases was automatic and 61 entries are unmapped. R3 unchanged.
