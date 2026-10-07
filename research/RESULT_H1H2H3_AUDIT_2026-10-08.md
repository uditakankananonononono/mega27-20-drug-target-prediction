# Result: descriptive audit H1+H3+H2 (plan frozen in PLAN_H1H2H3_AUDIT_2026-10-08.md, commit before this run). One run, no tuning. Script: src/targetscan/audit_h1h2h3.py (needs the pinned files in the working dir).
Raw output:
```
train structure <class 'list'> 5 [5010, 5009, 5009, 5009, 5009] test 5010
pool 25046 overlap pool/test 0
H1 fold (929, 0.3302680458785049)
dup protein entries 63
H1 random 100 splits hit count min/med/max 823 875.0 957
H3 additive MSE 0.5580259866682145 CI 0.7919917395583965
H3 ligand-only MSE 0.6447119218235383 CI 0.7434671586338354
H3 protein-only MSE 0.7189846018132262 CI 0.6637645259550102
H3 global-mean MSE 0.8015143274995329 CI 0.5 by construction
H2 floor pool 0.02931441500022997 test 0.0201090975671399
H2 ligand-only key 0.6192204727241335 protein-seq-only key 0.7027010193684247
H2 shuffled median 0.11282804983601426 label var 0.800491719138235
H5 dup-seq groups 18 ligand-groups differing 546 median diff 0.9586073148417746 max 4.745270023988987
```
Reading, measured only:
- Fold: DeepDTA setting 1, 5 train lists (~5009 each, pool 25046 rows) and a 5010-row test list. Pool and test do not overlap by row index.
- H1: 929 of 5010 test rows (18.5%) share an exact (ligand, sequence) key with a train row, because 63 protein entries duplicate another entry's sequence. 100 random splits of the same size give 823 / 875 / 957 (min/median/max), so this fold is inside the random range, not unusual. Those 929 rows differ from the train label for the same key by mean squared 0.330 pKd^2, so label lookup does not give a zero error.
- H3: additive ligand+protein model test MSE 0.558, CI 0.792; ligand-only 0.645 / 0.743; protein-only 0.719 / 0.664; global mean 0.802. So an additive no-interaction model already explains about 30% of test variance. Any interaction claim must be reported as the increment over this.
- H2: floor on train pool 0.0293, test 0.0201, all rows 0.0322. Ligand-only key floor 0.619, protein-sequence-only key floor 0.703, shuffled-label floor median 0.113 (label variance 0.800). A full-input floor of 0.03 is far below any reported DAVIS error and does not bind. The shuffled control shows the floor is low partly because most keys are unique (small groups), so a low floor is expected from group size, not only from label consistency.
- H5 (descriptive): 18 duplicate-sequence groups; all 546 conflicting exact-input groups come from them; median label difference within a conflict 0.959 pKd, max 4.745. These are assay-level differences between entries that share a sequence (e.g. variants or phosphorylated forms); not resolved here, DAVIS-complete covers it.
Caveats: label array Y still unverified against historical labels. Comparable published DAVIS numbers (e.g. DeepDTA) were not re-fetched here, so I make no comparison to them. Diagnostic only; no novelty, no benchmark claim, R3 unchanged.
