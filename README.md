# mega27-20-drug-target-prediction
Drug-target affinity prediction on DAVIS (68 ligands x 442 kinases) with preregistered rounds R1, R2, R3 (`PREREG_DAVIS_*.md`).

Status (Oct 7, 2026):
- R3 is closed negative; it is preserved, not restarted, and no threshold or test was changed.
- No novel result and no world-benchmark win is claimed. Mutant/domain ideas overlap prior art (DAVIS-complete); see `research/COLLISION_FLOOR_NOVELTY_2026-10-07.md`.
- `src/targetscan/collision_floor.py` (10 stdlib tests) computes an exact-input squared-error floor for deterministic ligand+sequence-only models. On the pinned public DAVIS mirror (`data/davis_pin/SOURCES_SHA256.md`) the floor is 0.0322 in pKd units (30056 rows, 546 conflicting input groups; 442 proteins, 379 unique sequences). It does not bind against typical DAVIS error. Diagnostic only, not a benchmark gate: `research/DAVIS_COLLISION_FLOOR_RESULT_2026-10-07.md`.
- Mirror ligand and protein files equal the historical data_cache files exactly. The label array Y is NOT verified against any historical label array (none in the repo). Pearson r 0.375 between old predictions and mirror pKd is a sanity check only, not verification.
- Descriptive audits (plans frozen before runs, no novelty claim): `research/RESULT_H1H2H3_AUDIT_2026-10-08.md` (fold overlap, additive baseline MSE 0.558, floors) and `research/RESULT_H4_FAMILY_HOLDOUT_2026-10-08.md` (family hold-out: ligand-mean MSE 0.647, sequence-kNN 0.702, kNN worse by 0.055).
- Draft hypothesis shortlist (no lock, no run): `research/HYPOTHESIS_SHORTLIST_2026-10-08.md`.
- pytest, torch and rdkit are not installed in the audit environment, so the full test suite was not run there.
- Paper: `paper/` (6 pages, scoped claims).
