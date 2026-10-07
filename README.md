# mega27-20-drug-target-prediction
Drug-target affinity prediction on DAVIS (68 ligands x 442 kinases) with preregistered rounds R1, R2, R3 (`PREREG_DAVIS_*.md`).

Status (Oct 7, 2026):
- R3 is closed negative; it is preserved, not restarted, and no threshold or test was changed.
- No novel result and no world-benchmark win is claimed. Mutant/domain ideas overlap prior art (DAVIS-complete); see `research/COLLISION_FLOOR_NOVELTY_2026-10-07.md`.
- `src/targetscan/collision_floor.py` (10 stdlib tests) computes an exact-input squared-error floor for deterministic ligand+sequence-only models. On the pinned public DAVIS mirror (`data/davis_pin/SOURCES_SHA256.md`) the floor is 0.0322 in pKd units (30056 rows, 546 conflicting input groups; 442 proteins, 379 unique sequences). It does not bind against typical DAVIS error. Diagnostic only, not a benchmark gate: `research/DAVIS_COLLISION_FLOOR_RESULT_2026-10-07.md`.
- Mirror labels have not yet been cross-checked against the historical DAVIS artifacts used in R1-R3.
- pytest, torch and rdkit are not installed in the audit environment, so the full test suite was not run there.
- Paper: `paper/` (6 pages, scoped claims).
