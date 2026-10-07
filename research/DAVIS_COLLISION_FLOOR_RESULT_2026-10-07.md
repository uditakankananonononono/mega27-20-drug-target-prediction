# DAVIS exact-input collision floor (pinned DeepDTA-mirror data), Oct 7
Inputs: sha256-pinned files in data/davis_pin/SOURCES_SHA256.md. Key = ligand SMILES + NUL + protein sequence bytes. All 68x442 = 30056 rows (no NaN).
- Unique ligands 68; proteins 442 but only 379 unique sequences (63 duplicate-sequence entries, e.g. phospho/mutant variants share sequence).
- Exact-input groups 25772; groups with conflicting labels 546.
- Empirical MSE floor, raw Kd (nM): 559622.46 (scale-dominated, not interpretable).
- Empirical MSE floor, pKd = -log10(Kd/1e9): 0.0322.
Scope: fixed observed rows, deterministic ligand+sequence-only models. Not a population bound, not a benchmark gate, not a novel result (prior art: Anderson & Bjarnadottir 2024, DAVIS-complete). A floor of 0.032 is far below typical reported DAVIS MSE, so it does not bind. This is a diagnostic on the public mirror, not a rerun of R3; R3 stays closed negative. Not yet checked: mirror Y vs our historical DAVIS artifacts.
