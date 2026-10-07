# DAVIS exact-input collision floor (pinned DeepDTA-mirror data), Oct 7
Inputs: sha256-pinned files in data/davis_pin/SOURCES_SHA256.md. Key = ligand SMILES + NUL + protein sequence bytes. All 68x442 = 30056 rows (no NaN).
- Unique ligands 68; proteins 442 but only 379 unique sequences (63 duplicate-sequence entries, e.g. phospho/mutant variants share sequence).
- Exact-input groups 25772; groups with conflicting labels 546.
- Empirical MSE floor, raw Kd (nM): 559622.46 (scale-dominated, not interpretable).
- Empirical MSE floor, pKd = -log10(Kd/1e9): 0.0322.
Scope: fixed observed rows, deterministic ligand+sequence-only models. Not a population bound, not a benchmark gate, not a novel result (prior art: Anderson & Bjarnadottir 2024, DAVIS-complete). A floor of 0.032 is far below typical reported DAVIS MSE, so it does not bind. This is a diagnostic on the public mirror, not a rerun of R3; R3 stays closed negative. Not yet checked: mirror Y vs our historical DAVIS artifacts.

## Cross-check update (Oct 8)
Historical data_cache/davis_ligands_can.txt and davis_proteins.txt equal the pinned mirror files exactly (same sha256 prefix, same ordered 68 ligands and 442 proteins). Labels Y remain unverified against any historical label array (none in repo); results/davis_full_pred.npy is predictions only.
Sanity only: Pearson r between davis_full_pred.npy and mirror pKd is 0.375 (68x442), sign-consistent with a weak real-label fit, not proof of label identity.
