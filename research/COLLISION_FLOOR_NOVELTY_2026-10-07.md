# Exact-input collision floor: novelty screen before implementation

This is new-hypothesis preparation, not an amendment or restart of R3. No locked predictive gate changes and no test reselection.

## Prior art checked

- Anderson and Bjarnadottir, *As good as it gets? A new approach to estimating possible prediction performance*, PLOS ONE, 2024. Primary page fetched: https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0296904 . Identical-input equal-prediction constraints, achievable-error diagnostics and an MSE extension are already discussed. The exact within-group squared-error decomposition is not a new theorem.
- Modification-aware DAVIS: https://arxiv.org/html/2512.00708 and the accompanying implementation https://github.com/ZhiGroup/DAVIS-complete . Curation of substitutions, insertions, deletions and phosphorylation precedes this lane. A diagnostic cannot claim modification-aware DAVIS as new.
- Broader protein-ligand bias work: https://pmc.ncbi.nlm.nih.gov/articles/PMC12552120/ . Dataset bias and leakage are already active benchmark topics, not an unoccupied novelty area.

## Candidate useful scope

An implementation-only diagnostic may partition a fixed observed set by the exact full input presented to a deterministic predictor and report within-group sum of squares divided by N. Completing the square proves that any prediction constant within each group has empirical MSE at least that quantity. This is a finite-sample oracle floor on those observed rows, not a generalization guarantee. A group key must include the ligand and all target inputs, not merely a target name or ligand. Hidden metadata, target IDs, sequence truncation and preprocessing alter equivalence classes. The grouping must be specified, not inferred from approximate similarity.

## Decision and remaining gaps

Proceed only as a tested engineering utility labeled replication/basic decomposition. No novelty credit, scientific completion count, measured DAVIS floor or benchmark superiority is claimed. Real-data use needs the exact input representation, raw-row provenance, fixed row scope and documented comparator. Censored affinity labels limit biological interpretation. The floor must never be used to soften existing thresholds or explain R3's failure without the corresponding evidence.
