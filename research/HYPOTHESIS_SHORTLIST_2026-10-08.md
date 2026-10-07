# Lane 20 hypothesis shortlist (draft for Main's pick). No lock, no run, no R3 restart.
Prior art already checked (see sources.md, COLLISION_FLOOR_NOVELTY): DAVIS-complete (arxiv 2512.00708, github ZhiGroup/DAVIS-complete), Anderson & Bjarnadottir 2024 (PLOS ONE 0296904), Nature MI 2025 bias paper (s42256-025-01124-5), PMC12552120. Anything below that touches those is a replication until a novelty screen says otherwise. None of these is claimed novel yet.

Data already pinned: DAVIS Y (68x442), ligands, proteins, DeepDTA fold setting 1 (sha256 in data/davis_pin). Label identity vs historical R1-R3 labels NOT verified.

## H1. Protein-duplicate leakage in the DeepDTA fold
Claim to test: 442 entries map to 379 unique sequences (63 duplicates), so some test rows share an exact (ligand, sequence) input with a training row; a model's reported error is partly memorized-label lookup.
- Needs: pinned fold files only (already have). Count test rows whose exact input appears in train, and the label disagreement among those pairs.
- Matched comparator: same count under random row folds; and DAVIS-complete's own curation for overlap (must read first; likely partly covered, so expect replication).
- Cost: seconds, stdlib. Output is descriptive, not a benchmark claim.
- Risk: probably known; screen DAVIS-complete and 2512.00708 before locking.

## H2. Collision-floor-aware reporting
Claim to test: report each published/own model's DAVIS MSE alongside the 0.0322 pKd floor and the fraction of rows in conflicting groups (546 groups).
- Needs: nothing new. Matched comparator: floor on shuffled-label control (should be near label variance) and on a ligand-only key (shows what each input omits).
- Status: diagnostic, prior art is Anderson & Bjarnadottir; value is documentation, no novelty credit.

## H3. Ligand-only and protein-only baselines as the matched comparator for any "interaction" claim
Claim to test: a ligand-mean plus protein-mean additive model (no interaction) recovers X% of explained variance on the DeepDTA test fold; any deep gain must be reported as the increment over it.
- Needs: Y and fold files (have). Matched comparator: additive model fit on the same train rows, same MSE and CI index metric. Cheap, CPU, stdlib or numpy.
- Prior art: additive/baseline-bias results exist in the bias literature (PMC12552120, Nature MI 2025); treat as replication. It is a useful floor for the paper.

## H4. Cold-protein split with family-level hold-out
Claim to test: error rises when entire kinase families are held out (setting 3-like), and the rise exceeds what ligand-only and additive baselines show.
- Needs: kinase family labels (UniProt/KLIFS, not in repo; pull with pinned checksum), pinned Y.
- Comparator: H3 additive model on the same family splits; random-fold same-size test.
- Risk: moderate prior art on cold-target splits. Needs a new lock and a novelty screen. Free CPU if the model is classical (ridge on fingerprints plus k-mer/Pfam features).

## H5. Mutant/phospho-entry sensitivity
Claim to test: DAVIS entries sharing a sequence with a wild-type kinase carry label disagreement larger than assay noise.
- Needs: DAVIS-complete modification annotations (fetch with checksum). Overlaps DAVIS-complete directly; not novel. Included only as a data-quality audit.

## Recommendation
Do H1 + H3 + H2 as one pre-lock descriptive audit (hours, CPU, no training). It gives the paper a real measured floor and baseline and settles the leakage question. H4 is the only one that could carry a new scientific claim and needs a lock plus a novelty screen first. Pick before anything is run. Nothing here changes R3, thresholds, or test selection.
