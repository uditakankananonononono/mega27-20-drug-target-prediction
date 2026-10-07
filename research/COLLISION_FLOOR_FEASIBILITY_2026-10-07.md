# Feasibility boundary, October 7

The current clone has tracked original and corrected sequence maps, Pfam crop metadata and ligand strings. A direct static audit of the tracked original map finds 54 identical WT/mutant entries, consistent with the historical summary. This is an input-cache consistency check, not a rerun of affinity outcomes.

The checked data_cache directory lacks the raw affinity matrix and fixed train/test fold files. Historical per-target summaries do not contain all label vectors and cannot reconstruct a squared-error floor. A new fetch would need a pinned source checksum and explicit scope; it cannot be represented as recovery of the historical R3 inputs without matching the original artifacts. Corrected sequences and domain crops also mean that original-map collision groups cannot explain R3's outcome.

No real-data floor has been calculated. Next requirements: exact input recipe and fixed-row manifest, raw label array checksum, explicit original-distribution diagnostic rather than R3 performance reanalysis, and independent provenance validation. Until these are recovered, retain the synthetic tested utility without empirical scientific claims. No restart, additional terminal test, threshold adjustment or test reselection follows from this gap.
