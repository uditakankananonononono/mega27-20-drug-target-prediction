# MEGA27-20 - drug-target prediction - AUDIT + COMPLETION PLAN
State from read API: src/targetscan package (data/models/train/metrics/external/chembl/audit), DAVIS ckpts,
mutant sensitivity+rescue, repurposing_candidates.json, ChEMBL corroboration, 20+ external tools registered in
docs/TOOLS_AND_DATASETS.md (counting rule: accession-level, used-vs-staged tiers), paper/main.tex+pdf exist.

## Audit checklist (after clone)
1. Run tests + reproduce DAVIS eval; record current CI/MSE.
2. Benchmark frame: DeepDTA DAVIS CI~0.878 MSE~0.261 (verify exact + split), GraphDTA CI~0.893 (verify).
   Beat = same split, same metrics, better numbers, seed variance reported.
3. Discovery audit: what in repurposing_candidates/mutant_rescue is NOVEL vs rediscovery (BI-2536->PLK1 is
   rediscovery - not a discovery). A discovery needs a named novel nomination + independent corroboration.
4. Gap count vs program bars: datasets, 200 tools, formulas, 50+ text-body pages, judge rounds 0/10.

## Completion gates: benchmark beat (locked comparator), 1+ novel corroborated nomination (discovery),
10+ judge rounds w/ fixes, 50+pp paper, negatives preserved (incl. corrected 3D5U ligand claim).
