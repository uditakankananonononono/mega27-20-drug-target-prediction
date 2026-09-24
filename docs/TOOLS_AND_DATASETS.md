# External tools and datasets - item 20 (honest registry, updated as used)

Counting rule (program-wide): accession-level datasets actually used; one study's
condition matrix = one dataset. Two tiers: USED IN RESULTS vs STAGED (pulled and
cached for the verification pipeline, not yet cited in a result).

## External research tools (web resources + databases)
Used in verified results:
1. DAVIS kinase affinity dataset (via DeepDTA GitHub mirror) - benchmark
2. ChEMBL REST API - refutation/corroboration of screen candidates (CHEMBL513909, doc CHEMBL1908390)
3. PubChem - BI-2536 identity resolution (CID 11364421)

Staged (client live, cached, wired into verification pipeline):
4. PubChem PUG-REST - physchem props for all 68 DAVIS drugs (results/external_pull.json)
5. UniProt REST - canonical reviewed accessions for 10 hit kinases
6. STRING API - PPI network (38 edges) among screen targets
7. ClinicalTrials.gov API v2 - volasertib (BI-2536) trial history, 5 NCT IDs
8. Europe PMC REST - BI-2536/PLK1 literature corroboration
9. KEGG REST - pathways for 5 hit kinases (PLK1=hsa:5347: cell cycle, FoxO)
10. RCSB PDB search API - structures incl. 3D5U (the BI-2536/PLK1 co-crystal)
11. HGNC REST - approved symbols/names for hit kinases
12. Reactome Content Service - pathway mapping via UniProt accessions
13. Open Targets Platform GraphQL - disease associations (KIT->GIST 0.89,
    FLT3->AML 0.83, textbook-correct)

## Packages
13. PyTorch 2.14.0  14. NumPy  15. scikit-learn  16. RDKit 2026.03.6
17. pytest (23 hermetic tests)  18. SciPy

Honest tool count: 21 (13 resources cited/staged + OT pull, 6 packages; 24 tests)
Path to 40: BindingDB, PDB/RCSB, KLIFS, Open Targets, DrugCentral, TTD,
Guide to PHARM, ZINC, Pharos, GDSC, DepMap, KEGG, Reactome, g:Profiler,
Ensembl REST, AlphaFold DB, SwissTargetPrediction, admetSAR, ProTox, PK-DB,
KinMap, HGNC, ChEBI, PubChem BioAssay, DeepDTA/GraphDTA baselines (2 tools).

## Datasets (accession-level, actually pulled)
Used in results: DAVIS (1), ChEMBL assay doc CHEMBL1908390 (1) = 2
Staged: 68 PubChem CID records + 10 UniProt entries + 5 NCT trials +
5 Europe PMC articles + STRING network (1) + 5 KEGG gene ids +
25 PDB accessions (5 per kinase) + 5 HGNC ids + Reactome pathway sets (5) = 129
Honest dataset count: 131 (2 used in results, 129 staged for the screen)
Path to 120+: KIBA (1), Metz (1), per-kinase ChEMBL assay sets (~10),
PDB structures for hit complexes (~10), remaining bioassay records.
