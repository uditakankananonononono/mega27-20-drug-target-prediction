# Lane-20 locked amendment queue (from PROVIDED ROUND 1, wamid...NzcwNQA=)

Foldback discipline: every item lands as (critique + NOVELTY CHANGE + commit ref).
Items marked PRIOR-PLANNED were already in the R2 prereg/work before the verdict;
they satisfy the critique but are not counted as verdict-driven novelty.

W1 model novelty -> mutation-aware protein encoder + uncertainty-aware head.
   PARTIALLY PRIOR-PLANNED: R2 = frozen ESM-2 + learned attention-pool (ec26870).
   NEW: uncertainty head (W7), explicit mutation representation (W13).
W2 auditor as formal contribution -> NEW: package benchmark auditor as standalone
   tool run on DAVIS + KIBA + BindingDB + PDBbind with a dataset-integrity score;
   headline paper contribution (builds on our committed 54-mutant defect finding,
   results/benchmark_audit_davis_mutants.json).
W3 DAVIS-clean corrected retrain -> PRIOR-PLANNED: corrected dataset 42/45
   (8300e5e) + R2 trainer running (ec26870). NEW: locked original-vs-corrected
   comparison table (R1 0.7533 vs R2 final, same splits).
W4 external benchmark validation -> NEW: evaluate R2 protocol on KIBA
   (load_kiba already in src/targetscan/data/davis_kiba.py).
W5 identical-condition baselines -> NEW: retrain DeepDTA-protocol baseline under
   our splits/budget; matched-protocol comparator table.
W6 harder splits -> NEW: drug-disjoint, protein-disjoint (kinase-family),
   scaffold-split evaluations of the final model.
W7 uncertainty estimation -> NEW: MC dropout (dropout already in head) +
   3-seed ensemble variance; confidence reported per prediction.
W8 interpretability -> NEW: residue attention maps from R2 attention-pool,
   saliency on drug atoms; case studies on known binding sites.
W9 attention-failure analysis -> NEW: analyze R1 attn (0.7060) vs CNN (0.7080);
   attention-weight correlation with kinase domains.
W10 range compression -> NEW: focal/quantile/pairwise-ranking loss arm;
   preregister before eval (weighted-loss failure already documented).
W11 modern PLM comparison -> PARTIALLY PRIOR-PLANNED: R2 ESM-2 vs R1 CNN is
   this comparison. NEW: ProtT5 arm if CPU-feasible, else documented scope cut.
W12 ligand representations -> NEW: Morgan-FP baseline exists (run_davis_lr.py);
   add Graph-Transformer/SMILES-Transformer arm as feasible, else documented.
W13 explicit mutation representation -> NEW: wt + position + substitution
   features on top of ESM-2 embeddings.
W14 mutation window -> PRIOR-PLANNED/RESOLVED: R2 embeds FULL-LENGTH corrected
   sequences (data_cache/esm2/*.npy [L,320], no 600-residue window).
   KIT(D816V)-class artifact eliminated; document in paper.
W15 prospective blind challenge -> NEW: train on older drugs, predict newer
   kinase inhibitors (temporal split).
W16 external auditor release -> NEW: run auditor on public datasets/repos,
   report additional defects discovered.
W17 biological/clinical connection -> NEW: connect predictions to resistance
   mutations (e.g. ABL1 T315I/F317I), approved drugs, clinical kinase relevance.
W18 cross-protein-class transfer -> NEW: transfer test on a non-kinase dataset
   (GPCR/ion-channel), documented scope.

Order of execution (locked): W3/W11/W14 land with R2 (running); W2 auditor
tool + W16; W6 splits; W5 baselines; W10 losses; W7/W8/W9 analyses; W4/W18
external; W13/W15; W12/W17 paper integration. Paper rebuild absorbs all.
