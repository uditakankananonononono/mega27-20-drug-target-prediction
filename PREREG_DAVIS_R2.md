# PREREG: DAVIS rung R2 (locked 2026-09-27 09:11 IST, BEFORE any R2 training)

Escalation path required by PREREG_DAVIS_R1.md falsification clause (triggered:
val_ci 0.7610 < 0.78 @ ep130; verdict results/DAVIS_R1_VERDICT.md @ aeb2f51)
+ user standing rule (WhatsApp 8:55:01 verified: negatives never terminal,
pivot within topic by asking an LLM) + Gemini direction consult (added surface,
verbatim isef_judge/davis_gemini_direction_capture.txt @ d56a851).

## What R2 changes (two coupled interventions, locked together)
1. DATASET CORRECTION (the pivot core): train/evaluate on
   data_cache/davis_proteins_corrected.json - the 45 mutant entries repaired
   by src/pull_uniprot_variants.py (20 direct + 17 offset + 5 phospho-reuse =
   42 corrected, UniProt-validated, provenance results/uniprot_variant_pull_provenance.json
   @ 8300e5e). RET wt + 3 mutants: distributed sequence lacks the kinase
   domain (second dataset defect found during correction); R2 uses FULL
   UniProt canonical RET (P07949) for RET wt and all RET mutants, and flags
   these 4 entries in provenance. WT entries otherwise byte-identical to R1.
2. PROTEIN ENCODER: replace the from-scratch 1D CNN with FROZEN pretrained
   ESM-2 (esm2_t6_8M_UR50D, 8M params, free download from
   dl.fbaipublicfiles.com via torch.hub facebookresearch/esm). Per-residue
   embeddings precomputed ONCE for all 442 corrected targets and cached;
   attention-pool (small learned pooler) -> existing drug tower + affinity head.
   No protein-LM fine-tuning (sandbox RAM); only pooler+drug tower+head train.

## Endpoints (all locked; evaluated on the same DeepDTA DAVIS test fold as R1)
- PRIMARY (benchmark): test CI > 0.878 AND test MSE < 0.261 (DeepDTA bar).
- SECONDARY (mutation sensitivity, the pivot's headline evidence): on mutant/WT
  drug pairs with |true delta pKd| >= 1.0 (323 pairs per committed audit),
  Spearman(predicted delta, true delta) >= 0.50. R1's published-equivalent
  model scores ~0 here BY CONSTRUCTION (mutant inputs were byte-identical to
  WT) - this contrast is the paper's key exhibit, not a new metric invented
  post hoc: the pair subset was enumerated in the audit BEFORE R2.
- FALSIFICATION GATE: if val_ci has not reached 0.80 by absolute R2 epoch 40
  (full-data-equivalent chunks counted as in R1), declare R2 an honest
  negative and escalate to R3 (next LLM redirection round).

## Sandbox contingencies (decision rules, locked now)
- If the ESM-2 weights download fails from this sandbox: fall back to
  intervention 1 alone (corrected sequences, R1 architecture retrained from
  scratch, same endpoints/gates) and record the fallback in the verdict.
- If ESM-2 embedding precompute exceeds sandbox budget: truncate protein
  inputs to the kinase-domain window used for the offset correction (recorded
  per entry in provenance) - decision recorded before training, not after.
- Chunk protocol as R1 (8k-pair chunks, Adam lr 1e-3, bs 256), batch-level
  checkpointing; all chunks logged; no cherry-picked restarts. Negatives stay.

## Honesty
R1's negative remains in the paper (compact limitations treatment per user
rule 8:27:49). The mutation-sensitivity contrast uses the pre-enumerated
323-pair subset. Any endpoint change after this lock voids the rung.
