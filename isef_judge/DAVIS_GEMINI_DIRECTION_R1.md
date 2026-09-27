# LANE-20 DIRECTION CONSULT (model: GEMINI free tier) - logged 2026-09-27 09:00 IST

Model surface: Google Gemini (free tier), her Google account via cloud browser
config-c saved profile. Conversation: https://gemini.google.com/app/3a0791f99c3b4e08
Verbatim response: isef_judge/davis_gemini_direction_capture.txt
Prompt: staged /tmp/gem_l20_01..04 + gem_l20_q (4 paper parts + question incl.
R1 falsification outcome ep130 val_ci 0.7610 < 0.78; direction-rule Q4 added).
Authority: user WhatsApp 8:55:01 AM verbatim ("...not to provide any failures
and always pivot to a new direction within that topic by asking an llm.",
wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhggQUMyRkY0M0JBODE5OEE3RjE4NTVBRDA3RUE5MkYyMjQA)
relayed by parent 08:55 with explicit Gemini mandate; Gemini surface authorized
12:11:28 AM (wamid...NDFGMzcA) + simultaneous-surfaces 12:13:07 AM.

JUDGE COUNT NOTE: per the bound master spec ChatGPT is THE mandatory judge.
This Gemini consult is an ADDED-SURFACE direction round; it does NOT increment
lane-20's 10-round ChatGPT judge count (still 0/10).

## Gemini's direction (headline; full text in capture)
1. Weaknesses: (a) CI gap vs DeepDTA under identical metrics; (b) 8k-chunk
   training protocol non-standard, needs ablation vs full-batch; (c) random
   split discipline - cold-target/cold-drug splits unverified.
2. Strongest R2 improvement: pretrained pLM (ESM-2 frozen embeddings) for
   protein + GATv2 drug graph; evidence = CI>=0.878/MSE<=0.261 parity +
   mutation-sensitivity Spearman>=0.50 on corrected sequences.
3. Most novel contribution: automated open-source DTA benchmark data-integrity
   audit suite (extends our committed 54-mutant byte-identical defect finding).
4. PIVOT ANGLE (headline): mutation-aware affinity modeling on corrected
   DAVIS-variant sequences (UniProt-sourced true mutant sequences, delta-pKd
   prediction, ABL1 T315I-style resistance shifts).

## Decision taken (adopted vs parked, verified against repo state)
- ADOPT pivot step 1 (corrected DAVIS-variant dataset via UniProt REST): free
  API, extends committed audit results/benchmark_audit_davis_mutants.json +
  mutant_sensitivity.json - strongest signal, becomes R2 prereg core.
- R2 prereg (PREREG_DAVIS_R2.md) will lock BEFORE any R2 training per RULE 6:
  ESM-2 (free downloadable) frozen embeddings + existing GNN, corrected-
  sequence dataset, delta-pKd evaluation subset, same falsification discipline.
- PARKED: full-epoch ablation (compute-infeasible in sandbox; honest
  limitations note), cold-split evaluation (R2 secondary if compute allows).
