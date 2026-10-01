# DAVIS R3: locked negative at the 40-chunk cap

Verified September 30, 2026. The locked R3 experiment stopped at chunk 40. This is experiment completion, not project completion or a benchmark beat.

## Measured result

- Final validation CI: 0.7337111904535643.
- Best validation CI: 0.7419244919292435.
- Preregistered validation falsification gate >= 0.80: failed.
- Single terminal test CI: 0.7233740663239825. Required > 0.878: failed.
- Single terminal test MSE: 1.19143545627594. Required < 0.261: failed.
- Primary gate MET: false.
- Stop reason: chunk_cap_40; checkpoint done=true, chunk=40, batch=0, stall=2.
- Final chunk train MSE: 2.2679; validation MSE: 1.1507.
- Original terminal output reports 348 domain-cropped proteins, 94 full-length fallbacks, and 38 context-capped at 1022.

The terminal test is recorded once in the original training output and equals the verdict JSON. It has not been rerun. The training output resumes from chunk 39 and records chunk 40 once. The complete 39-row incoming log prefix matches the previously authenticated chunk-39 artifact exactly.

## Authentication provenance

Source workflow run: 36675298621. Continuity artifact: 11080352957. Terminal evidence artifact: 11080367964. Both downloaded ZIP digests were checked locally against the reported artifact digests.

- Continuity ZIP SHA256: 225c476be9d85441620f46b923d9c4e38e46de9f3899b17cc4a996ea4542143a.
- Terminal ZIP SHA256: 91e78a353642c287f959a1724fad80c49058df0c435c964e4aa1e01f9929cfbe.
- Ciphertext SHA256: 7fa4aa958bcd01663fec55efcf58d8c2a24675b334ebaa5f70b20edcb8754ad8.
- Keyed ciphertext HMAC-SHA256: 566afa1f31d68dcacb1204fc1b3c0eb568c13b050d384ca1c361fe8035a3f6d1.
- Plaintext checkpoint SHA256: ca48cea69c2ce5dfec0087d5f75ac13cccb9b13898ff6d7eb5ed59b931eb59e6.

Ciphertext SHA and keyed HMAC were checked before decryption. The plaintext SHA, checkpoint state, status JSON and metrics log were then checked together. Evidence-only workflow commit 3090be0b3489e610ce70e5d34ae43731f00d65a5 was reported to add terminal-file upload; scientific trainer and protocol were not changed by this verification work.

## Interpretation and open gaps

R3 does not meet its validation or primary test gates. Small validation gains do not reverse this negative. No chunk 41, threshold change, restart or test reselection is justified by these data. Published matched-comparator evidence and full project criteria are not established by completing this run.

The negative should be included alongside the earlier R1/R2 results in the paper. Current-literature review and steering consultation remain open before selecting a different hypothesis and locking any new experiment. Full paper, novelty, dataset/tool/formula breadth, and benchmark criteria remain open. No public posting or remote push was performed for this audit.
