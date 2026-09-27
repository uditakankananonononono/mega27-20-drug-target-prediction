# PREREG_DAVIS_R3 ADDENDUM 1 - fallback context cap (locked 2026-09-28 ~00:51 IST,
# BEFORE the first locked R3 chunk; ratified via parent 00:49 IST)

## Deviation
PREREG_DAVIS_R3.md says full-length fallback kinases (no resolvable Pfam
PF00069 domain) "keep full length". 38 of the 94 fallbacks exceed the ESM-2
t6_8M hard context of 1024 tokens (1022 residues + BOS/EOS) and physically
cannot be encoded at full length. These 38 are CAPPED at their first 1022
residues. Exclusion was considered and rejected: it would silently change the
locked DAVIS benchmark set, which is worse.

## Convention match
Identical to the repo's existing src/precompute_esm2.py path (MAXLEN = 1022,
truncations recorded), so R2's frozen-embedding cache and R3's live encoder
treat over-length sequences the same way.

## The exact 38 capped targets (name | full length)
| target | full_len |
|---|---|
| MTOR | 2548 |
| MRCKA | 1731 |
| MRCKB | 1710 |
| PIK3C2B | 1633 |
| PIK3C2G | 1444 |
| MAP4K4 | 1272 |
| QSK | 1262 |
| ABL1(E255K) | 1166 |
| ABL1(F317I) | 1166 |
| ABL1(F317I)p | 1166 |
| ABL1(F317L) | 1166 |
| ABL1(F317L)p | 1166 |
| ABL1(H396P) | 1166 |
| ABL1(H396P)p | 1166 |
| ABL1(M351T) | 1166 |
| ABL1(Q252H) | 1166 |
| ABL1(Q252H)p | 1166 |
| ABL1(T315I) | 1166 |
| ABL1(T315I)p | 1166 |
| ABL1(Y253F) | 1166 |
| ABL1 | 1166 |
| ABL1p | 1166 |
| ABL2 | 1166 |
| TIE2 | 1156 |
| PIK3CG | 1101 |
| PIK3CB | 1069 |
| PIK3CA | 1068 |
| PIK3CA(C420R) | 1068 |
| PIK3CA(E542K) | 1068 |
| PIK3CA(E545A) | 1068 |
| PIK3CA(E545K) | 1068 |
| PIK3CA(H1047L) | 1068 |
| PIK3CA(H1047Y) | 1068 |
| PIK3CA(I800L) | 1068 |
| PIK3CA(M1043I) | 1068 |
| PIK3CA(Q546K) | 1068 |
| PIK3CD | 1043 |
| CDKL5 | 1029 |

## Unchanged locks
All other R3 locks stand verbatim: pKd labels, DeepDTA setting1 splits,
1200-pair rng-0 val split, 8k-pair chunks, batch 256 optimizer steps
(implemented as length-sorted adaptive micro-batches with gradient
accumulation - optimizer-step semantics identical), locked early stop
(no >= 0.002 val_ci improvement over 8 consecutive chunks), falsification
gate val_ci < 0.80 at the earlier of chunk 40 or early stop, single locked
test evaluation on the stopped checkpoint, primary gate CI > 0.878 AND
MSE < 0.261.
