#!/usr/bin/env bash
set -euo pipefail
umask 077
: "${R3_CKPT_TRANSFER_KEY:?missing transfer key}"
mkdir -p results .r3-private
base='https://github.com/uditakankananonononono/mega27-20-drug-target-prediction/releases/download/r3-ckpt-c23-b7168'
for spec in \
  'r3-ckpt-c23-b7168.enc-7d218b40.part_00 14f886fc63f305cb7b5309fc21f0e7033370c7b4f58db5b0766908bb3e6a27f9' \
  'r3-ckpt-c23-b7168.enc-bb723315.part_01 ded14d599ba5cb17a83d9c3559ba2e252e890d4709c1825985f0988c54233fb8' \
  'r3-ckpt-c23-b7168.enc-1a39d288.part_02 fef71c61181ebb0927c6210953bbff2f893c0779298cc9bc56ec8ff88009196a'; do
  read -r name digest <<< "$spec"
  curl -fsSL --retry 3 "$base/$name" -o ".r3-private/$name"
  printf '%s  %s\n' "$digest" ".r3-private/$name" | sha256sum -c - >/dev/null
done
cat .r3-private/*.part_00 .r3-private/*.part_01 .r3-private/*.part_02 > .r3-private/bootstrap.enc
printf '%s  %s\n' '8e09dff63f9e33438f821b9757882625d4d68612eb63c30c15cf4c697d394279' '.r3-private/bootstrap.enc' | sha256sum -c - >/dev/null
python3 - <<'PY'
import hashlib,hmac,os
k=os.environ['R3_CKPT_TRANSFER_KEY'].encode()
v=hmac.new(k,open('.r3-private/bootstrap.enc','rb').read(),hashlib.sha256).hexdigest()
if not hmac.compare_digest(v,'eac310faea176a74b9adfaa674a005348fcf4212d4fe6a44f101809e343dadc6'):
    raise SystemExit('checkpoint authentication failed')
PY
# Never echo commands carrying the key, and never print the checkpoint or environment.
mkdir -p "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints"
curl -fsSL --retry 3 https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t6_8M_UR50D.pt -o "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints/esm2_t6_8M_UR50D.pt"
printf '%s  %s\n' '46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928' "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints/esm2_t6_8M_UR50D.pt" | sha256sum -c - >/dev/null
printf '%s' "$R3_CKPT_TRANSFER_KEY" > .r3-private/key
openssl enc -d -aes-256-cbc -pbkdf2 -iter 600000 -in .r3-private/bootstrap.enc -out results/davis_r3_ckpt.pt -pass file:.r3-private/key
printf '%s  %s\n' 'f9bfbd53389144260f952301ea57191cf71b38d829af8a3c234034ac3868a2b7' 'results/davis_r3_ckpt.pt' | sha256sum -c - >/dev/null
python3 - <<'PY'
import torch
x=torch.load('results/davis_r3_ckpt.pt',map_location='cpu',weights_only=False)
assert x['chunk']==22 and x['batch_in_chunk']==7168 and x['stall']==0 and not x['done']
assert abs(float(x['best_val'])-0.7127747383282877)<1e-12
PY
# The benchmark dataset files are fetched by load_davis, then checked before training.
python3 - <<'PY'
import sys
sys.path.insert(0,'src')
from targetscan.data.davis_kiba import load_davis
load_davis()
PY
printf '%s\n' \
  '9855c4f234ec6dd18295b0e0cd6ee54214cafa908c6eb1d447b47d454fc0294a  data_cache/davis_Y' \
  'bb7a63b2178c40a4e2417bc701b95511407787fa5eabb6a43a695fe6bd722f56  data_cache/davis_folds_train_fold_setting1.txt' \
  'c97589fa2fae318b01870854af0ab6fdd155dd38524cc6ac1d0319418ac6d568  data_cache/davis_folds_test_fold_setting1.txt' | sha256sum -c - >/dev/null
ESM_REPO_REF=2b369911bb5b4b0dda914521b9475cad1656b2ac python3 run_davis_r3.py 1 > .r3-private/training.out 2>&1
python3 - <<'PY'
import json,torch
rows=[json.loads(x) for x in open('results/davis_r3_log.jsonl')]
assert rows[-1]['chunk']==23
x=torch.load('results/davis_r3_ckpt.pt',map_location='cpu',weights_only=False)
assert x['chunk']==23 and x['batch_in_chunk']==0
PY
# No plain checkpoint in an artifact: encrypt output for independent box comparison.
openssl enc -aes-256-cbc -salt -pbkdf2 -iter 600000 -in results/davis_r3_ckpt.pt -out .r3-private/comparison.enc -pass file:.r3-private/key
sha256sum .r3-private/comparison.enc > .r3-private/comparison.sha256
python3 - <<'PY'
import hashlib,hmac,os
k=os.environ['R3_CKPT_TRANSFER_KEY'].encode()
p='.r3-private/comparison.enc'
with open('.r3-private/comparison.hmac','w') as f:
    f.write(hmac.new(k,open(p,'rb').read(),hashlib.sha256).hexdigest()+'\n')
PY
cp results/davis_r3_log.jsonl .r3-private/comparison-log.jsonl
rm -f .r3-private/key results/davis_r3_ckpt.pt
