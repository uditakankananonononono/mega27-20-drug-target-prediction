#!/usr/bin/env bash
set -euo pipefail
umask 077
: "${R3_CKPT_TRANSFER_KEY:?missing key}"
: "${SOURCE_RUN_ID:?missing source run}"
: "${SOURCE_CIPHER_SHA:?missing cipher sha}"
: "${SOURCE_HMAC:?missing hmac}"
: "${SOURCE_PLAIN_SHA:?missing plain sha}"
: "${SOURCE_CHUNK:?missing chunk}"
: "${EXPECTED_NEXT_CHUNK:?missing next chunk}"
: "${GH_TOKEN:?missing scoped actions token}"
[[ "$SOURCE_RUN_ID" =~ ^[0-9]+$ && "$SOURCE_CHUNK" =~ ^[0-9]+$ && "$EXPECTED_NEXT_CHUNK" =~ ^[0-9]+$ ]]
[[ "$EXPECTED_NEXT_CHUNK" -eq $((SOURCE_CHUNK+1)) ]]
for v in "$SOURCE_CIPHER_SHA" "$SOURCE_HMAC" "$SOURCE_PLAIN_SHA"; do [[ "$v" =~ ^[a-f0-9]{64}$ ]]; done
mkdir -p results .r3-private
# The prior run's artifact is encrypted; no plaintext checkpoint travels between jobs.
gh run download "$SOURCE_RUN_ID" -R uditakankananonononono/mega27-20-drug-target-prediction -n r3-continuity-encrypted -D .r3-private
printf '%s  %s\n' "$SOURCE_CIPHER_SHA" '.r3-private/comparison.enc' | sha256sum -c - >/dev/null
python3 - <<'PY'
import hashlib,hmac,os
key=os.environ['R3_CKPT_TRANSFER_KEY'].encode()
p='.r3-private/comparison.enc'
expected=os.environ['SOURCE_HMAC']
stored=open('.r3-private/comparison.hmac').read().strip()
observed=hmac.new(key,open(p,'rb').read(),hashlib.sha256).hexdigest()
if not (hmac.compare_digest(observed,expected) and hmac.compare_digest(observed,stored)):
    raise SystemExit('checkpoint authentication failed')
PY
# Install/fetch the same ESM weight as the verified comparison run.
mkdir -p "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints"
curl -fsSL --retry 3 https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t6_8M_UR50D.pt -o "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints/esm2_t6_8M_UR50D.pt"
printf '%s  %s\n' '46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928' "${TORCH_HOME:-$HOME/.cache/torch}/hub/checkpoints/esm2_t6_8M_UR50D.pt" | sha256sum -c - >/dev/null
printf '%s' "$R3_CKPT_TRANSFER_KEY" > .r3-private/key
openssl enc -d -aes-256-cbc -pbkdf2 -iter 600000 -in .r3-private/comparison.enc -out results/davis_r3_ckpt.pt -pass file:.r3-private/key
printf '%s  %s\n' "$SOURCE_PLAIN_SHA" 'results/davis_r3_ckpt.pt' | sha256sum -c - >/dev/null
python3 - <<'PY'
import json,os,torch
x=torch.load('results/davis_r3_ckpt.pt',map_location='cpu',weights_only=False)
chunk=int(os.environ['SOURCE_CHUNK'])
assert x['chunk']==chunk and x['batch_in_chunk']==0 and not x['done']
rows=[json.loads(s) for s in open('.r3-private/comparison-log.jsonl')]
assert len(rows)==chunk and rows[-1]['chunk']==chunk
assert abs(x['best_val']-rows[-1]['best_val_ci'])<0.0001
assert x['stall']==rows[-1]['stall']
PY
cp .r3-private/comparison-log.jsonl results/davis_r3_log.jsonl
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
import json,os,torch
rows=[json.loads(s) for s in open('results/davis_r3_log.jsonl')]
assert len(rows)==int(os.environ['EXPECTED_NEXT_CHUNK'])
x=torch.load('results/davis_r3_ckpt.pt',map_location='cpu',weights_only=False)
assert x['chunk']==len(rows) and x['batch_in_chunk']==0
# The trainer's locked stop and single test evaluation decide finality; do not override.
with open('.r3-private/comparison-status.json','w') as f:
    json.dump({'chunk':x['chunk'],'stall':x['stall'],'done':x['done'],'best_val':x['best_val']},f)
PY
# Replace the incoming ciphertext with the new checkpoint; artifact contains only ciphertext and metrics.
openssl enc -aes-256-cbc -salt -pbkdf2 -iter 600000 -in results/davis_r3_ckpt.pt -out .r3-private/next.enc -pass file:.r3-private/key
mv .r3-private/next.enc .r3-private/comparison.enc
sha256sum .r3-private/comparison.enc > .r3-private/comparison.sha256
sha256sum results/davis_r3_ckpt.pt > .r3-private/plain-checkpoint.sha256
python3 - <<'PY'
import hashlib,hmac,os
k=os.environ['R3_CKPT_TRANSFER_KEY'].encode()
p='.r3-private/comparison.enc'
with open('.r3-private/comparison.hmac','w') as f:
    f.write(hmac.new(k,open(p,'rb').read(),hashlib.sha256).hexdigest()+'\n')
PY
cp results/davis_r3_log.jsonl .r3-private/comparison-log.jsonl
rm -f .r3-private/key results/davis_r3_ckpt.pt
