#!/usr/bin/env python3
"""R2 groundwork: correct the DAVIS mutant-sequence defect found by our audit.

DAVIS-as-distributed (DeepDTA mirror) ships all 54 mutant kinase entries with
sequences byte-identical to their wild types while labels differ (up to 4.49
pKd). This script builds corrected mutant sequences by applying the point
mutation encoded in each entry name (e.g. ABL1(T315I)) to the in-repo DAVIS
wild-type sequence, VALIDATED against the UniProt canonical sequence:

  1. direct:   davis_wt[pos-1] == original residue -> apply at pos
  2. offset:   residue matches in UniProt canonical at pos, and a unique 21-mer
               anchor aligns davis_wt inside canonical with a constant offset ->
               apply at pos-offset in davis_wt
  3. review:   neither holds -> entry goes to manual_review, sequence NOT made up

Entries suffixed 'p' (phospho state) reuse their base mutation's sequence
(phosphorylation is not sequence-encoded); recorded in provenance.

Outputs:
  data_cache/davis_proteins_corrected.json  (all 442 targets, mutants fixed)
  results/uniprot_variant_pull_provenance.json

Usage: python3 src/pull_uniprot_variants.py
"""
import json, re, time, urllib.request, urllib.parse, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data_cache/davis_proteins.txt"
OUT = ROOT / "data_cache/davis_proteins_corrected.json"
PROV = ROOT / "results/uniprot_variant_pull_provenance.json"
MUT_RE = re.compile(r"^([A-Za-z0-9]+)\(([A-Z])(\d+)([A-Z])\)(p?)$")

def fetch_canonical(gene):
    q = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true")
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=fasta&size=1"
    with urllib.request.urlopen(url, timeout=20) as r:
        fasta = r.read().decode()
    lines = fasta.strip().split("\n")
    acc = lines[0].split("|")[1] if "|" in lines[0] else lines[0]
    return acc, "".join(lines[1:])

def kmer_offset(anchor_seq, canonical, k=21):
    """Dominant constant offset of anchor_seq inside canonical, by k-mer vote.

    Constructs may extend or internally differ from canonical, so a single
    unique-anchor rule is too strict: take k-mers across the sequence, keep
    offsets with >=3 agreeing votes, require the top offset to hold >=70% of
    all votes. Returns (offset, n_votes_for_top) or (None, 0)."""
    from collections import Counter
    votes = Counter()
    for i in range(0, len(anchor_seq) - k + 1):
        j = canonical.find(anchor_seq[i:i + k])
        if j >= 0:
            votes[i - j] += 1
    if not votes:
        return None, 0
    off, top = votes.most_common(1)[0]
    if top >= 3 and top >= 0.7 * sum(votes.values()):
        return off, top
    return None, 0

def main():
    prots = json.loads(SRC.read_text())
    prov = {"script": "src/pull_uniprot_variants.py", "source": str(SRC.name),
            "uniprot_api": "rest.uniprot.org/uniprotkb (reviewed, human canonical)",
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "entries": [], "manual_review": [], "counts": {}}
    corrected, cache = dict(prots), {}
    mutants = [(k, *MUT_RE.match(k).groups()) for k in prots if MUT_RE.match(k)]
    counts = {"direct": 0, "offset": 0, "p_suffix_reuse": 0, "manual_review": 0}
    for name, wt, orig, pos_s, new, p in mutants:
        pos = int(pos_s)
        entry = {"mutant": name, "wt": wt, "mutation": f"{orig}{pos}{new}"}
        if p:
            base_name = name[:-1]
            if base_name in corrected and base_name != name and corrected[base_name] != prots[base_name]:
                corrected[name] = corrected[base_name]
                entry.update(method="p_suffix_reuse", note="phospho state not sequence-encoded; reused " + base_name)
                counts["p_suffix_reuse"] += 1
                prov["entries"].append(entry); continue
        davis_wt = prots.get(wt)
        if davis_wt is None:
            entry.update(method="manual_review", reason="no WT entry in davis_proteins")
            counts["manual_review"] += 1; prov["entries"].append(entry); prov["manual_review"].append(name); continue
        if pos <= len(davis_wt) and davis_wt[pos - 1] == orig:
            corrected[name] = davis_wt[:pos - 1] + new + davis_wt[pos:]
            entry.update(method="direct", note="davis_wt residue matches at stated position")
            counts["direct"] += 1; prov["entries"].append(entry); continue
        if wt not in cache:
            cache[wt] = fetch_canonical(wt); time.sleep(0.34)  # <=3 req/s
        acc, canon = cache[wt]
        entry["uniprot_accession"] = acc
        canon_ok = pos <= len(canon) and canon[pos - 1] == orig
        off, nvotes = kmer_offset(davis_wt, canon) if canon_ok else (None, 0)
        applied = pos + off
        if canon_ok and off is not None and 1 <= applied <= len(davis_wt) and davis_wt[applied - 1] == orig:
            p2 = applied
            corrected[name] = davis_wt[:p2 - 1] + new + davis_wt[p2:]
            entry.update(method="offset", offset=off, applied_position=p2, anchor_votes=nvotes,
                         note="UniProt canonical matches at stated position; offset via dominant 21-mer vote")
            counts["offset"] += 1; prov["entries"].append(entry); continue
        entry.update(method="manual_review",
                     reason=f"residue mismatch (davis_wt[{pos}]={davis_wt[pos-1] if pos<=len(davis_wt) else '?'}, canonical_ok={canon_ok}, anchor_offset={off}, votes={nvotes})")
        counts["manual_review"] += 1; prov["entries"].append(entry); prov["manual_review"].append(name)
    prov["counts"] = counts
    prov["n_targets_out"] = len(corrected)
    prov["md5"] = hashlib.md5(json.dumps(corrected, sort_keys=True).encode()).hexdigest()
    OUT.write_text(json.dumps(corrected))
    PROV.write_text(json.dumps(prov, indent=1))
    print("counts:", counts, "| manual_review:", prov["manual_review"], "| out:", OUT.name)

if __name__ == "__main__":
    main()
