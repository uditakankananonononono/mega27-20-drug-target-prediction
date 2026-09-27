"""Build the frozen Pfam PF00069 (Pkinase) domain crop table for the 442 DAVIS
proteins (PREREG_DAVIS_R3.md section R3.1). Sequence-only evidence: UniProt
reviewed human entry per base gene -> 'Protein kinase' domain feature
start/end. Names that do not resolve -> keep full length, declared.
Output: data_cache/davis_pfam_crop.json  {davis_name: {base, start, end, source, uniprot}}
1-based inclusive boundaries, mapped onto the DAVIS sequence; if the DAVIS
sequence is not found to contain the UniProt canonical subsequence, falls
back to direct coordinates only when len matches, else full length + flag.
"""
import json, re, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent.parent
OUT = ROOT / "data_cache" / "davis_pfam_crop.json"
UP = "https://rest.uniprot.org/uniprotkb/search"


def fetch_domain(gene):
    q = (f"{UP}?query=gene_exact:{gene}+AND+organism_id:9606+AND+reviewed:true"
         f"&fields=accession,gene_names,protein_name,ft_domain,sequence&format=json&size=1")
    with urllib.request.urlopen(q, timeout=30) as r:
        d = json.load(r)
    if not d.get("results"):
        return None
    e = d["results"][0]
    acc = e["primaryAccession"]
    seq = e["sequence"]["value"]
    doms = [f for f in e.get("features", []) if f.get("type") == "Domain"
            and "kinase" in f.get("description", "").lower()]
    if not doms:
        return {"uniprot": acc, "seq": seq, "domain": None}
    f = max(doms, key=lambda x: x["location"]["end"]["value"] - x["location"]["start"]["value"])
    return {"uniprot": acc, "seq": seq,
            "domain": (f["location"]["start"]["value"], f["location"]["end"]["value"], f["description"])}


def base_of(name):
    b = re.sub(r"\(.*?\)", "", name)
    return b[:-1] if b.endswith("p") else b


def main():
    prots = json.load(open(ROOT / "data_cache" / "davis_proteins_corrected.json"))
    table = json.load(open(OUT)) if OUT.exists() else {}
    genes = sorted({base_of(k) for k in prots} - {v["base"] for v in table.values()})
    cache_p = ROOT / "data_cache" / "uniprot_domain_cache.json"
    cache = json.load(open(cache_p)) if cache_p.exists() else {}
    todo = [g for g in genes if g not in cache]
    print(f"to fetch: {len(todo)}", flush=True)
    for g in todo:
        try:
            r = fetch_domain(g)
        except Exception as ex:
            r = {"error": str(ex)[:80]}
        cache[g] = r
        time.sleep(0.15)
        if len(cache) % 25 == 0:
            json.dump(cache, open(cache_p, "w"))
            print(f"cached {len(cache)}", flush=True)
    json.dump(cache, open(cache_p, "w"))
    # assemble
    for name, seq in prots.items():
        if name in table:
            continue
        b = base_of(name)
        c = cache.get(b)
        rec = {"base": b, "start": 1, "end": len(seq), "source": "full_length", "uniprot": None}
        if c and not c.get("error") and c.get("domain"):
            us, ue, desc = c["domain"]
            useq = c["seq"]
            if useq == seq or useq in seq:
                off = seq.find(useq)
                rec = {"base": b, "start": off + us, "end": off + ue,
                       "source": f"uniprot:{desc}", "uniprot": c["uniprot"]}
            elif len(useq) == len(seq):
                rec = {"base": b, "start": us, "end": ue,
                       "source": f"uniprot-lenmatch:{desc}", "uniprot": c["uniprot"]}
            else:
                rec["source"] = f"full_length:seq_mismatch({c['uniprot']})"
        table[name] = rec
    json.dump(table, open(OUT, "w"), indent=1)
    n_crop = sum(1 for v in table.values() if v["source"].startswith("uniprot"))
    print(f"table: {len(table)} proteins, {n_crop} cropped, {len(table)-n_crop} full-length", flush=True)


if __name__ == "__main__":
    main()
