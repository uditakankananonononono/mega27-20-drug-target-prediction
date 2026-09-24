"""Second live pull: KEGG, RCSB PDB, HGNC, Reactome for the hit kinases."""
import json, sys
sys.path.insert(0, "src")
from targetscan import external

HITS = {"PLK1": "P53350", "KIT": "P10721", "FLT3": "P36888", "ABL1": "P00519",
        "EGFR": "P00533"}
out = {"kegg": {}, "pdb": {}, "hgnc": {}, "reactome": {}}
for g, acc in HITS.items():
    try: out["kegg"][g] = external._cached("kegg_" + g, lambda: external.kegg_pathways(g))
    except Exception as e: out["kegg"][g] = {"error": str(e)[:80]}
    try: out["pdb"][g] = external._cached("pdb_" + g, lambda: external.rcsb_search(g))
    except Exception as e: out["pdb"][g] = {"error": str(e)[:80]}
    try: out["hgnc"][g] = external._cached("hgnc_" + g, lambda: external.hgnc_symbol(g))
    except Exception as e: out["hgnc"][g] = {"error": str(e)[:80]}
    try: out["reactome"][g] = external._cached("react_" + g, lambda: external.reactome_pathways(acc))
    except Exception as e: out["reactome"][g] = {"error": str(e)[:80]}
json.dump(out, open("results/external_pull2.json", "w"), default=str)
for k in out:
    ok = sum(1 for v in out[k].values() if not (isinstance(v, dict) and "error" in v))
    print(k, "ok:", ok, "/", len(HITS))
print("PLK1 kegg:", out["kegg"]["PLK1"])
print("PLK1 pdb:", out["pdb"]["PLK1"])
