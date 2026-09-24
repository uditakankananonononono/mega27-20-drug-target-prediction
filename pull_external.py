"""Live pulls from external research tools -> data_cache/external/.
Run outside CI. Bounded: DAVIS drugs (PubChem), top kinases (UniProt),
BI-2536 literature + trials, STRING network of hit targets."""
import json, sys, time
sys.path.insert(0, "src")
from targetscan.data.davis_kiba import load_davis
from targetscan import external

ds = load_davis()
out = {"pubchem": {}, "uniprot": {}, "trials": None, "literature": None, "string": None}
t0 = time.time()
for i, smi in enumerate(ds.smiles):
    try:
        out["pubchem"][smi] = external._cached("pc_%d" % i, lambda s=smi: external.pubchem_props(s))
    except Exception as e:
        out["pubchem"][smi] = {"error": str(e)[:120]}
    if time.time() - t0 > 45:
        print("time bound at drug", i); break
print("pubchem done:", len(out["pubchem"]))
for g in ["PLK1", "KIT", "FLT3", "ABL1", "EGFR", "SRC", "CDK2", "MET", "BRAF", "JAK2"]:
    try:
        out["uniprot"][g] = external._cached("up_" + g, lambda: external.uniprot_entry(g))
    except Exception as e:
        out["uniprot"][g] = {"error": str(e)[:120]}
print("uniprot done:", len(out["uniprot"]))
out["trials"] = external._cached("trials_volasertib", lambda: external.clinicaltrials_search("volasertib"))
out["literature"] = external._cached("epmc_bi2536", lambda: external.europepmc_search('"BI-2536" AND PLK1'))
out["string"] = external._cached("string_hits", lambda: external.string_interactions(["PLK1", "KIT", "FLT3", "ABL1"]))
json.dump(out, open("results/external_pull.json", "w"), default=str)
print("OK", {k: len(v) if hasattr(v, "__len__") else v for k, v in out.items()})
