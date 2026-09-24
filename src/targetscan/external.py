"""External research-tool clients for verification and enrichment.
Each client hits one public resource and caches raw JSON to disk.
Live calls only outside CI (hermetic tests mock _get)."""
import json, os, time, urllib.parse, urllib.request

CACHE = os.path.join(os.path.dirname(__file__), "..", "..", "data_cache", "external")


def _get(url, timeout=20):
    req = urllib.request.Request(url, headers={"User-Agent": "targetscan/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _cached(name, fetcher):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    data = fetcher()
    json.dump(data, open(path, "w"))
    time.sleep(0.2)  # be polite to public APIs
    return data


def pubchem_props(smiles):
    """PubChem PUG-REST: physchem properties for a SMILES."""
    enc = urllib.parse.quote(smiles, safe="")
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
           f"{enc}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,"
           "HBondAcceptorCount,RotatableBondCount/JSON")
    d = _get(url)
    return d["PropertyTable"]["Properties"][0]


def uniprot_entry(gene, organism="9606"):
    """UniProt REST: canonical human entry for a gene name."""
    q = urllib.parse.quote(f"gene:{gene} AND organism_id:{organism} AND reviewed:true")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           "&fields=accession,id,protein_name,length,cc_function&format=json&size=1")
    d = _get(url)
    return d.get("results", [{}])[0]


def string_interactions(genes, species=9606, limit=10):
    """STRING API: PPI network among given gene symbols."""
    ids = "%0d".join(genes)
    url = (f"https://string-db.org/api/json/network?identifiers={ids}"
           f"&species={species}&limit={limit}")
    return _get(url)


def clinicaltrials_search(query, page_size=5):
    """ClinicalTrials.gov API v2: trials matching a query term."""
    q = urllib.parse.quote(query)
    url = (f"https://clinicaltrials.gov/api/v2/studies?query.term={q}"
           f"&pageSize={page_size}&fields=NCTId,BriefTitle,OverallStatus")
    return _get(url).get("studies", [])


def europepmc_search(query, page_size=5):
    """Europe PMC REST: literature hits for a query."""
    q = urllib.parse.quote(query)
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query={q}&format=json&pageSize={page_size}")
    return _get(url).get("resultList", {}).get("result", [])
