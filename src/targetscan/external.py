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


def kegg_pathways(symbol):
    """KEGG REST: pathways for a human gene symbol (resolves hsa id first)."""
    import urllib.request as u
    def _txt(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=20) as r:
            return r.read().decode()
    hits = [l for l in _txt(f"https://rest.kegg.jp/find/genes/{symbol}").strip().split("\n")
            if l.split("\t")[0].startswith("hsa:")]
    # exact symbol match first (field 2 is "SYM, ALIAS, ...; description")
    exact = [l for l in hits if l.split("\t")[1].split(";")[0].split(",")[0].strip() == symbol]
    pick = (exact or hits or [None])[0]
    if pick is None:
        return {"kegg_id": None, "pathways": []}
    kegg_id = pick.split("\t")[0]
    txt = _txt(f"https://rest.kegg.jp/link/pathway/{kegg_id}")
    return {"kegg_id": kegg_id,
            "pathways": [l.split("\t")[1] for l in txt.strip().split("\n") if "\t" in l]}


def rcsb_search(query):
    """RCSB PDB search API (GET): structures matching a free-text query."""
    import urllib.parse
    q = {"query": {"type": "terminal", "service": "full_text",
                   "parameters": {"value": query}},
         "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 5}}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(q))
    import urllib.request as u
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=25) as r:
        d = json.loads(r.read().decode())
    return [e["identifier"] for e in d.get("result_set", [])]


def hgnc_symbol(symbol):
    """HGNC REST: approved symbol/name for a gene."""
    url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "targetscan/0.1"})
    with urllib.request.urlopen(req, timeout=20) as r:
        d = json.loads(r.read().decode())
    docs = d.get("response", {}).get("docs", [])
    return docs[0] if docs else {}


def chebi_entry(chebi_id):
    """ChEBI record via EBI OLS4 REST (official ChEBI ontology route)."""
    import urllib.request as u
    url = f"https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?obo_id={chebi_id}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=20) as r:
        d = json.loads(r.read().decode())
    terms = d.get("_embedded", {}).get("terms", [])
    if not terms:
        raise KeyError(chebi_id)
    t0 = terms[0]
    return {"chebi_id": chebi_id, "label": t0.get("label"),
            "description": (t0.get("description") or [""])[0][:300]}


def reactome_pathways(uniprot_acc):
    """Reactome Content Service: pathways for a UniProt accession."""
    url = f"https://reactome.org/ContentService/data/mapping/UniProt/{uniprot_acc}/pathways"
    import urllib.request as u
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=20) as r:
        return json.loads(r.read().decode())


def opentargets_associations(ensembl_id, size=5):
    """Open Targets Platform GraphQL: top disease associations for a target."""
    q = {"query": '{ target(ensemblId: "%s") { approvedSymbol associatedDiseases(page: {index: 0, size: %d}) { rows { disease { name } score } } } }' % (ensembl_id, size)}
    import urllib.request as u
    req = u.Request("https://api.platform.opentargets.org/api/v4/graphql",
                    data=json.dumps(q).encode(),
                    headers={"Content-Type": "application/json", "User-Agent": "targetscan/0.1"})
    with u.urlopen(req, timeout=25) as r:
        d = json.loads(r.read().decode())
    t = d.get("data", {}).get("target") or {}
    rows = (t.get("associatedDiseases") or {}).get("rows", [])
    return [{"disease": r["disease"]["name"], "score": r["score"]} for r in rows]


def dgidb_interactions(gene, max_rows=10):
    """DGIdb GraphQL: known drug interactions for a gene."""
    import urllib.request as u
    q = {"query": '{ genes(names: ["%s"]) { nodes { name interactions { drug { name } interactionScore } } } }' % gene}
    req = u.Request("https://dgidb.org/api/graphql", data=json.dumps(q).encode(),
                    headers={"Content-Type": "application/json", "User-Agent": "targetscan/0.1"})
    with u.urlopen(req, timeout=30) as r:
        d = json.loads(r.read().decode())
    nodes = d["data"]["genes"]["nodes"]
    if not nodes:
        return {"gene": gene, "interactions": []}
    rows = [{"drug": i["drug"]["name"], "score": i["interactionScore"]}
            for i in nodes[0]["interactions"][:max_rows]]
    return {"gene": gene, "interactions": rows}


def pdbe_entry(pdb_id):
    """PDBe API: entry summary + entity inventory for a PDB id."""
    import urllib.request as u
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=20) as r:
            return json.loads(r.read().decode())
    pid = pdb_id.lower()
    s = _get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pid}")[pid][0]
    mols = _get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pid}").get(pid, [])
    return {"pdb_id": pid, "title": s["title"],
            "method": s.get("experimental_method", [None])[0],
            "entities": [{"type": m.get("molecule_type"),
                          "name": (m.get("molecule_name") or [""])[0]} for m in mols]}


def interpro_domains(uniprot_acc, page_size=25):
    """InterPro API: integrated domain/family entries for a UniProt protein."""
    import urllib.request as u
    url = (f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/"
           f"{uniprot_acc}/?page_size={page_size}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    return {"accession": uniprot_acc, "count": d["count"],
            "domains": [{"accession": x["metadata"]["accession"],
                         "name": x["metadata"]["name"],
                         "type": x["metadata"]["type"]} for x in d["results"]]}


def mobidb_entry(uniprot_acc):
    """MobiDB API: disorder + Pfam domains for a UniProt accession."""
    import urllib.request as u
    url = f"https://mobidb.org/api/download?acc={uniprot_acc}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    if not d:
        raise KeyError(uniprot_acc)
    m = d[0]
    dis = m.get("prediction-disorder-priority", {})
    pfam = m.get("homology-domain-pfam", {})
    return {"accession": uniprot_acc, "gene": m.get("gene"), "length": m.get("length"),
            "disorder_fraction": dis.get("content_fraction"),
            "pfam_domains": pfam.get("regions_names", [])[:4]}


def alphafold_prediction(uniprot_acc):
    """AlphaFold DB API: predicted-structure metadata for a UniProt accession."""
    import urllib.request as u
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    if not d:
        raise KeyError(uniprot_acc)
    m = d[0]
    return {"accession": uniprot_acc, "model": m["modelEntityId"],
            "gene": m.get("gene"), "mean_plddt": m["globalMetricValue"],
            "frac_very_high": m.get("fractionPlddtVeryHigh")}


def quickgo_annotations(uniprot_acc, aspect="biological_process", limit=100):
    """QuickGO (EBI GOA): GO annotations for a UniProt accession."""
    import urllib.request as u
    url = (f"https://www.ebi.ac.uk/QuickGO/services/annotation/search"
           f"?geneProductId=UniProtKB:{uniprot_acc}&goAspect={aspect}&limit={limit}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    return {"accession": uniprot_acc, "n_hits": d["numberOfHits"],
            "go_ids": sorted({x["goId"] for x in d["results"]})}


def monarch_gene_diseases(symbol):
    """Monarch Initiative v3 API: search + causal disease associations."""
    import urllib.request as u
    from urllib.parse import quote
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
            return json.loads(r.read().decode())
    s = _get(f"https://api.monarchinitiative.org/v3/api/search?q={quote(symbol)}&limit=3")
    hit = next((i for i in s["items"] if i["name"] == symbol and i["category"] == "biolink:Gene"), s["items"][0])
    a = _get("https://api.monarchinitiative.org/v3/api/association"
             f"?category=biolink:CausalGeneToDiseaseAssociation&entity={hit['id']}&limit=10")
    return {"symbol": symbol, "monarch_id": hit["id"], "xrefs": hit.get("xref", []),
            "causal_diseases": [x.get("object_label") or x.get("object") for x in a["items"]]}


def ncbi_gene(symbol, organism="human"):
    """NCBI E-utilities: GeneID + summary for a gene symbol."""
    import urllib.request as u
    from urllib.parse import quote
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
            return json.loads(r.read().decode())
    q = quote(f"{symbol}[sym] AND {organism}[orgn]")
    s = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={q}&retmode=json")
    ids = s["esearchresult"]["idlist"]
    if not ids:
        raise KeyError(symbol)
    summ = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={ids[0]}&retmode=json")
    doc = summ["result"][ids[0]]
    return {"symbol": symbol, "gene_id": ids[0], "name": doc.get("name"),
            "description": doc.get("description"), "chromosome": doc.get("chromosome")}


def ensembl_lookup(symbol):
    """Ensembl REST: gene id + biotype for a human symbol."""
    import urllib.request as u
    url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
           "?content-type=application/json")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        return json.loads(r.read().decode())


def pubchem_bioassays(gene_id, max_aids=200):
    """PubChem PUG REST: bioassay AIDs targeting a GeneID."""
    import urllib.request as u
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/"
           f"{gene_id}/aids/JSON?MaxRecords={max_aids}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    aids = d.get("IdentifierList", {}).get("AID", [])
    return {"gene_id": gene_id, "aids": aids, "n_aids": len(aids)}


def klifs_structures(symbol, species="Human"):
    """KLIFS API v2: PDB structure list for a kinase symbol (human)."""
    import urllib.request as u
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=60) as r:
            return json.loads(r.read().decode())
    names = _get("https://klifs.net/api_v2/kinase_names")
    hits = [k for k in names if k.get("gene_name") == symbol and k.get("species") == species]
    if not hits:
        raise KeyError(symbol)
    kin = hits[0]
    structs = _get(f"https://klifs.net/api_v2/structures_list?kinase_ID={kin['kinase_ID']}")
    return {"symbol": symbol, "kinase_id": kin["kinase_ID"], "uniprot": kin.get("accession"),
            "n_structures": len(structs),
            "pdbs": sorted({s.get("pdb") for s in structs if s.get("pdb")})}


def gprofiler_enrichment(symbols, sources=("GO:BP",), organism="hsapiens"):
    """g:Profiler g:GOSt functional enrichment for a gene list (POST)."""
    import urllib.request as u
    body = json.dumps({"organism": organism, "query": list(symbols),
                       "sources": list(sources), "no_evidences": True}).encode()
    req = u.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                    data=body, headers={"Content-Type": "application/json",
                                        "User-Agent": "targetscan/0.1"})
    with u.urlopen(req, timeout=60) as r:
        d = json.loads(r.read().decode())
    terms = [{"id": t.get("native"), "name": t.get("name"),
              "p_value": t.get("p_value"), "source": t.get("source")}
             for t in d.get("result", [])]
    return {"query": list(symbols), "n_terms": len(terms), "terms": terms}


def hpa_search(symbol, columns="g,gs,eg"):
    """Human Protein Atlas search_download API: gene search, gzipped JSON."""
    import gzip, io
    import urllib.request as u
    from urllib.parse import quote
    url = (f"https://www.proteinatlas.org/api/search_download.php"
           f"?search={quote(symbol)}&format=json&columns={columns}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        raw = r.read()
    try:
        raw = gzip.decompress(raw)
    except OSError:
        pass  # some responses are plain JSON
    rows = json.loads(raw.decode())
    hits = [x for x in rows if x.get("Gene") == symbol]
    return {"query": symbol, "n_rows": len(rows), "hits": hits}


def openfda_label(brand_name):
    """openFDA drug label API: indications for a brand-name drug."""
    import urllib.request as u
    from urllib.parse import quote
    url = (f"https://api.fda.gov/drug/label.json"
           f"?search=openfda.brand_name:%22{quote(brand_name)}%22&limit=1")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    res = d.get("results", [])
    if not res:
        return {"brand": brand_name, "found": False}
    r0 = res[0]
    return {"brand": brand_name, "found": True,
            "generic": (r0.get("openfda", {}).get("generic_name") or [None])[0],
            "indications": (r0.get("indications_and_usage") or [""])[0][:400],
            "set_id": r0.get("set_id")}


def omnipath_interactions(symbol, fields="sources"):
    """OmniPath REST: signed directed interactions for a gene symbol."""
    import urllib.request as u
    from urllib.parse import quote
    url = (f"https://omnipathdb.org/interactions?genesymbols=1"
           f"&partners={quote(symbol)}&fields={fields}&format=json")
    with u.urlopen(u.Request(url, headers={"User-Agent": "targetscan/0.1"}), timeout=30) as r:
        rows = json.loads(r.read().decode())
    partners = sorted({(x.get("source_genesymbol") if x.get("target_genesymbol") == symbol
                        else x.get("target_genesymbol")) for x in rows})
    dbs = sorted({s for x in rows for s in x.get("sources", [])})
    return {"symbol": symbol, "n_interactions": len(rows),
            "n_partners": len(partners), "partners_sample": partners[:10],
            "databases": dbs}


def wikidata_drug_targets(drug_label):
    """Wikidata SPARQL: drug QID + 'interacts with' (P129) target labels."""
    import urllib.request as u
    from urllib.parse import quote
    q = (
        "SELECT ?drug ?drugLabel ?target ?targetLabel WHERE {"
        f" ?drug rdfs:label \"{drug_label}\"@en ."
        " OPTIONAL { ?drug wdt:P129 ?target . }"
        " SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. }"
        " } LIMIT 50")
    url = ("https://query.wikidata.org/sparql?query=" + quote(q)
           + "&format=json")
    req = u.Request(url, headers={"Accept": "application/sparql-results+json",
                                  "User-Agent": "targetscan/0.1 (research)"})
    with u.urlopen(req, timeout=40) as r:
        d = json.loads(r.read().decode())
    rows = d.get("results", {}).get("bindings", [])
    qid = rows[0]["drug"]["value"].rsplit("/", 1)[-1] if rows else None
    targets = sorted({x["targetLabel"]["value"] for x in rows if "targetLabel" in x})
    return {"drug": drug_label, "qid": qid, "n_target_rows": len(rows),
            "interacts_with": targets}


def ebi_proteins_variants(accession, positions=None, size=400):
    """EBI Proteins variation API: variant features for a UniProt accession,
    optionally filtered to 1-based positions (string compare)."""
    import urllib.request as u
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}?size={size}"
    with u.urlopen(u.Request(url, headers={"Accept": "application/json",
                                           "User-Agent": "targetscan/0.1"}), timeout=40) as r:
        d = json.loads(r.read().decode())
    feats = d.get("features", [])
    out = []
    for f in feats:
        if positions and f.get("begin") not in positions:
            continue
        clin = f.get("clinicalSignificances") or []
        out.append({"pos": f.get("begin"), "wt": f.get("wildType"),
                    "mut": f.get("mutatedType"), "type": f.get("type"),
                    "consequence": f.get("consequenceType"),
                    "clinical": [c.get("type") for c in clin]})
    return {"accession": accession, "n_features_total": len(feats),
            "variants": out}
