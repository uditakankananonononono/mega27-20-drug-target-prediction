"""Hermetic tests for external clients (all HTTP mocked)."""
import json
from unittest import mock
from targetscan import external


def fake_get(payload):
    return lambda url, timeout=20: payload


def test_pubchem_props_parses(monkeypatch):
    payload = {"PropertyTable": {"Properties": [{"MolecularWeight": 521.7, "XLogP": 2.8}]}}
    monkeypatch.setattr(external, "_get", fake_get(payload))
    d = external.pubchem_props("CCO")
    assert d["MolecularWeight"] == 521.7


def test_uniprot_entry_first_result(monkeypatch):
    monkeypatch.setattr(external, "_get",
                        fake_get({"results": [{"primaryAccession": "P53350"}]}))
    assert external.uniprot_entry("PLK1")["primaryAccession"] == "P53350"


def test_string_network_edges(monkeypatch):
    edges = [{"preferredName_A": "PLK1", "preferredName_B": "CDK1", "score": 0.99}]
    monkeypatch.setattr(external, "_get", fake_get(edges))
    out = external.string_interactions(["PLK1", "CDK1"])
    assert out[0]["score"] == 0.99


def test_clinicaltrials_fields(monkeypatch):
    monkeypatch.setattr(external, "_get",
                        fake_get({"studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT1"}}}]}))
    assert len(external.clinicaltrials_search("volasertib")) == 1


def test_europepmc_results(monkeypatch):
    monkeypatch.setattr(external, "_get",
                        fake_get({"resultList": {"result": [{"id": "1", "title": "t"}]}}))
    assert external.europepmc_search("BI-2536 PLK1")[0]["id"] == "1"


def test_cache_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(external, "CACHE", str(tmp_path))
    calls = {"n": 0}
    def fetch():
        calls["n"] += 1
        return {"x": 1}
    external._cached("t1", fetch)
    external._cached("t1", fetch)
    assert calls["n"] == 1
    assert json.load(open(tmp_path / "t1.json"))["x"] == 1


def test_cli_eval_reports_committed_log(tmp_path, capsys):
    import json
    from targetscan.cli import main
    log = tmp_path / "davis_log.jsonl"
    log.write_text(json.dumps({"epoch": 45, "test_ci": 0.7048, "test_mse": 0.7377}) + "\n")
    main(["--ckpt", str(tmp_path / "x.pt"), "eval"])
    out = json.loads(capsys.readouterr().out)
    assert out["test_ci"] == 0.7048 and out["source_file"].endswith("davis_log.jsonl")


def test_kegg_and_hgnc(monkeypatch):
    import targetscan.external as ex
    class R:
        def __init__(self, body): self.body = body
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return self.body
    calls = {"n": 0}
    def fake_urlopen(req, timeout=20):
        calls["n"] += 1
        url = req.full_url if hasattr(req, "full_url") else req
        if "find/genes" in url:
            return R(b"hsa:5347\tPLK1; serine/threonine-protein kinase PLK1\n")
        return R(b"hsa:5347\tpath:hsa04110\n")
    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    out = ex.kegg_pathways("PLK1")
    assert out["kegg_id"] == "hsa:5347" and out["pathways"] == ["path:hsa04110"]


def test_rcsb_search(monkeypatch):
    import targetscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b'{"result_set": [{"identifier": "4J52"}]}'
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=25: R())
    assert ex.rcsb_search("PLK1") == ["4J52"]


def test_opentargets(monkeypatch):
    import targetscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b'{"data": {"target": {"approvedSymbol": "PLK1", "associatedDiseases": {"rows": [{"disease": {"name": "X"}, "score": 0.5}]}}}}'
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=25: R())
    out = ex.opentargets_associations("ENSG00000166851")
    assert out == [{"disease": "X", "score": 0.5}]


def test_chebi_entry(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"_embedded": {"terms": [{"label": "imatinib", "description": ["A benzamide."]}]}}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=20: R())
    out = ex.chebi_entry("CHEBI:45783")
    assert out["label"] == "imatinib"


def test_dgidb_interactions(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"data": {"genes": {"nodes": [{"name": "FLT3", "interactions": [
                {"drug": {"name": "MIDOSTAURIN"}, "interactionScore": 0.9}]}]}}}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.dgidb_interactions("FLT3")
    assert out["interactions"][0]["drug"] == "MIDOSTAURIN"


def test_pdbe_entry(monkeypatch):
    import targetscan.external as ex, json as j
    calls = iter([
        {"3d5u": [{"title": "Plk1 catalytic domain", "experimental_method": ["X-ray diffraction"]}]},
        {"3d5u": [{"molecule_type": "polypeptide(L)", "molecule_name": ["PLK"]}]},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=20: R(next(calls)))
    out = ex.pdbe_entry("3D5U")
    assert out["entities"][0]["name"] == "PLK"


def test_interpro_domains(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"count": 1, "results": [{"metadata": {
                "accession": "IPR000719", "name": "Protein kinase domain", "type": "domain"}}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.interpro_domains("P10721")
    assert out["domains"][0]["accession"] == "IPR000719"


def test_mobidb_entry(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps([{"acc": "P10721", "gene": "KIT", "length": 976,
                             "prediction-disorder-priority": {"content_fraction": 0.06},
                             "homology-domain-pfam": {"regions_names": ["Ig-like", "Kinase"]}}]).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.mobidb_entry("P10721")
    assert out["gene"] == "KIT" and out["pfam_domains"] == ["Ig-like", "Kinase"]


def test_alphafold_prediction(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps([{"modelEntityId": "AF-P10721-F1", "gene": "KIT",
                             "globalMetricValue": 80.0, "fractionPlddtVeryHigh": 0.5}]).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.alphafold_prediction("P10721")
    assert out["gene"] == "KIT"


def test_quickgo_annotations(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"numberOfHits": 1, "results": [{"goId": "GO:0000086"}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.quickgo_annotations("P53350")
    assert out["go_ids"] == ["GO:0000086"]


def test_monarch_gene_diseases(monkeypatch):
    import targetscan.external as ex, json as j
    calls = iter([
        {"items": [{"id": "HGNC:6186", "name": "KIT", "category": "biolink:Gene", "xref": ["OMIM:164920"]}]},
        {"items": [{"object_label": "gastrointestinal stromal tumor"}]},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(next(calls)))
    out = ex.monarch_gene_diseases("KIT")
    assert out["causal_diseases"] == ["gastrointestinal stromal tumor"]


def test_ncbi_gene(monkeypatch):
    import targetscan.external as ex, json as j
    calls = iter([
        {"esearchresult": {"idlist": ["3815"]}},
        {"result": {"3815": {"name": "KIT", "description": "KIT proto-oncogene", "chromosome": "4"}}},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(next(calls)))
    out = ex.ncbi_gene("KIT")
    assert out["gene_id"] == "3815" and out["chromosome"] == "4"


def test_ensembl_lookup(monkeypatch):
    import targetscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"id": "ENSG00000157404", "biotype": "protein_coding", "display_name": "KIT"}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.ensembl_lookup("KIT")
    assert out["id"] == "ENSG00000157404"


def test_pubchem_bioassays_parses_aid_list(monkeypatch):
    from targetscan.external import pubchem_bioassays
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"IdentifierList": {"AID": [101, 202, 303]}}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = pubchem_bioassays(5347)
    assert out["gene_id"] == 5347
    assert out["aids"] == [101, 202, 303]
    assert out["n_aids"] == 3


def test_klifs_structures_parses(monkeypatch):
    from targetscan.external import klifs_structures
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payloads = iter([
        J.dumps([{"kinase_ID": 311, "gene_name": "PLK1", "species": "Human",
                  "accession": "P53350"}]).encode(),
        J.dumps([{"pdb": "3D5U"}, {"pdb": "4J52"}, {"pdb": "3D5U"}]).encode(),
    ])
    monkeypatch.setattr("urllib.request.urlopen",
                        lambda req, timeout=60: R(next(payloads)))
    out = klifs_structures("PLK1")
    assert out["kinase_id"] == 311
    assert out["n_structures"] == 3
    assert out["pdbs"] == ["3D5U", "4J52"]


def test_audit_detects_erased_variants():
    from targetscan.audit import audit, find_erased_variants
    seqs = {"KIT": "AAAA", "KIT(D816V)": "AAAA", "EGFR": "MMMM", "EGFR(L858R)": "MMMN"}
    erased = find_erased_variants(seqs)
    assert [e["mutant"] for e in erased] == ["KIT(D816V)"]
    import numpy as np
    aff = np.array([[1.0, 1.0, 2.0, 3.0], [1.0, 3.5, 2.0, 3.2]])
    rep = audit(seqs, aff, ["KIT", "KIT(D816V)", "EGFR", "EGFR(L858R)"], threshold=1.0)
    assert rep["label_variance"]["n_mutant_entries"] == 1
    assert rep["label_variance"]["n_pairs_ge_threshold"] == 1
    assert rep["verdict"].startswith("FAIL")
    clean = audit({"A": "AAAA", "B": "CCCC"})
    assert clean["verdict"].startswith("PASS")


def test_gprofiler_enrichment_parses(monkeypatch):
    from targetscan.external import gprofiler_enrichment
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"result": [
        {"native": "GO:0000278", "name": "mitotic cell cycle",
         "p_value": 1e-9, "source": "GO:BP"}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=60: R(payload))
    out = gprofiler_enrichment(["PLK1", "CDK1"])
    assert out["n_terms"] == 1
    assert out["terms"][0]["id"] == "GO:0000278"


def test_hpa_search_handles_gzip_and_exact_match(monkeypatch):
    from targetscan.external import hpa_search
    import gzip, io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = gzip.compress(J.dumps([
        {"Gene": "PLK1", "Ensembl": "ENSG00000166851"},
        {"Gene": "KIZ", "Ensembl": "ENSG00000088970"}]).encode())
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = hpa_search("PLK1")
    assert out["n_rows"] == 2
    assert out["hits"] == [{"Gene": "PLK1", "Ensembl": "ENSG00000166851"}]
