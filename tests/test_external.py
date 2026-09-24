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
