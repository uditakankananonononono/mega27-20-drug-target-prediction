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
