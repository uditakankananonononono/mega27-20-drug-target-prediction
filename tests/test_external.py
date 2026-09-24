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
