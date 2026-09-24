"""ChEMBL public API client for bioactivity corroboration (free, no key).

Used to independently check model-predicted drug-kinase pairs: a pair is
'corroborated' when ChEMBL reports a binding assay (Ki/Kd/IC50) <= 1 uM
between the drug (by InChIKey) and the kinase (by gene name).
"""
from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DIR = os.environ.get(
    "TARGETSCAN_CACHE",
    os.path.join(os.path.dirname(__file__), "..", "..", "data_cache"),
)


class ChemblError(RuntimeError):
    pass


def _get(path: str, params: dict = None, timeout: int = 60):
    url = f"{API}{path}.json"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:
            if attempt == 2:
                raise ChemblError(f"{url}: {exc}") from exc
            time.sleep(2 ** attempt)


def _cached(name, fn):
    os.makedirs(CACHE_DIR, exist_ok=True)
    p = os.path.join(CACHE_DIR, name)
    if os.path.exists(p):
        return json.load(open(p))
    out = fn()
    json.dump(out, open(p, "w"))
    return out


def molecule_by_inchikey(inchikey: str):
    def _f():
        r = _get("/molecule", {"pref_name__isnull": "false",
                               "molecule_structures__standard_inchi_key": inchikey})
        return r.get("molecules", [])
    return _cached(f"mol_{inchikey}.json", _f)


def target_by_gene(gene: str, organism: str = "Homo sapiens"):
    def _f():
        r = _get("/target", {"target_components__target_component_synonyms__component_synonym": gene,
                             "target_type": "SINGLE PROTEIN"})
        return r.get("targets", [])
    return _cached(f"tgt_{gene}.json", _f)


def binding_activities(molecule_chembl_id: str, target_chembl_id: str,
                       limit: int = 50):
    """Assay activities between one molecule and one target."""
    def _f():
        r = _get("/activity", {
            "molecule_chembl_id": molecule_chembl_id,
            "target_chembl_id": target_chembl_id,
            "limit": limit,
        })
        acts = r.get("activities", [])
        return [{"type": a.get("standard_type"), "value": a.get("standard_value"),
                 "units": a.get("standard_units"), "assay_type": a.get("assay_type"),
                 "pchembl": a.get("pchembl_value")} for a in acts]
    return _cached(f"act_{molecule_chembl_id}_{target_chembl_id}.json", _f)


def corroborate(inchikey: str, gene: str, max_nm: float = 1000.0):
    """Return (status, evidence): corroborated / no_data / not_found."""
    mols = molecule_by_inchikey(inchikey)
    if not mols:
        return "not_found", "molecule not in ChEMBL"
    tgts = [t for t in target_by_gene(gene)
            if t.get("organism") == "Homo sapiens"] or target_by_gene(gene)
    if not tgts:
        return "not_found", "target not in ChEMBL"
    best = None
    for t in tgts[:3]:
        for a in binding_activities(mols[0]["molecule_chembl_id"],
                                    t["target_chembl_id"]):
            if a["type"] in ("Ki", "Kd", "IC50", "EC50") and a["value"]:
                try:
                    v = float(a["value"])
                except (TypeError, ValueError):
                    continue
                if a["units"] == "nM" and (best is None or v < best):
                    best = v
    if best is None:
        return "no_data", "no binding assay in ChEMBL"
    verdict = "corroborated" if best <= max_nm else "refuted"
    return verdict, f"best binding {best:.1f} nM"
