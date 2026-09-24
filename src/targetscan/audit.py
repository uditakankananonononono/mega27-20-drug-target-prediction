"""DTI dataset auditor - name/sequence/label consistency checks.

Built from the DAVIS mutation-erasure finding: a dataset may label targets
as mutants while shipping WT sequences, making mutation effects unlearnable.
This module flags such defects in any {name: sequence} + affinity-matrix
dataset. Pure functions over plain data; no model, no network.
"""
import hashlib


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()


def find_erased_variants(sequences, wt_of=None):
    """Entries whose name marks a variant (e.g. 'KIT(D816V)') but whose
    sequence is byte-identical to the named wild type.

    sequences: dict {name: aa_sequence}
    wt_of: optional callable name -> wt name; default strips '(...)'.
    Returns list of {mutant, wt, seq_md5}.
    """
    if wt_of is None:
        wt_of = lambda n: n.split("(")[0] if "(" in n else None
    out = []
    for n, s in sequences.items():
        wt = wt_of(n)
        if wt and wt in sequences and sequences[wt] == s:
            out.append({"mutant": n, "wt": wt, "seq_md5": md5(s)})
    return out


def label_variance_report(sequences, affinity, names, threshold=1.0):
    """For each erased variant, quantify the label spread vs WT that no
    sequence-based model can represent.

    affinity: (n_drugs, n_targets) array-like; names: target order.
    Returns dict with totals and per-target stats.
    """
    import numpy as np
    Y = np.asarray(affinity, dtype=float)
    idx = {n: i for i, n in enumerate(names)}
    erased = find_erased_variants({n: sequences[n] for n in names if n in sequences})
    per, tot, big = [], 0, 0
    for e in erased:
        d = np.abs(Y[:, idx[e["wt"]]] - Y[:, idx[e["mutant"]]])
        b = int((d >= threshold).sum())
        tot += len(d)
        big += b
        per.append({"mutant": e["mutant"], "wt": e["wt"],
                    "mean_abs_delta": round(float(d.mean()), 4),
                    "n_abs_delta_ge_threshold": b,
                    "max_abs_delta": round(float(d.max()), 3)})
    per.sort(key=lambda r: -r["n_abs_delta_ge_threshold"])
    return {"n_mutant_entries": len(erased), "n_pairs": tot,
            "n_pairs_ge_threshold": big,
            "pct_pairs_ge_threshold": round(100.0 * big / tot, 2) if tot else 0.0,
            "threshold": threshold, "per_target": per}


def audit(sequences, affinity=None, names=None, threshold=1.0):
    """Full audit: erased variants + (optionally) label variance."""
    rep = {"erased_variants": find_erased_variants(sequences)}
    if affinity is not None and names is not None:
        rep["label_variance"] = label_variance_report(sequences, affinity, names, threshold)
    rep["verdict"] = ("FAIL: mutant labels with WT sequences - mutation effects unlearnable"
                      if rep["erased_variants"] else "PASS: no name/sequence inconsistencies")
    return rep
