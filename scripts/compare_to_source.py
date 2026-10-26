#!/usr/bin/env python3
"""Compare NUCLEON with the IAEA LiveChart export it was built from.

For every row of source-data/IAEA nuclide livechart.csv (ground states), check
that NUCLEON has the nuclide and that its values match the source:

  existence, stability flag, half-life in seconds (:halfLife against
  half_life_sec, or :decayWidthValue/:decayWidthUnit for resonances given as a
  width), the original value and unit (:halfLifeValue / :halfLifeUnit),
  natural abundance, and decay modes with branching ratios.

This is the automated counterpart to the manual audit: it tests whether
anything went wrong between the source file and the ontology. Isomers are
not in the export and are covered by the manual audit sample instead.

Usage: python3 scripts/compare_to_source.py [NUCLEON.ttl] [csv] [outdir]
"""
import csv
import math
import os
import sys
from collections import Counter

import pandas as pd
import rdflib
from rdflib import RDF, RDFS, Namespace

N = Namespace("https://rsna.org/radlex/nucleon#")
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
WIDTH_UNITS = {"eV", "keV", "MeV"}
UNIT_ALIASES = {"Y": "a", "y": "a", "M": "m", "D": "d", "H": "h", "S": "s", "MS": "ms", "US": "µs", "us": "µs",
                "NS": "ns", "PS": "ps", "FS": "fs", "AS": "as", "ZS": "zs", "YS": "ys"}
REL_TOL = 1e-3


def local(u):
    return str(u).split("#", 1)[1]


def close(a, b, tol=REL_TOL):
    if a is None or b is None:
        return a is None and b is None
    if a == b:
        return True
    return abs(a - b) <= tol * max(abs(a), abs(b))


def main(ttl="NUCLEON.ttl", src="source-data/IAEA nuclide livechart.csv", outdir="validation"):
    g = rdflib.Graph()
    g.parse(ttl)
    code_to_mode = {}
    for m in g.subjects(RDF.type, N.DecayMode):
        for a in g.objects(m, SKOS.altLabel):
            if str(a).isupper() or any(c.isdigit() for c in str(a)) or str(a) in ("A", "N", "P"):
                code_to_mode.setdefault(str(a), local(m))
    code_to_mode["SF+EC+B+"] = "BetaPlusDecay"   # bundled rate, documented in branch comments

    d = pd.read_csv(src)
    issues, n_rows = [], 0
    counts = Counter()

    def flag(name, field, nucleon, source, note=""):
        issues.append({"nuclide": name, "field": field, "nucleon": nucleon, "source": source, "note": note})
        counts[field] += 1

    for _, r in d.iterrows():
        n_rows += 1
        name = f"{r.symbol}-{int(r.z + r.n)}"
        u = N[name]
        if (u, N.massNumber, None) not in g:
            flag(name, "missing from NUCLEON", "", "present")
            continue
        stable_src = str(r.half_life).strip().upper() == "STABLE"
        stable_n = (u, N.isStable, rdflib.Literal(True)) in g
        if stable_src != stable_n:
            flag(name, "isStable", stable_n, stable_src)

        # half-life or decay width
        unit = None if pd.isna(r.unit_hl) else str(r.unit_hl).strip()
        if not stable_src and unit:
            hl_n = g.value(u, N.halfLife)
            if unit in WIDTH_UNITS:
                wv, wu = g.value(u, N.decayWidthValue), g.value(u, N.decayWidthUnit)
                if wv is None or not close(float(wv), float(r.half_life)) or str(wu) != unit:
                    flag(name, "decay width", f"{wv} {wu}", f"{r.half_life} {unit}")
            else:
                src_s = None if pd.isna(r.half_life_sec) else float(r.half_life_sec)
                if hl_n is None and src_s is not None:
                    flag(name, "halfLife (s)", "", src_s, "absent in NUCLEON")
                elif hl_n is not None and src_s is not None and not close(float(hl_n), src_s):
                    ratio = float(hl_n) / src_s
                    note = ""
                    for lab, f in (("h vs s", 3600), ("min vs s", 60), ("d vs s", 86400), ("a vs s", 31557600)):
                        if math.isclose(ratio * f, 1, rel_tol=5e-3) or math.isclose(ratio / f, 1, rel_tol=5e-3):
                            note = f"off by a factor of {f}: likely unit error ({lab})"
                    flag(name, "halfLife (s)", float(hl_n), src_s, note or f"ratio {ratio:.4g}")
                v, un = g.value(u, N.halfLifeValue), g.value(u, N.halfLifeUnit)
                if v is not None and (not close(float(v), float(r.half_life))
                                      or UNIT_ALIASES.get(unit, unit) != str(un)):
                    flag(name, "halfLifeValue/Unit", f"{v} {un}", f"{r.half_life} {unit}")

        # abundance
        if not pd.isna(r.abundance):
            ab = g.value(u, N.naturalAbundance)
            if ab is None or not close(float(ab), float(r.abundance)):
                flag(name, "naturalAbundance", ab, r.abundance)

        # decay modes and ratios
        src_modes = []
        for i in (1, 2, 3):
            code = r[f"decay_{i}"]
            if pd.isna(code):
                continue
            pct = r[f"decay_{i}_%"]
            code = str(code).strip()
            mode = code_to_mode.get(code)
            if mode is None and any(c.isalpha() for c in code):
                mode = "ClusterDecay"  # the source names the emitted cluster: 14C, 24NE, {+22}Ne, 34SI ...
            src_modes.append((mode or f"?{code}",
                              None if pd.isna(pct) else float(pct)))
        n_modes = []
        for b in g.objects(u, N.hasDecayBranch):
            br = g.value(b, N.branchingRatio)
            n_modes.append((local(g.value(b, N.hasDecayMode)), None if br is None else float(br)))
        if sorted(m for m, _ in src_modes) != sorted(m for m, _ in n_modes):
            flag(name, "decay modes", sorted(n_modes), sorted(src_modes))
        else:
            for m, p in src_modes:
                ns = [q for mm, q in n_modes if mm == m]
                if not any(close(q, p, 1e-6) for q in ns):
                    flag(name, "branching ratio", f"{m} {ns}", f"{m} {p}")

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "source_comparison.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["nuclide", "field", "nucleon", "source", "note"])
        w.writeheader()
        w.writerows(issues)
    print(f"source rows: {n_rows}; discrepancies: {len(issues)} in {len({i['nuclide'] for i in issues})} nuclides")
    for k, v in counts.most_common():
        print(f"  {v:5}  {k}")
    print("likely unit errors:", [(i["nuclide"], i["nucleon"], i["source"]) for i in issues if "unit error" in i["note"]])
    return issues


if __name__ == "__main__":
    main(*sys.argv[1:])
