#!/usr/bin/env python3
"""Cross-check NUCLEON's gamma emissions (ICRP-107) against ENSDF.

For every NUCLEON gamma-ray emission with yield >= 5% per decay, look for a
gamma of the same energy (within 0.5 keV or 0.1%) in the ENSDF decay datasets
of the same parent (paceENSDF 0.6.3: alpha, beta-minus and EC/beta-plus
datasets; isomeric-transition gammas, e.g. Tc-99m's 140.5 keV, cannot be
checked). ICRP-107 was compiled from ENSDF evaluations of about 2005, so this
mainly tests the transfer of the data into NUCLEON.

Usage: python3 scripts/compare_emissions.py [NUCLEON.ttl] [outdir]
"""
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict, Counter

import rdflib
from rdflib import RDF, Namespace

N = Namespace("https://rsna.org/radlex/nucleon#")
MIN_YIELD = 0.05


def local(u):
    return str(u).split("#", 1)[1]


def main(ttl="NUCLEON.ttl", outdir="validation"):
    g = rdflib.Graph()
    g.parse(ttl)
    import paceENSDF
    root = os.path.join(os.path.dirname(paceENSDF.__file__), "ENSDF_JSON")
    gammas = defaultdict(list)   # (symbol, A, excited?) -> [(E, dataset)]
    for f in glob.glob(os.path.join(root, "*.json")):
        d = json.load(open(f, encoding="utf-8"))
        sym = re.match(r"[A-Za-z]+", d["parentID"]).group(0)
        key = (sym, int(d["parentAtomicMass"]), d["levelEnergyParentDecay"] != "0")
        for lev in d["levelScheme"]:
            for gm in lev.get("gammaDecay") or []:
                try:
                    gammas[key].append((float(gm["gammaEnergy"]), os.path.basename(f)))
                except (TypeError, ValueError):
                    pass   # unknown energy ("X") or none
    rows = []
    for n in sorted({local(s) for s in g.subjects(N.hasEmission, None)}):
        sym, rest = n.split("-")
        a = int(re.match(r"\d+", rest).group(0))
        key = (sym, a, (N[n], N.metastableState, None) in g)
        for e in g.objects(N[n], N.hasEmission):
            if g.value(e, N.hasRadiationType) != N.GammaRay or float(g.value(e, N.emissionYield)) < MIN_YIELD:
                continue
            en = float(g.value(e, N.emissionEnergy))
            cand = gammas.get(key)
            if not cand:
                verdict, best = "not checkable: no ENSDF decay dataset with gammas", None
            else:
                best = min(cand, key=lambda c: abs(c[0] - en))
                verdict = "agree" if abs(best[0] - en) <= max(0.5, 0.001 * en) else "no ENSDF gamma at this energy"
            rows.append({"nuclide": n, "emission": local(e), "energy_keV": en,
                         "yield": float(g.value(e, N.emissionYield)),
                         "nearest_ensdf_keV": best[0] if best else "", "verdict": verdict})
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "emission_comparison.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    c = Counter(r["verdict"] for r in rows)
    print(f"gamma emissions with yield >= {MIN_YIELD:.0%}: {len(rows)}")
    for k, v in c.most_common():
        print(f"  {v:5}  {k}")
    print("no match:", [(r["nuclide"], r["energy_keV"], r["nearest_ensdf_keV"]) for r in rows
                        if r["verdict"].startswith("no ENSDF")])
    print("not checkable:", sorted({r["nuclide"] for r in rows if r["verdict"].startswith("not")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
