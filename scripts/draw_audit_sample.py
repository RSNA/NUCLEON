#!/usr/bin/env python3
"""Draw a stratified random sample of nuclides for manual audit.

Every nuclide is assigned to exactly one stratum (first match in the order
below), then a fixed number is drawn from each with a fixed seed, so the sample
is reproducible. The reviewer compares the NUCLEON values in the sheet with the
IAEA LiveChart of Nuclides web interface (https://nds.iaea.org/livechart/).

Outputs (in validation/):
  audit_sample.csv      the review sheet; reviewer columns are blank
  audit_sample_key.csv  which sampled nuclides the automated comparisons
                        flagged; keep from the reviewer until the audit ends, so
                        the audit stays blind

Usage: python3 scripts/draw_audit_sample.py [NUCLEON.ttl] [outdir] [seed]
"""
import csv
import os
import random
import sys

import rdflib
from rdflib import RDFS, Namespace

N = Namespace("https://rsna.org/radlex/nucleon#")

STRATA = [  # (name, target n, description)
    ("dicom", 20, "DICOM-coded radionuclide (CID 18 or CID 4020)"),
    ("isomer", 15, "metastable isomer"),
    ("stable", 10, "flagged stable"),
    ("heavy", 10, "Z >= 89 (actinides and superheavy)"),
    ("limit", 5, "half-life given by the source as a limit or approximation"),
    ("short", 10, "half-life < 1 s or given as a decay width"),
    ("multi", 15, "two or more decay branches"),
    ("single", 15, "one decay branch, half-life >= 1 s"),
]


def local(u):
    return str(u).split("#", 1)[1]


def main(ttl="NUCLEON.ttl", outdir="validation", seed="20261021"):
    g = rdflib.Graph()
    g.parse(ttl)
    nuclides = sorted(local(s) for s in g.subjects(N.massNumber, None))

    def stratum(name):
        u = N[name]
        hl = g.value(u, N.halfLife)
        nb = len(list(g.objects(u, N.hasDecayBranch)))
        comment = " ".join(str(c) for c in g.objects(u, RDFS.comment)).lower()
        if (u, N.dicomCode, None) in g:
            return "dicom"
        if (u, N.metastableState, None) in g:
            return "isomer"
        if (u, N.isStable, rdflib.Literal(True)) in g:
            return "stable"
        if int(g.value(u, N.atomicNumber)) >= 89:
            return "heavy"
        if "limit" in comment or "approximat" in comment:
            return "limit"
        if (u, N.decayWidthValue, None) in g or (hl is not None and float(hl) < 1):
            return "short"
        if nb >= 2:
            return "multi"
        if nb == 1:
            return "single"
        return None

    pools = {s: [] for s, _, _ in STRATA}
    for n in nuclides:
        s = stratum(n)
        if s:
            pools[s].append(n)

    rng = random.Random(int(seed))
    sample = []
    for s, k, desc in STRATA:
        pick = sorted(rng.sample(pools[s], min(k, len(pools[s]))))
        sample += [(s, desc, len(pools[s]), n) for n in pick]

    # automated flags, for the separate key file
    flagged = {}
    for fname, col in (("source_comparison.csv", "nuclide"), ("halflife_vs_icrp107.csv", "nuclide"),
                       ("branching_vs_icrp107.csv", "nuclide")):
        p = os.path.join(outdir, fname)
        if os.path.exists(p):
            for r in csv.DictReader(open(p, encoding="utf-8")):
                if fname == "source_comparison.csv" or r.get("flag"):
                    flagged.setdefault(r[col], set()).add(f"{fname}:{r.get('field') or r.get('mode_family') or 'half-life'}")

    rows, key = [], []
    for i, (s, desc, pool, n) in enumerate(sample, 1):
        u = N[n]
        val = lambda p: "" if g.value(u, p) is None else str(g.value(u, p))
        rows.append({
            "id": i, "stratum": s, "nuclide": n,
            "nucleon_half_life": f"{val(N.halfLifeValue)} {val(N.halfLifeUnit)}".strip(),
            "nucleon_half_life_s": val(N.halfLife),
            "nucleon_decay_width": f"{val(N.decayWidthValue)} {val(N.decayWidthUnit)}".strip(),
            "nucleon_is_stable": val(N.isStable),
            "nucleon_natural_abundance_pct": val(N.naturalAbundance),
            "nucleon_decay_summary": val(N.decaySummary),
            "nucleon_dicom_code": val(N.dicomCode),
            "iaea_half_life": "", "half_life_match": "",
            "iaea_decay_modes_and_ratios": "", "modes_and_ratios_match": "",
            "daughters_correct": "", "abundance_match": "",
            "reviewer": "", "notes": "",
        })
        key.append({"id": i, "nuclide": n, "automated_flags": "; ".join(sorted(flagged.get(n, [])))})

    os.makedirs(outdir, exist_ok=True)
    for fname, data in (("audit_sample.csv", rows), ("audit_sample_key.csv", key)):
        with open(os.path.join(outdir, fname), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys()))
            w.writeheader()
            w.writerows(data)
    print(f"seed {seed}; sample {len(rows)}")
    for s, k, desc in STRATA:
        print(f"  {s:7} {sum(1 for r in rows if r['stratum'] == s):3} of {len(pools[s]):5}  {desc}")
    print(f"sampled nuclides with an automated flag: {sum(1 for k in key if k['automated_flags'])}")


if __name__ == "__main__":
    main(*sys.argv[1:])
