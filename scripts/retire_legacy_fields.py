#!/usr/bin/env python3
"""Replace the last values that came from the unattributed 2017 spreadsheet.

After scripts/reload_from_source.py (ground-state half-lives, widths,
abundances) and scripts/rebuild_isomers.py (all isomers from NUBASE2020),
four kinds of ground-state data still came from the 2017 spreadsheet. This
script replaces or removes them:

  1. :spinParity of ground states: from the IAEA export (jp); where the
     export has none, NUBASE2020; otherwise removed.
  2. Half-life of the ground states that carried a spreadsheet half-life or
     width but that the IAEA export does not cover, or covers without a
     half-life: from NUBASE2020, with a note where NUBASE gives an
     estimate from systematics (#). A spreadsheet decay width is removed when
     NUBASE gives a time.
  3. :decayEnergy (Q-value, MeV) on the four branches that carried one: from
     the IAEA export's Q-values (qbm, qa, qec; keV).
  4. :decayModeOrAbundance (the spreadsheet's free text): removed, with its
     declaration. Decay modes, branching ratios and abundances are all
     asserted as structured data from cited sources.

Run scripts/add_display_annotations.py afterwards.

Usage: python3 scripts/retire_legacy_fields.py NUCLEON.ttl [--write]
"""
import os
import re
import sys
from decimal import Decimal

import pandas as pd
import rdflib
from rdflib import RDF, RDFS, Namespace

sys.path.insert(0, os.path.dirname(__file__))
import nubase  # noqa: E402
from reload_from_source import LIMIT_PREFIX, dbl, num, set_lines  # noqa: E402

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)
SRC = "source-data/IAEA nuclide livechart.csv"
OP = {">": "a lower limit (>)", ">=": "a lower limit (>=)", "<": "an upper limit (<)",
      "<=": "an upper limit (<=)", "~": "an approximate value (~)"}
Q_COLUMN = {"BetaMinusDecay": "qbm", "AlphaDecay": "qa", "ElectronCapture": "qec",
            "BetaPlusDecay": "qec", "PositronEmission": "qec"}
SYSTEMATICS = "NUBASE2020 gives this half-life as an estimate from systematics (#)."


def local(u):
    return str(u).split("#", 1)[1]


def main(path, write=False):
    g = rdflib.Graph()
    g.parse(path)
    text = open(path, encoding="utf-8").read()
    zsym = {int(g.value(e, N.atomicNumber)): str(g.value(e, RDFS.label))
            for e in g.subjects(RDF.type, N.Element) if (e, N.atomicNumber, None) in g}
    nb = nubase.parse(nubase_path(), zsym)
    d = pd.read_csv(SRC, dtype=str)
    d["name"] = d.symbol + "-" + (d.z.astype(int) + d.n.astype(int)).astype(str)
    src = d.set_index("name")
    gs = sorted(local(s) for s in g.subjects(N.massNumber, None) if (s, N.metastableState, None) not in g)
    counts = {"spin-parity from IAEA": 0, "spin-parity from NUBASE2020": 0, "spin-parity removed": 0,
              "half-life from NUBASE2020": 0, "decayEnergy from IAEA Q-value": 0}

    def edit(name, updates):
        nonlocal text
        s = text.index(f"###  {NS}{name}\n")
        e = text.index("\n\n\n", s) + 1
        text = text[:s] + set_lines(text[s:e], name, updates) + text[e:]

    for name in gs:
        u = N[name]
        updates = {}
        # 1. spin-parity
        jp = None
        if name in src.index and not pd.isna(src.at[name, "jp"]):
            jp, k = " ".join(str(src.at[name, "jp"]).split()), "spin-parity from IAEA"
        elif name in nb and nb[name]["spin_parity"]:
            jp, k = nb[name]["spin_parity"], "spin-parity from NUBASE2020"
        else:
            k = "spin-parity removed"
        cur = g.value(u, N.spinParity)
        if (str(cur) if cur is not None else None) != jp:
            updates[":spinParity"] = f'"{jp}"' if jp else None
            if jp or cur is not None:
                counts[k] += 1
        # 2. half-life for ground states the export does not cover
        uncovered = name not in src.index or pd.isna(src.at[name, "half_life"])
        had_legacy = (u, N.halfLife, None) in g or (u, N.decayWidthValue, None) in g
        if uncovered and had_legacy and name in nb and nb[name]["hl_s"] is not None \
                and (u, N.isStable, rdflib.Literal(True)) not in g:
            r = nb[name]
            updates.update({":halfLife": dbl(r["hl_s"]), ":halfLifeValue": num(r["hl_value"]),
                            ":halfLifeUnit": f'"{nubase.UNIT_NUCLEON.get(r["hl_unit"], r["hl_unit"])}"',
                            ":decayWidthValue": None, ":decayWidthUnit": None})
            note = []
            if r["hl_op"]:
                note.append(f"{LIMIT_PREFIX}{OP[r['hl_op']]}.")
            if r["hl_systematic"]:
                note.append(SYSTEMATICS)
            # set_lines handles one limit comment; systematics note is appended separately
            updates["rdfs:comment"] = f'"{note[0]}"@en' if note and note[0].startswith(LIMIT_PREFIX) else None
            counts["half-life from NUBASE2020"] += 1
            if updates:
                edit(name, updates)
            if r["hl_systematic"] and SYSTEMATICS not in str(list(g.objects(u, RDFS.comment))):
                s = text.index(f"###  {NS}{name}\n")
                e = text.index("\n\n\n", s) + 1
                blk = text[s:e]
                pad = " " * (len(name) + 2)
                text = text[:s] + blk[:-3] + f' ;\n{pad}rdfs:comment "{SYSTEMATICS}"@en .\n' + text[e:]
            continue
        if updates:
            edit(name, updates)

    # 3. decayEnergy from IAEA Q-values
    for b, q in list(g.subject_objects(N.decayEnergy)):
        parent = local(next(g.subjects(N.hasDecayBranch, b)))
        col = Q_COLUMN.get(local(g.value(b, N.hasDecayMode)))
        if col and parent in src.index and not pd.isna(src.at[parent, col]):
            mev = format((Decimal(src.at[parent, col]) / 1000).normalize(), "f")
            bname = local(b)
            s = text.index(f"###  {NS}{bname}\n")
            e = text.index("\n\n\n", s) + 1
            text = text[:s] + set_lines(text[s:e], bname, {":decayEnergy": mev}) + text[e:]
            counts["decayEnergy from IAEA Q-value"] += 1
            print(f"  {bname}: {q} -> {mev} MeV ({col})")

    # 4. remove :decayModeOrAbundance
    lit = r'"(?:[^"\\]|\\.)*"'
    n_before = text.count(":decayModeOrAbundance ")
    text = re.sub(rf"\n[ ]+:decayModeOrAbundance {lit} ;(?=\n)", "", text)
    text = re.sub(rf" ;\n[ ]+:decayModeOrAbundance {lit} \.\n", " .\n", text)
    s = text.index(f"###  {NS}decayModeOrAbundance\n")
    e = text.index("\n\n\n", s) + 3
    text = text[:s] + text[e:]
    left = len(re.findall(r"^[ ]+:decayModeOrAbundance ", text, flags=re.M))
    assert left == 0, f"{left} :decayModeOrAbundance assertions left"
    counts[":decayModeOrAbundance removed"] = n_before - 1

    for k, v in counts.items():
        print(f"{v:6}  {k}")
    if write:
        open(path, "w", encoding="utf-8").write(text)
        print(f"wrote {path}")


def nubase_path():
    return "source-data/nubase_3.mas20.txt"


if __name__ == "__main__":
    main(sys.argv[1], "--write" in sys.argv)
