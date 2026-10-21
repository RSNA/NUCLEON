#!/usr/bin/env python3
"""Reload half-lives, decay widths and natural abundances from the IAEA export.

NUCLEON's decay branches were built from source-data/IAEA nuclide
livechart.csv, but its half-lives, decay widths and abundances were carried
over from an earlier spreadsheet build (see validation/README.md). This script
replaces them with the export's values for every ground state the export
covers:

  :halfLife          export half_life_sec (xsd:double); for an entry given
                     as a decay width, hbar*ln2/Gamma as documented on
                     :halfLife
  :halfLifeValue,    the export's value and unit as given (unit codes mapped
  :halfLifeUnit      to NUCLEON's: Y -> a, us -> µs)
  :decayWidthValue,  for entries given as a width
  :decayWidthUnit
  rdfs:comment       the half-life limit note (>, <, >=, ~), rewritten from
                     the export's operator_hl
  :naturalAbundance  export abundance (percent)

It does not change a nuclide that NUCLEON flags stable while the export gives
a half-life, or the reverse; those are reported. Values the export lacks are
left as they are. Np-242m, an isomer and so not in the export, takes the
5.5-min half-life that the earlier build had attached to the Np-242 ground
state (agreeing with ICRP-107).

Edits the Turtle text in place. Run scripts/add_display_annotations.py
afterwards to refresh the class-level mirrors.

Usage: python3 scripts/reload_from_source.py NUCLEON.ttl [csv] [--write]
"""
import csv
import math
import re
import sys
from decimal import Decimal

import pandas as pd
import rdflib
from rdflib import Namespace

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)
HBAR_EV_S = 6.582119569e-16
WIDTH_EV = {"eV": 1.0, "keV": 1e3, "MeV": 1e6}
UNIT = {"Y": "a", "us": "µs"}
OPERATOR = {"GT": "a lower limit (>)", "GE": "a lower limit (>=)", "LT": "an upper limit (<)",
            "AP": "an approximate value (~)"}
LIMIT_PREFIX = "The source reports this half-life as "


def num(s):
    """Turtle numeric literal for a number given as text, without changing its value."""
    s = str(s).strip()
    if re.fullmatch(r"[-+]?\d+", s):
        return s
    if re.fullmatch(r"[-+]?\d*\.\d+", s):
        return s
    if re.fullmatch(r"[-+]?(\d+\.?\d*|\.\d+)[eE][-+]?\d+", s):
        return s  # Turtle reads exponent notation as xsd:double; keep the source's form
    return repr(float(s))


def dec(x):
    """xsd:decimal literal (no exponent)."""
    t = format(Decimal(str(x)).normalize(), "f")
    return t if "." in t else t + ".0"


def dbl(x):
    return f'"{float(x)!r}"^^xsd:double'


def set_lines(block, name, updates):
    """updates: {predicate_text: value_text or None}. Remove existing lines for
    each predicate and append the new ones before the block's final ' .'."""
    pad = " " * (len(name) + 2)
    for pred in updates:
        if pred == "rdfs:comment":
            pat = rf'\n[ ]+rdfs:comment "{re.escape(LIMIT_PREFIX)}[^"]*"@en'
        else:
            pat = rf"\n[ ]+{re.escape(pred)} [^\n]*?"
        block = re.sub(pat + r" ;(?=\n)", "", block)
        block = re.sub(r"( ;)?" + pat + r" \.\n$", " .\n", block)
    new = [f"{pad}{p} {v}" for p, v in updates.items() if v is not None]
    if new:
        assert block.endswith(" .\n"), name
        block = block[:-3] + " ;\n" + " ;\n".join(new) + " .\n"
    return block


def main(path, src="source-data/IAEA nuclide livechart.csv", write=False):
    g = rdflib.Graph()
    g.parse(path)
    d = pd.read_csv(src, dtype=str)
    text = open(path, encoding="utf-8").read()
    log, skipped = [], []

    def edit(name, updates, why):
        nonlocal text
        start = text.index(f"###  {NS}{name}\n")
        end = text.index("\n\n\n", start) + 1
        text = text[:start] + set_lines(text[start:end], name, updates) + text[end:]
        log.append({"nuclide": name, "change": why,
                    "fields": "; ".join(f"{p}={v}" for p, v in updates.items())})

    for _, r in d.iterrows():
        name = f"{r.symbol}-{int(r.z) + int(r.n)}"
        u = N[name]
        if (u, N.massNumber, None) not in g:
            continue
        stable_src = str(r.half_life).strip().upper() == "STABLE"
        stable_n = (u, N.isStable, rdflib.Literal(True)) in g
        updates, why = {}, []

        if not pd.isna(r.abundance):
            old = g.value(u, N.naturalAbundance)
            if old is None or Decimal(str(old)) != Decimal(r.abundance):
                updates[":naturalAbundance"] = dec(r.abundance)
                why.append("abundance")

        if not stable_src and not pd.isna(r.half_life) and not pd.isna(r.unit_hl):
            if stable_n:
                skipped.append((name, "NUCLEON flags stable; export gives "
                                f"{r.half_life} {r.unit_hl}"))
            else:
                unit = r.unit_hl.strip()
                if unit in WIDTH_EV:
                    gamma = float(r.half_life) * WIDTH_EV[unit]
                    hl = HBAR_EV_S * math.log(2) / gamma
                    upd = {":halfLife": dbl(hl), ":decayWidthValue": dec(r.half_life),
                           ":decayWidthUnit": f'"{unit}"', ":halfLifeValue": None, ":halfLifeUnit": None}
                else:
                    hl = float(r.half_life_sec)
                    upd = {":halfLife": dbl(hl), ":halfLifeValue": num(r.half_life),
                           ":halfLifeUnit": f'"{UNIT.get(unit, unit)}"',
                           ":decayWidthValue": None, ":decayWidthUnit": None}
                op = None if pd.isna(r.operator_hl) else OPERATOR.get(r.operator_hl.strip())
                upd["rdfs:comment"] = f'"{LIMIT_PREFIX}{op}."@en' if op else None
                # only record a change if something differs
                cur = {":halfLife": g.value(u, N.halfLife), ":halfLifeValue": g.value(u, N.halfLifeValue),
                       ":halfLifeUnit": g.value(u, N.halfLifeUnit), ":decayWidthValue": g.value(u, N.decayWidthValue)}
                changed = (cur[":halfLife"] is None or not math.isclose(float(cur[":halfLife"]), hl, rel_tol=1e-9)
                           or (":halfLifeValue" in upd and upd[":halfLifeValue"] is not None and
                               (cur[":halfLifeValue"] is None or float(cur[":halfLifeValue"]) != float(r.half_life)
                                or str(cur[":halfLifeUnit"]) != UNIT.get(unit, unit)))
                           or (unit in WIDTH_EV and (cur[":decayWidthValue"] is None or
                                                     float(cur[":decayWidthValue"]) != float(r.half_life))))
                has_limit = any(str(c).startswith(LIMIT_PREFIX) for c in g.objects(u, rdflib.RDFS.comment))
                if changed or bool(op) != has_limit:
                    updates.update(upd)
                    why.append("half-life" if unit not in WIDTH_EV else "decay width")
        elif stable_src and not stable_n:
            skipped.append((name, "export says STABLE; NUCLEON flags radioactive"))
        elif not stable_src and pd.isna(r.half_life) and (u, N.halfLife, None) in g:
            skipped.append((name, "export gives no half-life; NUCLEON value kept"))

        if updates:
            edit(name, updates, ", ".join(why))

    # Np-242m: the earlier build had the isomer's 5.5 min on the ground state
    edit("Np-242m", {":halfLife": dbl(330.0), ":halfLifeValue": "5.5", ":halfLifeUnit": '"m"'},
         "half-life swapped back from ground state (5.5 min; ICRP-107 agrees)")

    c = {}
    for e in log:
        for w in e["change"].split(", "):
            c[w] = c.get(w, 0) + 1
    print(f"nuclides changed: {len(log)}; by field: {c}")
    print(f"not changed ({len(skipped)}):")
    for s in skipped:
        print("  ", *s)
    if write:
        open(path, "w", encoding="utf-8").write(text)
        with open("validation/reload_from_source_log.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["nuclide", "change", "fields"])
            w.writeheader()
            w.writerows(log)
        print(f"wrote {path} and validation/reload_from_source_log.csv")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--write"]
    main(*args, write="--write" in sys.argv)
