#!/usr/bin/env python3
"""Rebuild NUCLEON's metastable isomers from NUBASE2020.

Replaces every isomer, and its decay branches, with the isomers of NUBASE2020
(source-data/nubase_3.mas20.txt): every state with index i = 1-6 that has a
half-life (or is observationally stable). Isobaric analog states (i = 8, 9)
and levels without a half-life are not isomers and are excluded.

Naming follows NUBASE's own state index, so each IRI traces to one NUBASE row:
i = 1 (m) -> X-Am, i = 2 (n) -> X-Am2, i = 3 (p) -> X-Am3, ... i = 6 (x) -> X-Am6.
:metastableState is i. NUBASE numbers states by excitation energy, so a
long-lived isomer above shorter-lived ones has a higher suffix (e.g. the
160-d Lu-177 isomer is Lu-177m3).

For each isomer: atomic, mass and neutron number; :excitationEnergy (keV);
half-life in seconds plus value and unit as given; limit and systematics
notes; spin-parity; labels; and one :DecayBranch per decay mode in the
NUBASE decay string, with daughters derived as for ground states (isomeric
transition -> ground state of the same nuclide; cluster decay -> parent
minus the named cluster; fission -> no daughter). Observationally stable
isomers (Ta-180m) get :isStable true and their natural abundance instead.

DICOM codes, DICOM context groups and skos:exactMatch on an existing isomer
are carried to the NUBASE state with the closest half-life in the same
nucleus. Every existing isomer's fate is written to
validation/isomer_mapping.csv.

Run scripts/add_display_annotations.py afterwards.

Usage: python3 scripts/rebuild_isomers.py NUCLEON.ttl [--write]
"""
import collections
import csv
import itertools
import math
import os
import re
import sys
from decimal import Decimal

import rdflib
from rdflib import RDF, RDFS, Namespace

sys.path.insert(0, os.path.dirname(__file__))
import nubase  # noqa: E402

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")
NUBASE_FILE = "source-data/nubase_3.mas20.txt"
SAME_STATE = math.log10(3)   # half-lives within 3x: treat as the same state for code carry-over
OP_HL = {">": "a lower limit (>)", ">=": "a lower limit (>=)", "<": "an upper limit (<)",
         "<=": "an upper limit (<=)", "~": "an approximate value (~)"}
OP_BR = {"~": "approximate (~{v}%)", "<": "an upper limit (<{v}%); the limit is stored",
         "<=": "an upper limit (<={v}%); the limit is stored",
         ">": "a lower limit (>{v}%); the limit is stored", ">=": "a lower limit (>={v}%); the limit is stored"}

ELEMENT_NAMES = {}   # filled from the ontology: symbol -> English name


def local(u):
    return str(u).split("#", 1)[1]


def dec(x):
    t = format(Decimal(str(x)).normalize(), "f")
    return t


def block_span(text, name):
    start = text.index(f"###  {NS}{name}\n")
    end = text.index("\n\n\n", start) + 3
    return start, end


def render(name, pairs):
    """Protégé-style block: pairs = [(predicate, [objects])]."""
    head = f"###  {NS}{name}\n:{name} "
    pad = " " * (len(name) + 2)
    lines = []
    for k, (p, objs) in enumerate(pairs):
        lead = head if k == 0 else pad
        objpad = " " * (len(lead.split("\n")[-1]) + len(p) + 1)
        lines.append(lead + p + " " + (" ,\n" + objpad).join(objs))
    return " ;\n".join(lines) + " .\n\n\n"


def main(path, write=False):
    g = rdflib.Graph()
    g.parse(path)
    text = open(path, encoding="utf-8").read()

    zsym, gs_names = {}, set()
    for e in g.subjects(RDF.type, N.Element):
        if (e, N.atomicNumber, None) in g:
            sym = str(g.value(e, RDFS.label))
            zsym[int(g.value(e, N.atomicNumber))] = sym
            en = [str(l) for l in g.objects(e, SKOS.altLabel) if getattr(l, "language", None) == "en"]
            ELEMENT_NAMES[sym] = min(en, key=len) if en else sym
    for n in g.subjects(N.massNumber, None):
        if (n, N.metastableState, None) not in g:
            gs_names.add(local(n))
    za_of = {}
    for n in gs_names:
        za_of[(int(g.value(N[n], N.atomicNumber)), int(g.value(N[n], N.massNumber)))] = n

    nb = nubase.parse(NUBASE_FILE, zsym)
    iso = {k: v for k, v in nb.items() if v["i"] > 0}
    mode_label = {}
    for b in g.subjects(RDF.type, N.DecayBranch):
        mode_label.setdefault(local(g.value(b, N.decayMode)), str(g.value(b, RDFS.label)).split(" ", 1)[1])
    mode_label.setdefault("ClusterDecay", "cluster branch")
    dz = {local(m): g.value(m, N.changeInAtomicNumber) for m in g.subjects(RDF.type, N.DecayMode)}
    da = {local(m): g.value(m, N.changeInMassNumber) for m in g.subjects(RDF.type, N.DecayMode)}

    # ---- map existing isomers to NUBASE states (for carrying codes) ----------------------
    old = sorted(local(s) for s in g.subjects(N.metastableState, None))
    by_nuc = collections.defaultdict(list)
    for k in old:
        by_nuc[re.sub(r"m\d?$", "", k)].append(k)
    nb_nuc = collections.defaultdict(list)
    for k, v in iso.items():
        nb_nuc[f"{v['symbol']}-{v['A']}"].append(k)

    def lr(o, c):
        ho, hc = g.value(N[o], N.halfLife), iso[c]["hl_s"]
        return None if ho is None or hc is None else nubase.log_ratio(float(ho), hc)

    mapping = {}
    for nuc, olds in by_nuc.items():
        cands = nb_nuc.get(nuc, [])
        k = min(len(cands), len(olds))
        best = None
        for sub in itertools.permutations(olds, k):
            for perm in itertools.permutations(cands, k):
                pairs = [(o, c, lr(o, c)) for o, c in zip(sub, perm)]
                cost = sum(min(x if x is not None else 1.0, 3) for *_, x in pairs)
                if best is None or cost < best[0]:
                    best = (cost, pairs)
        for o, c, x in (best[1] if best else []):
            mapping[o] = (c, x)
    carried = collections.defaultdict(list)   # new name -> [(pred, obj_text)]
    rows = []
    for o in old:
        c, x = mapping.get(o, (None, None))
        same = c is not None and (x is None or x <= SAME_STATE)
        codes = []
        for p, pt in ((SKOS.exactMatch, "skos:exactMatch"), (N.dicomCode, ":dicomCode"),
                      (N.dicomContextGroup, ":dicomContextGroup")):
            for obj in g.objects(N[o], p):
                codes.append((pt, obj.n3() if isinstance(obj, rdflib.URIRef) else f'"{obj}"'))
        if codes and not same:
            raise SystemExit(f"{o} carries codes but has no reliable NUBASE match; resolve by hand")
        if same and codes:
            carried[c] += codes
        status = ("removed: no NUBASE2020 state" if c is None else
                  "same IRI" if c == o else f"renamed to {c}")
        if c is not None and not same:
            status += "; 2017 half-life differed from NUBASE2020 by more than 3x (NUBASE2020 value used)"
        rows.append({"old_iri": o, "old_half_life_s": g.value(N[o], N.halfLife),
                     "new_iri": c or "",
                     "nubase_label": iso[c]["nubase_label"].strip() if c else "",
                     "nubase_half_life": f"{iso[c]['hl_value'] or ''} {iso[c]['hl_unit'] or ''}".strip() if c else "",
                     "status": status, "codes_carried": "; ".join(f"{p} {v}" for p, v in codes)})
    for r in rows:   # flag old IRIs whose name is now used by a different state
        if r["old_iri"] in iso and r["new_iri"] != r["old_iri"]:
            r["status"] += f"; IRI now denotes NUBASE state {iso[r['old_iri']]['nubase_label'].strip()}"

    # ---- remove old isomer blocks and their branches ------------------------------------
    for o in old:
        for b in list(g.objects(N[o], N.hasDecayBranch)):
            s, e = block_span(text, local(b))
            text = text[:s] + text[e:]
        s, e = block_span(text, o)
        text = text[:s] + text[e:]

    # ---- render new isomers ----------------------------------------------------------------
    unknown_modes = collections.Counter()
    new_by_gs = collections.defaultdict(list)
    n_branches = 0
    for name, r in sorted(iso.items(), key=lambda kv: (kv[1]["Z"], kv[1]["A"], kv[1]["i"])):
        sym, a, z, i = r["symbol"], r["A"], r["Z"], r["i"]
        gs = za_of[(z, a)]
        suf = nubase.suffix(i)
        branches, comments = [], []
        if r["hl_op"]:
            comments.append(f"The source reports this half-life as {OP_HL[r['hl_op']]}.")
        if r["hl_systematic"]:
            comments.append("NUBASE2020 gives this half-life as an estimate from systematics (#).")
        if r["exc_systematic"]:
            comments.append("NUBASE2020 gives this excitation energy as an estimate from systematics (#).")
        abundance = None
        if r["stable"]:
            m = re.search(r"IS=([0-9.]+)", r["decay"])
            abundance = m.group(1) if m else None
            comments.append("Observationally stable in NUBASE2020; no decay has been observed.")
        else:
            for k, (mode, op, val) in enumerate(nubase.decay_modes(r["decay"]), 1):
                sign_note = None
                if mode == "?B":
                    # NUBASE gives "B" (beta decay, sign unspecified): take the ground state's beta mode
                    gs_modes = {local(g.value(b, N.decayMode)) for b in g.objects(N[gs], N.hasDecayBranch)}
                    mode = "BetaMinusDecay" if "BetaMinusDecay" in gs_modes else "BetaPlusDecay"
                    sign_note = ("NUBASE2020 gives this mode as 'B' (beta decay, sign unspecified); "
                                 "assigned from the ground state's beta mode.")
                if mode.startswith("?"):
                    unknown_modes[mode] += 1
                    continue
                cluster = None
                if mode.startswith("CD:"):
                    cluster = mode[3:]
                    cm = nubase.CLUSTER.match(cluster)
                    ca, csym = int(cm.group(1)), cm.group(2)
                    cz = next(zz for zz, ss in zsym.items() if ss == csym)
                    mode, d_z, d_a = "ClusterDecay", -cz, -ca
                else:
                    d_z, d_a = dz.get(mode), da.get(mode)
                bid = f"{name}_b{len(branches) + 1}"
                props = [("rdf:type", ["owl:NamedIndividual", ":DecayBranch"]), (":decayMode", [f":{mode}"])]
                if d_z is not None and d_a is not None:
                    dn = za_of.get((z + int(d_z), a + int(d_a)))
                    if dn:
                        props.append((":decayProduct", [f":{dn}"]))
                if val is not None:
                    props.append((":branchingRatio", [dec(val)]))
                bc = [sign_note] if sign_note else []
                if cluster:
                    bc.append(f"Emitted cluster: {cluster}.")
                if op in OP_BR and val is not None:
                    bc.append("Branching ratio reported by the source as " + OP_BR[op].format(v=dec(val)) + ".")
                for c in bc:
                    props.append(("rdfs:comment", [f'"{c}"@en']))
                props.append(("rdfs:label", [f'"{name} {mode_label[mode]}"@en']))
                branches.append((bid, props))
        n_branches += len(branches)

        p = [("rdf:type", ["owl:Class", "owl:NamedIndividual", ":Nuclide"]),
             ("rdfs:subClassOf", [":Nuclide"])]
        if branches:
            p.append((":hasDecayBranch", [f":{b}" for b, _ in branches]))
        p.append((":atomicNumber", [str(z)]))
        if r["exc_keV"] is not None:
            p.append((":excitationEnergy", [dec(r["exc_keV"])]))
        if r["hl_s"] is not None and not r["stable"]:
            p.append((":halfLife", [f'"{r["hl_s"]!r}"^^xsd:double']))
        p.append((":isStable", ["true" if r["stable"] else "false"]))
        p.append((":massNumber", [str(a)]))
        p.append((":metastableState", [str(i)]))
        if abundance:
            p.append((":naturalAbundance", [dec(abundance)]))
        p.append((":neutronNumber", [str(a - z)]))
        for c in comments:
            p.append(("rdfs:comment", [f'"{c}"@en']))
        p.append(("rdfs:label", [f'"{name}"']))
        el = ELEMENT_NAMES[sym]
        p.append(("skos:altLabel", [f'"{el}-{a}{suf}"@en', f'"({a}{suf}){sym}"', f'"({a}{suf}){el}"@en']))
        by_pred = collections.defaultdict(list)
        for pt, v in carried.get(name, []):
            by_pred[pt].append(v)
        for pt in ("skos:exactMatch", ":dicomCode", ":dicomContextGroup"):
            if by_pred[pt]:
                p.append((pt, sorted(set(by_pred[pt]))))
        if r["hl_value"] is not None and not r["stable"]:
            p.append((":halfLifeUnit", [f'"{nubase.UNIT_NUCLEON.get(r["hl_unit"], r["hl_unit"])}"']))
            p.append((":halfLifeValue", [r["hl_value"]]))
        if r["spin_parity"]:
            p.append((":spinParity", [f'"{r["spin_parity"]}"']))
        block = render(name, p) + "".join(render(b, bp) for b, bp in branches)
        new_by_gs[gs].append(block)

    # ---- insert each nucleus's isomers after its ground state and branches -------------
    for gs, blocks in new_by_gs.items():
        s, e = block_span(text, gs)
        while True:
            nxt = text.find("###  " + NS, e)
            if nxt == -1:
                break
            nm = text[nxt + 5 + len(NS):text.index("\n", nxt)]
            if not nm.startswith(gs + "_b"):
                break
            e = text.index("\n\n\n", nxt) + 3
        text = text[:e] + "".join(blocks) + text[e:]

    # ---- element :hasNuclide lists ----------------------------------------------------------
    members = collections.defaultdict(list)
    for n in gs_names:
        members[re.match(r"[A-Za-z]+", n).group(0)].append(n)
    for n in iso:
        members[iso[n]["symbol"]].append(n)
    for sym, ms in members.items():
        s, e = block_span(text, sym)
        blk = text[s:e]
        m = re.search(r"\n([ ]+):hasNuclide (.*?)( [;.])\n", blk, flags=re.S)
        pad = m.group(1) + " " * len(":hasNuclide ")
        new = "\n" + m.group(1) + ":hasNuclide " + (" ,\n" + pad).join(f":{x}" for x in sorted(ms)) + m.group(3) + "\n"
        blk = blk[:m.start()] + new + blk[m.end():]
        text = text[:s] + blk + text[e:]

    # ---- :excitationEnergy declaration --------------------------------------------------------
    decl_iri = f"###  {NS}excitationEnergy\n"
    if decl_iri not in text:
        anchor = f"###  {NS}groupNumber\n"
        decl = (decl_iri + ":excitationEnergy rdf:type owl:DatatypeProperty ;\n"
                "                  rdfs:domain :Nuclide ;\n"
                "                  rdfs:range xsd:decimal ;\n"
                '                  rdfs:comment "Excitation energy of a metastable isomer above the ground state, in keV, from NUBASE2020. Distinguishes isomers of the same nuclide; :metastableState gives the NUBASE state index (1 = m, 2 = n, ...)."@en ;\n'
                '                  rdfs:label "excitation energy (keV)"@en .\n\n\n')
        text = text.replace(anchor, decl + anchor, 1)

    st = collections.Counter(r["status"].split(";")[0] for r in rows)
    print(f"NUBASE2020 isomers written: {len(iso)}; their decay branches: {n_branches}")
    print(f"existing isomers: {len(old)} -> {dict(st)}")
    print(f"old IRIs now naming a different state: {sum(1 for r in rows if 'IRI now denotes' in r['status'])}")
    print(f"isomers carrying codes: {len(carried)} {sorted(carried)}")
    print(f"unknown decay codes skipped: {dict(unknown_modes)}")
    if write:
        open(path, "w", encoding="utf-8").write(text)
        with open("validation/isomer_mapping.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path} and validation/isomer_mapping.csv")


if __name__ == "__main__":
    main(sys.argv[1], "--write" in sys.argv)
