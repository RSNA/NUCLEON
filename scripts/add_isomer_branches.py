#!/usr/bin/env python3
"""Add decay branches to metastable isomers in NUCLEON.ttl.

The IAEA ground-state export used to build NUCLEON's decay branches does not
cover isomers, so 713 of 716 isomers were released without a :DecayBranch.
Most carry the source's decay text in :decayModeOrAbundance, e.g.
"%IT=99.9963 , % beta =0.0037". This script parses that text, creates one
:DecayBranch per reported mode, derives the daughter from the mode's
:changeInAtomicNumber / :changeInMassNumber exactly as for ground states,
and writes the result back into the Turtle file as text, so the existing
Protégé-style layout is preserved and the diff shows only what changed.

Conventions (same as the ground-state build):
  * The daughter is the ground state of nuclide (Z + dZ, A + dA); an isomeric
    transition therefore resolves to the ground state of the same nuclide.
  * Fission modes get no :decayProduct.
  * A daughter absent from the ontology is left unasserted.
  * "?" ratios are omitted. Approximate ("~") and limit ("<", ">") ratios are
    stored as their numeric value with an rdfs:comment saying so.

Isomers that already have branches are re-derived as well; the script
reports any difference and keeps existing extra assertions such as
:decayEnergy.

Usage: python3 scripts/add_isomer_branches.py NUCLEON.ttl [--write]
"""
import re
import sys
from decimal import Decimal

import rdflib
from rdflib import RDF, RDFS, Namespace

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)

# Source mode text (after the leading %) -> NUCLEON DecayMode local name
MODE_MAP = {
    "IT": "IsomericTransition",
    "K+% positron": "BetaPlusDecay",
    "K+ positron": "BetaPlusDecay",
    "beta": "BetaMinusDecay",
    "alpha": "AlphaDecay",
    "K": "ElectronCapture",
    "Kp": "BetaPlusDelayedProtonEmission",
    "positron P": "BetaPlusDelayedProtonEmission",
    "K alpha": "BetaPlusDelayedAlphaEmission",
    "SF": "SpontaneousFission",
    "beta n": "BetaMinusDelayedNeutronEmission",
    "beta alpha": "BetaMinusDelayedAlphaEmission",
    "beta SF": "BetaMinusDelayedFission",
    "p": "ProtonEmission",
    "n": "NeutronEmission",
}

TOKEN = re.compile(
    r"%\s*(?P<mode>[A-Za-z][A-Za-z+% ]*?)\s*(?P<op>=|<|>|~)\s*"
    r"(?P<val>\?|[0-9]*\.?[0-9]+(?:[Ee][-+]?[0-9]+)?)"
)

OP_NOTE = {
    "~": "Branching ratio reported by the source as approximate (~{v}%).",
    "<": "Branching ratio reported by the source as an upper limit (<{v}%); the limit is stored.",
    ">": "Branching ratio reported by the source as a lower limit (>{v}%); the limit is stored.",
}


def local(u):
    return str(u).split("#", 1)[1]


def parse(text):
    """Return [(mode_local, op, value_str_or_None)] in source order."""
    out, unknown = [], []
    for m in TOKEN.finditer(text):
        key = " ".join(m.group("mode").split())
        mode = MODE_MAP.get(key)
        if mode is None:
            unknown.append(key)
            continue
        val = None if m.group("val") == "?" else format(Decimal(m.group("val")).normalize(), "f")
        out.append((mode, m.group("op"), val))
    return out, unknown


def main(path, write):
    g = rdflib.Graph()
    g.parse(path)

    sym = {int(g.value(e, N.atomicNumber)): str(g.value(e, RDFS.label))
           for e in g.subjects(RDF.type, N.Element) if (e, N.atomicNumber, None) in g}
    nuclides = {local(n) for n in g.subjects(RDF.type, N.Nuclide)}
    dz = {local(m): g.value(m, N.changeInAtomicNumber) for m in g.subjects(RDF.type, N.DecayMode)}
    da = {local(m): g.value(m, N.changeInMassNumber) for m in g.subjects(RDF.type, N.DecayMode)}
    mode_label = {}
    for b in g.subjects(RDF.type, N.DecayBranch):
        mode_label.setdefault(local(g.value(b, N.decayMode)),
                              str(g.value(b, RDFS.label)).split(" ", 1)[1])
    mode_label.setdefault("BetaPlusDelayedProtonEmission", "beta+-delayed-proton branch")

    text = open(path, encoding="utf-8").read()
    report = {"added_branches": 0, "isomers_patched": 0, "no_text": [], "unparsed": [],
              "unknown_modes": {}, "no_daughter": [], "existing_diff": []}

    isomers = sorted(local(n) for n in g.subjects(N.metastableState, None))
    for name in isomers:
        n = N[name]
        src = g.value(n, N.decayModeOrAbundance)
        if src is None:
            report["no_text"].append(name)
            continue
        branches, unknown = parse(str(src))
        for k in unknown:
            report["unknown_modes"].setdefault(k, []).append(name)
        if not branches:
            report["unparsed"].append((name, str(src)))
            continue

        z, a = int(g.value(n, N.atomicNumber)), int(g.value(n, N.massNumber))
        existing = sorted(g.objects(n, N.hasDecayBranch), key=str)

        new_blocks, ids = [], []
        for i, (mode, op, val) in enumerate(branches, 1):
            bid = f"{name}_b{i}"
            ids.append(bid)
            daughter = None
            if dz[mode] is not None and da[mode] is not None:
                z2, a2 = z + int(dz[mode]), a + int(da[mode])
                cand = f"{sym.get(z2, '?')}-{a2}"
                if cand in nuclides:
                    daughter = cand
                else:
                    report["no_daughter"].append((bid, cand))
            props = [f":decayMode :{mode}"]
            if daughter:
                props.append(f":decayProduct :{daughter}")
            if val is not None:
                props.append(f":branchingRatio {val}")
            # keep extra assertions (e.g. :decayEnergy) from an existing branch of the same mode
            for eb in existing:
                if local(g.value(eb, N.decayMode)) == mode:
                    for p, o in g.predicate_objects(eb):
                        if p == N.decayEnergy:
                            props.append(f":decayEnergy {o}")
                    old_ratio = g.value(eb, N.branchingRatio)
                    if old_ratio is not None and val is not None and Decimal(str(old_ratio)) != Decimal(val):
                        report["existing_diff"].append((local(eb), str(old_ratio), val))
            if op in OP_NOTE and val is not None:
                props.append(f'rdfs:comment "{OP_NOTE[op].format(v=val)}"@en')
            props.append(f'rdfs:label "{name} {mode_label[mode]}"@en')
            head = f":{bid} rdf:type owl:NamedIndividual ,\n"
            ind = " " * (len(bid) + 2)
            block = (f"###  {NS}{bid}\n" + head + " " * (len(bid) + 2 + 9) + ":DecayBranch ;\n"
                     + " ;\n".join(ind + p for p in props) + " .\n")
            new_blocks.append(block)
        report["added_branches"] += len(ids) - len(existing)

        # --- rewrite the nuclide block -------------------------------------------------
        start = text.index(f"###  {NS}{name}\n")
        end = text.index("\n\n\n", start) + 3
        block = text[start:end]
        pad = " " * (len(name) + 2)
        block = re.sub(r"\n\s*:hasDecayBranch [^;]*;", "", block)
        hdb = pad + ":hasDecayBranch " + (" ,\n" + pad + " " * 16).join(f":{b}" for b in ids) + " ;"
        block = block.replace(f"rdfs:subClassOf :Nuclide ;\n", f"rdfs:subClassOf :Nuclide ;\n{hdb}\n", 1)

        # drop existing branch blocks for this nuclide (they follow immediately)
        rest = text[end:]
        for eb in existing:
            eb_start = rest.index(f"###  {NS}{local(eb)}\n")
            eb_end = rest.index("\n\n\n", eb_start) + 3
            rest = rest[:eb_start] + rest[eb_end:]
        text = text[:start] + block + "\n\n".join(new_blocks) + "\n\n" + rest
        report["isomers_patched"] += 1

    print(f"isomers: {len(isomers)}; patched: {report['isomers_patched']}; "
          f"net new branches: {report['added_branches']}")
    print(f"isomers with no decay text: {len(report['no_text'])}")
    print(f"unparsed text: {report['unparsed']}")
    print(f"unknown mode strings: { {k: len(v) for k, v in report['unknown_modes'].items()} }")
    print(f"branches with daughter outside ontology: {len(report['no_daughter'])} {report['no_daughter'][:10]}")
    print(f"existing branch ratios changed: {report['existing_diff']}")
    if write:
        open(path, "w", encoding="utf-8").write(text)
        print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv[1], "--write" in sys.argv)
