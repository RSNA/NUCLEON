#!/usr/bin/env python3
"""Mirror key physical values onto the class reading of each punned IRI.

NUCLEON asserts its physics on the individual reading of each element and
nuclide IRI. Class-oriented viewers such as NCBO BioPortal, and libraries
that read a punned IRI only as a class (owlready2), therefore show no values.
OWL 2 DL forbids one property IRI from being both a datatype property and an
annotation property, so this script adds a parallel set of annotation
properties, each a sub-property of :displayValue and linked to the datatype
property it mirrors by :mirrorOf, and copies every value across verbatim.
It also adds a one-line :decaySummary to every nuclide, e.g.
"IT 99.9963% → Tc-99; β− 0.0037% → Ru-99".

The datatype properties remain authoritative. The mirrors are generated, never
edited by hand; V6 in nucleon_queries.rq checks that they match, and
scripts/validate.py regenerates every :decaySummary and compares.

The script edits the Turtle text in place, preserving the Protégé layout, and
is idempotent: re-running it replaces earlier output.

Usage: python3 scripts/add_display_annotations.py NUCLEON.ttl [--write]
"""
import re
import sys

import rdflib
from rdflib import RDF, RDFS, Namespace

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")

# (datatype property, display annotation property, label)
MIRRORS = [
    ("atomicNumber", "atomicNumberDisplay", "atomic number"),
    ("massNumber", "massNumberDisplay", "mass number"),
    ("neutronNumber", "neutronNumberDisplay", "neutron number"),
    ("halfLife", "halfLifeDisplay", "half-life in seconds"),
    ("isStable", "isStableDisplay", "is stable"),
    ("metastableState", "metastableStateDisplay", "metastable state"),
    ("naturalAbundance", "naturalAbundanceDisplay", "natural abundance in percent"),
    ("decayWidthValue", "decayWidthValueDisplay", "decay width value"),
    ("decayWidthUnit", "decayWidthUnitDisplay", "decay width unit"),
    ("groupNumber", "groupNumberDisplay", "group number"),
    ("periodNumber", "periodNumberDisplay", "period number"),
]
GENERATED = [d for _, d, _ in MIRRORS] + ["decaySummary"]

DECL_START = "###  " + NS + "displayValue\n"
DECL_END_MARK = "#    Object Properties"


def local(u):
    return str(u).split("#", 1)[1]


def mode_symbol(g, mode):
    """Short NUBASE-style symbol for a decay mode, taken from its altLabels."""
    alts = sorted(str(a) for a in g.objects(mode, SKOS.altLabel))
    greek = [a for a in alts if any(ord(c) > 127 for c in a)]
    if greek:
        return min(greek, key=lambda a: (len(a), a))
    for sym in ("IT", "SF", "CD", "p", "n", "2p", "2n", "e+"):
        if sym in alts:
            return sym
    return str(g.value(mode, RDFS.label))


def decay_summary(g, nuclide):
    """One-line summary of a nuclide's decay, or None if there is nothing to say."""
    branches = sorted(g.objects(nuclide, N.hasDecayBranch),
                      key=lambda b: int(str(b).rsplit("_b", 1)[1]))
    if not branches:
        return "stable" if (nuclide, N.isStable, rdflib.Literal(True)) in g else None
    parts = []
    for b in branches:
        mode = g.value(b, N.decayMode)
        s = mode_symbol(g, mode)
        ratio = g.value(b, N.branchingRatio)
        comment = str(g.value(b, RDFS.comment) or "")
        op = ("~" if "approximate" in comment else
              "<" if "upper limit" in comment else
              ">" if "lower limit" in comment else "")
        if ratio is None:
            s += " ?%"
        else:
            r = float(ratio)
            txt = f"{r:.2g}" if 0 < r < 1e-3 else str(ratio)
            s += f" {op}{txt}%"
        d = g.value(b, N.decayProduct)
        if d is not None:
            s += f" → {local(d)}"
        elif (mode, N.changeInAtomicNumber, None) not in g:
            s += " (variable products)"
        parts.append(s)
    return "; ".join(parts)


def declarations():
    out = [
        DECL_START
        + ":displayValue rdf:type owl:AnnotationProperty ;\n"
        + '              rdfs:comment "Parent of the annotation properties that mirror, on the class reading of a punned element or nuclide IRI, values asserted on its individual reading, so that class-oriented viewers such as BioPortal can display them. Generated from the datatype properties by scripts/add_display_annotations.py; the datatype properties are authoritative. Each sub-property names the property it mirrors with :mirrorOf."@en ;\n'
        + '              rdfs:label "display value"@en .\n',
        "###  " + NS + "mirrorOf\n"
        + ":mirrorOf rdf:type owl:AnnotationProperty ;\n"
        + '          rdfs:comment "Relates a :displayValue annotation property to the datatype property whose values it mirrors."@en ;\n'
        + '          rdfs:label "mirror of"@en .\n',
        "###  " + NS + "decaySummary\n"
        + ":decaySummary rdf:type owl:AnnotationProperty ;\n"
        + "              rdfs:subPropertyOf :displayValue ;\n"
        + '              rdfs:comment "One-line summary of a nuclide\'s decay branches, generated from :hasDecayBranch: mode symbol, branching ratio in percent (? where the source gives none; ~, < or > where the source gives an approximation or limit) and daughter. \\"stable\\" for a stable nuclide with no branch."@en ;\n'
        + '              rdfs:label "decay summary"@en .\n',
    ]
    for src, disp, label in MIRRORS:
        pad = " " * (len(disp) + 2)
        out.append(
            f"###  {NS}{disp}\n"
            f":{disp} rdf:type owl:AnnotationProperty ;\n"
            f"{pad}rdfs:subPropertyOf :displayValue ;\n"
            f"{pad}:mirrorOf :{src} ;\n"
            f'{pad}rdfs:comment "Mirror of :{src} on the class reading, for display. Generated; do not edit."@en ;\n'
            f'{pad}rdfs:label "{label} (display)"@en .\n'
        )
    return "\n\n".join(out) + "\n\n"


def main(path, write):
    g = rdflib.Graph()
    g.parse(path)
    text = open(path, encoding="utf-8").read()

    # 1. remove earlier output (idempotency)
    gen = "|".join(GENERATED)
    text = re.sub(rf"\n[ ]+:(?:{gen}) [^\n]*? ;(?=\n)", "", text)
    text = re.sub(rf"( ;)\n[ ]+:(?:{gen}) [^\n]*? \.\n", " .\n", text)
    if DECL_START in text:
        s = text.index(DECL_START)
        e = text.index("#################################################################\n" + DECL_END_MARK)
        # keep everything after our blocks that is not ours
        blocks = text[s:e].split("\n\n\n")
        keep = [b for b in blocks if not any(f"{NS}{p}\n" in b for p in GENERATED + ["displayValue", "mirrorOf"])]
        text = text[:s] + "\n\n\n".join(keep) + text[e:]

    # 2. declarations, at the end of the annotation-property section
    e = text.index("#################################################################\n" + DECL_END_MARK)
    text = text[:e] + declarations() + text[e:]

    # 3. per-entity mirrors
    entities = [s for s in set(g.subjects(RDF.type, N.Nuclide)) | set(g.subjects(RDF.type, N.Element))
                if (s, RDF.type, rdflib.OWL.Class) in g]
    added = summaries = 0
    for ent in sorted(entities, key=str):
        name = local(ent)
        start = text.index(f"###  {NS}{name}\n")
        end = text.index("\n\n\n", start) + 1
        block = text[start:end]
        pad = " " * (len(name) + 2)
        lines = []
        for src, disp, _ in MIRRORS:
            m = re.search(rf"^[ ]+:{src} (.+?) [;.]$", block, flags=re.M)
            if m:
                lines.append(f"{pad}:{disp} {m.group(1)}")
        if (ent, RDF.type, N.Nuclide) in g:
            summ = decay_summary(g, ent)
            if summ:
                lines.append(f'{pad}:decaySummary "{summ}"@en')
                summaries += 1
        if not lines:
            continue
        assert block.endswith(" .\n"), name
        block = block[:-3] + " ;\n" + " ;\n".join(lines) + " .\n"
        text = text[:start] + block + text[end:]
        added += len(lines)

    print(f"entities: {len(entities)}; annotation assertions: {added}; decay summaries: {summaries}")
    if write:
        open(path, "w", encoding="utf-8").write(text)
        print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv[1], "--write" in sys.argv)
