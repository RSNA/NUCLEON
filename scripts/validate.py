#!/usr/bin/env python3
"""Run every query in nucleon_queries.rq against NUCLEON.ttl and print the
row counts, then the content counts reported in Table 1 of the manuscript and
a set of structural checks. Exits non-zero if any validation query (V*) or
structural check reports a violation.

Usage: python3 scripts/validate.py [NUCLEON.ttl] [nucleon_queries.rq]
"""
import re
import sys

import rdflib
from rdflib import OWL, RDF, RDFS, Namespace

sys.path.insert(0, __import__("os").path.dirname(__file__))
from add_display_annotations import decay_summary  # noqa: E402

N = Namespace("https://rsna.org/radlex/nucleon#")
SKOS = Namespace("http://www.w3.org/2004/02/skos/core#")


def split_queries(rq):
    parts = re.split(r"^### (\S+) — (.*?) ###\s*$", rq, flags=re.M)
    for i in range(1, len(parts), 3):
        qid, title, body = parts[i], parts[i + 1], parts[i + 2]
        body = body.split("\n####")[0]
        sparql = "\n".join(l for l in body.splitlines() if not l.lstrip().startswith("#"))
        yield qid, title, sparql.strip()


def main(ttl="NUCLEON.ttl", rq="nucleon_queries.rq"):
    g = rdflib.Graph()
    g.parse(ttl)
    failures = 0

    print("== Queries ==")
    for qid, title, q in split_queries(open(rq, encoding="utf-8").read()):
        rows = list(g.query(q))
        if len(rows) == 1 and len(rows[0]) == 1 and "COUNT" in q.upper() and "violations" in q:
            n = int(rows[0][0])
        else:
            n = len(rows)
        distinct = ""
        if qid in ("Q1", "Q4", "Q5"):
            distinct = f" ({len({r[0] for r in rows})} distinct)"
        print(f"{qid:4} {n:6}{distinct}  {title}")
        if qid.startswith("V") and n:
            failures += 1

    print("\n== Table 1 counts ==")
    nuclides = {s for s in g.subjects(RDF.type, N.Nuclide) if (s, RDF.type, OWL.NamedIndividual) in g}
    elements = {s for s in g.subjects(RDF.type, N.Element) if (s, RDF.type, OWL.NamedIndividual) in g}
    branches = set(g.subjects(RDF.type, N.DecayBranch))
    isomers = {n for n in nuclides if (n, N.metastableState, None) in g}
    stable = {n for n in nuclides if (n, N.isStable, rdflib.Literal(True)) in g}
    with_branch = {n for n in nuclides if (n, N.hasDecayBranch, None) in g}
    classes = set(g.subjects(RDF.type, OWL.Class))
    indivs = set(g.subjects(RDF.type, OWL.NamedIndividual))
    counts = {
        "triples": len(g),
        "named classes": len(classes),
        "named individuals": len(indivs),
        "punned": len(classes & indivs),
        "elements": len(elements),
        "nuclides": len(nuclides),
        "  metastable isomers": len(isomers),
        "  flagged stable": len(stable),
        "  with half-life": sum(1 for n in nuclides if (n, N.halfLife, None) in g),
        "  with spin-parity": sum(1 for n in nuclides if (n, N.spinParity, None) in g),
        "  with >=1 decay branch": len(with_branch),
        "  with >=2 decay branches": sum(1 for n in nuclides if len(list(g.objects(n, N.hasDecayBranch))) >= 2),
        "  isomers with a branch": len(with_branch & isomers),
        "  radioactive, no branch (ground state)": len(nuclides - stable - with_branch - isomers),
        "  radioactive, no branch (isomer)": len((isomers - stable) - with_branch),
        "decay branches": len(branches),
        "  with decay product": sum(1 for b in branches if (b, N.decayProduct, None) in g),
        "  with branching ratio": sum(1 for b in branches if (b, N.branchingRatio, None) in g),
        "  fission (no product by design)": sum(
            1 for b in branches
            if (b, N.decayProduct, None) not in g
            and (g.value(b, N.decayMode), N.changeInAtomicNumber, None) not in g),
        "decay modes": len(set(g.subjects(RDF.type, N.DecayMode))),
        "skos:exactMatch": len(list(g.triples((None, SKOS.exactMatch, None)))),
        "nuclides with dicomCode": len(set(g.subjects(N.dicomCode, None))),
        "language-tagged element names": sum(
            1 for e in elements for l in g.objects(e, SKOS.altLabel) if getattr(l, "language", None)),
    }
    for k, v in counts.items():
        print(f"{k:42} {v:>8,}")

    print("\n== Structural checks ==")
    checks = {
        "branch without a parent": [b for b in branches if (None, N.hasDecayBranch, b) not in g],
        "branch without a decay mode": [b for b in branches if (b, N.decayMode, None) not in g],
        "nuclide without an element": [n for n in nuclides if (None, N.hasNuclide, n) not in g],
        "branching ratio outside 0-100": [b for b in branches for r in g.objects(b, N.branchingRatio)
                                          if not 0 <= float(r) <= 100],
        "deprecated :isMedicalIsotope present": list(g.triples((None, N.isMedicalIsotope, None))),
        ":decaySummary differs from its branches": [
            n for n in nuclides
            if (str(g.value(n, N.decaySummary)) if g.value(n, N.decaySummary) is not None else None)
            != decay_summary(g, n)],
    }
    for k, v in checks.items():
        print(f"{k:42} {len(v):>8}")
        failures += bool(v)

    # Warnings: reported, not failures. The IAEA source lists observed or
    # predicted double-beta / alpha modes for some primordial nuclides that
    # NUCLEON flags as stable (e.g. Bi-209, half-life 2.01e19 y).
    warn = sorted(str(n).split("#")[1] for n in stable if (n, N.hasDecayBranch, None) in g)
    print(f"{'WARN stable-flagged nuclide with a branch':42} {len(warn):>8}  {' '.join(warn)}")

    print("\nPASS" if not failures else f"\nFAIL ({failures})")
    return failures


if __name__ == "__main__":
    sys.exit(1 if main(*sys.argv[1:]) else 0)
