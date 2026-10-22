#!/usr/bin/env python3
"""Check NUCLEON with two OWL 2 DL reasoners: HermiT and Pellet.

Both reasoners ship with owlready2 (pip install owlready2; Java required).
The ontology is converted to RDF/XML, then:

  Pellet 2.3.1: `consistency`, then `unsat` (classifies the ontology and
                lists unsatisfiable classes); about 10 s each.
  HermiT 1.3.8: consistency and classification through owlready2;
                about 15 min for version 1.5.

Writes validation/reasoners.txt.

Usage: python3 scripts/check_reasoners.py [NUCLEON.ttl] [--skip-hermit]
"""
import glob
import os
import subprocess
import sys
import tempfile
import time

import owlready2
import rdflib


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ttl = args[0] if args else "NUCLEON.ttl"
    skip_hermit = "--skip-hermit" in sys.argv
    tmp = tempfile.mkdtemp()
    owl = os.path.join(tmp, "nucleon.owl")
    g = rdflib.Graph()
    g.parse(ttl)
    g.serialize(owl, format="xml")
    lines = [f"Ontology: {ttl} ({len(g):,} triples)", f"owlready2 {owlready2.VERSION}", ""]

    pellet_dir = os.path.join(os.path.dirname(owlready2.__file__), "pellet")
    pellet_jar = os.path.basename(glob.glob(os.path.join(pellet_dir, "pellet-*.jar"))[0])
    cp = os.path.join(pellet_dir, "*")
    for cmd in ("consistency", "unsat"):
        t0 = time.time()
        out = subprocess.run(["java", "-Xmx6000M", "-cp", cp, "pellet.Pellet", cmd, f"file://{owl}"],
                             capture_output=True, text=True).stdout
        keep = [l for l in out.splitlines() if l.startswith(("Consistent", "Found", "Unsatisfiable"))
                or (cmd == "unsat" and l.strip() and "%" not in l and "Finding" not in l)]
        lines.append(f"Pellet ({pellet_jar}) {cmd}: " + "; ".join(keep) + f"  [{time.time() - t0:.0f} s]")

    if not skip_hermit:
        t0 = time.time()
        w = owlready2.World()
        o = w.get_ontology(f"file://{owl}").load()
        try:
            with o:
                owlready2.sync_reasoner_hermit(w, infer_property_values=False, debug=0)
            unsat = list(w.inconsistent_classes())
            lines.append(f"HermiT (HermiT.jar 1.3.8): consistent; unsatisfiable classes: {len(unsat)}"
                         f"  [{time.time() - t0:.0f} s]")
        except owlready2.OwlReadyInconsistentOntologyError:
            lines.append("HermiT (HermiT.jar 1.3.8): INCONSISTENT")

    report = "\n".join(lines) + "\n"
    print(report)
    os.makedirs("validation", exist_ok=True)
    open(os.path.join("validation", "reasoners.txt"), "w").write(report)
    sys.exit(0 if "INCONSISTENT" not in report and "Consistent: No" not in report else 1)


if __name__ == "__main__":
    main()
