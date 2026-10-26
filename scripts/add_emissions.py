#!/usr/bin/env python3
"""Add radiation emissions from ICRP Publication 107 to NUCLEON.

Scope: the nuclides DICOM codes (CID 18, CID 4020) and every radioactive
member of their decay chains (following :hasDecayBranch / :hasDecayProduct).
Source: ICRP Publication 107, Nuclear Decay Data for Dosimetric Calculations
(Ann. ICRP 38(3), 2008), file ICRP-07.RAD, read from the RadData R package
(data/ICRP_07.RAD.rda; https://github.com/markhogue/RadData).

ICRP-107 gives, for each radionuclide, every emitted radiation with its
energy and yield per decay (summed over the nuclide's decay branches). Each
emission with yield >= 0.1% per decay and energy >= 1 keV becomes an
:Emission individual linked from the nuclide by :hasEmission, with
:hasRadiationType (one of 14 :RadiationType individuals, ICRP-107's radiation
codes), :emissionEnergy (keV; for beta particles the mean energy of the
spectrum) and :emissionYield (per decay; may exceed 1, e.g. two annihilation
photons per positron). Complete spectra, including weaker lines, are in
ICRP-107 itself.

Each nuclide class also gets :emissionSummary, a display annotation listing
its principal photons (yield >= 1%) and particles (yield >= 10%).

ICRP-107 licence: the data files may be used, copied and distributed for
educational, research and not-for-profit purposes, with the authors'
copyright notice (source-data/LICENSE.ICRP-07).

Idempotent: re-running replaces earlier output.

Usage: python3 scripts/add_emissions.py NUCLEON.ttl RAD.rda [--write]
"""
import collections
import math
import re
import sys
from decimal import Decimal

import pyreadr
import rdflib
from rdflib import Namespace

NS = "https://rsna.org/radlex/nucleon#"
N = Namespace(NS)
MIN_YIELD = 1e-3          # 0.1% per decay
MIN_ENERGY_MEV = 1e-3     # 1 keV

# ICRP-107 radiation code -> (individual, label, symbol, comment)
RTYPES = {
    "G":  ("GammaRay", "gamma ray", "γ", "Photon emitted in a nuclear transition."),
    "PG": ("PromptFissionGammaRay", "prompt gamma ray from spontaneous fission", "γ(SF, prompt)", None),
    "DG": ("DelayedFissionGammaRay", "delayed gamma ray from spontaneous fission", "γ(SF, delayed)", None),
    "X":  ("XRay", "X-ray", "X", "Characteristic photon emitted as the atom's electron shells relax."),
    "AQ": ("AnnihilationPhoton", "annihilation photon", "γ±", "511-keV photon from positron annihilation; two per positron."),
    "B+": ("Positron", "positron (beta-plus particle)", "β+", "Energy given is the mean of the continuous spectrum."),
    "B-": ("BetaMinusParticle", "beta-minus particle", "β−", "Energy given is the mean of the continuous spectrum."),
    "BD": ("DelayedFissionBetaParticle", "delayed beta particle from spontaneous fission", "β(SF)", "Energy given is the mean of the continuous spectrum."),
    "IE": ("InternalConversionElectron", "internal conversion electron", "ce", None),
    "AE": ("AugerElectron", "Auger electron", "e(Auger)", None),
    "A":  ("AlphaParticle", "alpha particle", "α", None),
    "AR": ("AlphaRecoilNucleus", "alpha recoil nucleus", "α recoil", "The daughter nucleus recoiling from alpha emission."),
    "FF": ("FissionFragment", "fission fragment", "FF", None),
    "N":  ("Neutron", "neutron", "n", None),
}
ORDER = ["G", "X", "AQ", "PG", "DG", "B-", "B+", "BD", "IE", "AE", "A", "AR", "FF", "N"]
PHOTONS = {"G", "X", "AQ", "PG", "DG"}

DECL_IRIS = ["Emission", "RadiationType", "hasEmission", "hasRadiationType", "emissionEnergy",
             "emissionYield", "emissionSummary"] + [v[0] for v in RTYPES.values()]


def local(u):
    return str(u).split("#", 1)[1]


def sig(x, n=6):
    """Decimal literal with n significant digits, no exponent."""
    d = Decimal(f"{x:.{n}g}")
    t = format(d.normalize(), "f")
    return t if "." in t else t + ".0"


def block(name, pairs):
    head = f"###  {NS}{name}\n:{name} "
    pad = " " * (len(name) + 2)
    lines = []
    for k, (p, objs) in enumerate(pairs):
        lead = head if k == 0 else pad
        objpad = " " * (len(lead.split("\n")[-1]) + len(p) + 1)
        lines.append(lead + p + " " + (" ,\n" + objpad).join(objs))
    return " ;\n".join(lines) + " .\n\n\n"


def remove_previous(text):
    # emission individual blocks and declarations
    for m in list(re.finditer(rf"###  {re.escape(NS)}([^\n]+)\n", text))[::-1]:
        name = m.group(1)
        if re.search(r"_e\d+$", name) or name in DECL_IRIS:
            e = text.index("\n\n\n", m.start()) + 3
            text = text[:m.start()] + text[e:]
    text = text.replace("###  http://www.w3.org/2004/02/skos/core#notation\nskos:notation rdf:type owl:AnnotationProperty .\n\n\n", "")
    # :hasEmission and :emissionSummary lines in nuclide blocks
    lit = r'"(?:[^"\\]|\\.)*"@en'
    text = re.sub(r"\n[ ]+:hasEmission :[^\n]*?(?: ,\n[ ]+:[^\n]*?)* ;(?=\n)", "", text)
    text = re.sub(r" ;\n[ ]+:hasEmission :[^\n]*?(?: ,\n[ ]+:[^\n]*?)* \.\n", " .\n", text)
    text = re.sub(rf"\n[ ]+:emissionSummary {lit} ;(?=\n)", "", text)
    text = re.sub(rf" ;\n[ ]+:emissionSummary {lit} \.\n", " .\n", text)
    return text


def insert_before(text, header, chunk):
    marker = "#################################################################\n#    " + header
    i = text.index(marker)
    return text[:i] + chunk + text[i:]


def main(path, rad_path, write=False):
    g = rdflib.Graph()
    g.parse(path)
    text = remove_previous(open(path, encoding="utf-8").read())

    rad = pyreadr.read_r(rad_path)["ICRP_07.RAD"]
    rad = rad[(rad.prob >= MIN_YIELD) & (rad.E_MeV >= MIN_ENERGY_MEV)]

    # scope: DICOM-coded nuclides and radioactive chain members
    scope, stack = set(), [local(s) for s in g.subjects(N.dicomCode, None)]
    while stack:
        n = stack.pop()
        if n in scope:
            continue
        scope.add(n)
        for b in g.objects(N[n], N.hasDecayBranch):
            d = g.value(b, N.hasDecayProduct)
            if d is not None and (d, N.isStable, rdflib.Literal(True)) not in g:
                stack.append(local(d))

    # ICRP names isomers m, n in its own order; check each named state by half-life
    import radioactivedecay as rd
    # Match each NUCLEON nuclide to the ICRP-107 state of the same nucleus with the
    # closest half-life. ICRP labels isomers m, n in its own order, and older
    # evaluations sometimes swap ground state and isomer (ICRP Ta-178m = NUCLEON Ta-178).
    icrp_all = {str(x) for x in rd.DEFAULTDATA.nuclides}
    icrp_names = set(rad.RN)
    icrp_of, missing, remapped, newer = {}, [], [], []
    for n in sorted(scope):
        base = re.sub(r"m\d?$", "", n)
        cands = [c for c in icrp_all if re.sub(r"[mn]$", "", c) == base and c in icrp_names]
        hl_n = g.value(N[n], N.halfLife)
        if not cands:
            missing.append(n)
            continue
        if hl_n is None:
            best = n if n in cands else None
        else:
            best = min(cands, key=lambda c: abs(math.log10(float(hl_n) / rd.Nuclide(c).half_life("s"))))
            if abs(math.log10(float(hl_n) / rd.Nuclide(best).half_life("s"))) > math.log10(1.5):
                if n in cands:          # same state, newer half-life evaluation
                    newer.append(n)
                    best = n
        if best is None:
            missing.append(n)
            continue
        if best != n:
            remapped.append((n, best))
        icrp_of[n] = best

    # ---- declarations --------------------------------------------------------------------
    obj = (block("hasEmission", [("rdf:type", ["owl:ObjectProperty"]), ("rdfs:domain", [":Nuclide"]),
                                ("rdfs:range", [":Emission"]),
                                ("rdfs:comment", ['"Relates a radionuclide to a radiation it emits, from ICRP Publication 107."@en']),
                                ("rdfs:label", ['"has emission"@en'])])
           + block("hasRadiationType", [("rdf:type", ["owl:ObjectProperty", "owl:FunctionalProperty"]),
                                     ("rdfs:domain", [":Emission"]), ("rdfs:range", [":RadiationType"]),
                                     ("rdfs:comment", ['"The kind of radiation an emission is (gamma ray, X-ray, beta-minus particle, ...)."@en']),
                                     ("rdfs:label", ['"radiation type"@en'])]))
    dat = (block("emissionEnergy", [("rdf:type", ["owl:DatatypeProperty", "owl:FunctionalProperty"]),
                                    ("rdfs:domain", [":Emission"]), ("rdfs:range", ["xsd:decimal"]),
                                    ("rdfs:comment", ['"Energy of the emitted radiation in keV, from ICRP Publication 107. For beta particles (beta-minus, positron, delayed fission beta) it is the mean energy of the continuous spectrum; for all other radiations, the discrete energy."@en']),
                                    ("rdfs:label", ['"emission energy (keV)"@en'])])
           + block("emissionYield", [("rdf:type", ["owl:DatatypeProperty", "owl:FunctionalProperty"]),
                                     ("rdfs:domain", [":Emission"]), ("rdfs:range", ["xsd:decimal"]),
                                     ("rdfs:comment", ['"Number of such radiations emitted per decay of the nuclide (all decay branches together), from ICRP Publication 107. A fraction, not a percentage; it can exceed 1 (two annihilation photons per positron). Emissions with yield below 0.001 or energy below 1 keV are not included."@en']),
                                     ("rdfs:label", ['"emission yield (per decay)"@en'])]))
    ann = block("emissionSummary", [("rdf:type", ["owl:AnnotationProperty"]), ("rdfs:subPropertyOf", [":displayValue"]),
                                    ("rdfs:comment", ['"One-line summary of a radionuclide\'s principal emissions for display: photons with yield of at least 1% and particles with yield of at least 10%, as energy in keV and yield in percent. Generated from :hasEmission."@en']),
                                    ("rdfs:label", ['"emission summary"@en'])])
    cls = (block("Emission", [("rdf:type", ["owl:Class"]),
                              ("rdfs:comment", ['"A radiation emitted in the decay of a radionuclide, with its type, energy and yield per decay. Emission data are from ICRP Publication 107 (Nuclear Decay Data for Dosimetric Calculations, 2008) and are given for the radionuclides DICOM codes (CID 18, CID 4020) and the radioactive members of their decay chains."@en']),
                              ("rdfs:label", ['"emission"@en'])])
           + block("RadiationType", [("rdf:type", ["owl:Class"]),
                                     ("rdfs:comment", ['"A kind of ionizing radiation emitted in radioactive decay. The instances are the radiation types of ICRP Publication 107; each carries ICRP-107\'s radiation code as skos:notation."@en']),
                                     ("rdfs:label", ['"radiation type"@en'])]))
    ind = ""
    for code in ORDER:
        iri, label, sym, comment = RTYPES[code]
        p = [("rdf:type", ["owl:NamedIndividual", ":RadiationType"])]
        if comment:
            p.append(("rdfs:comment", [f'"{comment}"@en']))
        p += [("rdfs:label", [f'"{label}"@en']), ("skos:altLabel", [f'"{sym}"']),
              ("skos:notation", [f'"{code}"'])]
        ind += block(iri, p)
    text = insert_before(text, "Data properties", obj)
    text = insert_before(text, "Classes", dat)
    ann += "###  http://www.w3.org/2004/02/skos/core#notation\nskos:notation rdf:type owl:AnnotationProperty .\n\n\n"
    text = insert_before(text, "Object Properties", ann)
    text = insert_before(text, "Individuals", cls)
    text = insert_before(text, "General axioms", ind)

    # ---- emissions per nuclide ------------------------------------------------------------
    total, counts = 0, collections.Counter()
    for n in sorted(icrp_of):
        rows = rad[rad.RN == icrp_of[n]].copy()
        rows["o"] = rows.code_AN.map(ORDER.index)
        rows = rows.sort_values(["o", "E_MeV"])
        ems, summ_ph, summ_pa = [], [], []
        for k, r in enumerate(rows.itertuples(), 1):
            iri, label, sym, _ = RTYPES[r.code_AN]
            e_kev = sig(r.E_MeV * 1000)
            y = sig(r.prob)
            name = f"{n}_e{k}"
            ems.append(block(name, [("rdf:type", ["owl:NamedIndividual", ":Emission"]),
                                    (":hasRadiationType", [f":{iri}"]),
                                    (":emissionEnergy", [e_kev]), (":emissionYield", [y]),
                                    ("rdfs:label", [f'"{n} {label} {e_kev.rstrip("0").rstrip(".")} keV"@en'])]))
            counts[r.code_AN] += 1
            pct = r.prob * 100
            if r.code_AN in PHOTONS and pct >= 1:
                summ_ph.append((pct, f"{sym} {float(e_kev):.4g} keV ({pct:.3g}%)"))
            elif r.code_AN in ("A", "B-", "B+") and pct >= 10:
                summ_pa.append((pct, f"{sym} {'mean ' if r.code_AN != 'A' else ''}{float(e_kev):.4g} keV ({pct:.3g}%)"))
        if not ems:
            continue
        total += len(ems)
        summary = "; ".join(s for _, s in sorted(summ_ph, reverse=True)[:6] + sorted(summ_pa, reverse=True)[:3])
        # nuclide block: add :hasEmission and :emissionSummary
        s = text.index(f"###  {NS}{n}\n")
        e = text.index("\n\n\n", s) + 1
        blk = text[s:e]
        pad = " " * (len(n) + 2)
        objpad = pad + " " * len(":hasEmission ")
        add = f"{pad}:hasEmission " + (" ,\n" + objpad).join(f":{n}_e{k}" for k in range(1, len(ems) + 1))
        if summary:
            add += f' ;\n{pad}:emissionSummary "{summary}"@en'
        blk = blk[:-3] + " ;\n" + add + " .\n"
        text = text[:s] + blk + text[e:]
        # emission blocks after the nuclide's own block and its branch blocks
        s = text.index(f"###  {NS}{n}\n")
        e = text.index("\n\n\n", s) + 3
        while True:
            nxt = text.find("###  " + NS, e)
            nm = text[nxt + 5 + len(NS):text.index("\n", nxt)] if nxt != -1 else ""
            if not nm.startswith(n + "_b"):
                break
            e = text.index("\n\n\n", nxt) + 3
        text = text[:e] + "".join(ems) + text[e:]

    print(f"nuclides in scope: {len(scope)}; with emissions: {len(scope) - len(missing)}; "
          f"not in ICRP-107: {missing}")
    print(f"emissions: {total}; by type: {dict(counts)}")
    print(f"ICRP-107 state used under another name: {remapped}")
    print(f"same state, half-life re-evaluated since ICRP-107 by >1.5x: {newer}")
    if write:
        open(path, "w", encoding="utf-8").write(text)
        print(f"wrote {path}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--write"]
    main(*args, write="--write" in sys.argv)
