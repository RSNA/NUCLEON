#!/usr/bin/env python3
"""Compare every NUCLEON parent–daughter link with independent nuclear data.

NUCLEON derives each daughter from the parent's Z and A and the decay mode's
change in Z and A (validation query V1 confirms only that this arithmetic was
applied consistently). This script checks the derived links against two
sources that list daughters explicitly:

  * ICRP Publication 107 (2008): 1,252 radionuclides, read from the
    `radioactivedecay` package (dataset icrp107_ame2020_nubase2020).
    Used for every parent it covers: the nuclear-medicine subset.
  * ENSDF decay datasets, read from the `paceENSDF` package. Used for
    parents ICRP-107 does not cover. The package holds alpha, beta-minus and
    electron-capture/beta-plus datasets only, so other modes of those
    parents cannot be checked.

Comparison is by daughter nucleus (Z, A). NUCLEON resolves every daughter to
its ground state; where the reference names an isomeric daughter (e.g.
Mo-99 -> Tc-99m) the pair is counted as agreeing and flagged.

The IAEA LiveChart export that NUCLEON was built from is itself derived in
part from ENSDF, so ENSDF agreement is a weaker test of the source data than
ICRP-107 agreement; both are fully independent tests of the derivation code.

Also compares, for ICRP-107 parents, half-lives (tests unit normalisation to
seconds) and branching ratios summed per decay-mode family.

Usage: python3 scripts/compare_daughters.py [NUCLEON.ttl] [outdir]
Requires: pip install radioactivedecay==0.6.1 paceENSDF==0.6.3 rdflib
"""
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

import rdflib
from rdflib import RDF, RDFS, Namespace

N = Namespace("https://rsna.org/radlex/nucleon#")

FAMILY = {  # NUCLEON decay mode -> family comparable across sources
    "AlphaDecay": "alpha",
    "BetaMinusDecay": "beta-",
    "BetaPlusDecay": "EC/beta+", "ElectronCapture": "EC/beta+", "PositronEmission": "EC/beta+",
    "IsomericTransition": "IT",
    "SpontaneousFission": "SF",
}
ICRP_FAMILY = {"α": "alpha", "β-": "beta-", "β+ & EC": "EC/beta+", "EC": "EC/beta+", "β+": "EC/beta+",
               "IT": "IT", "SF": "SF"}
ENSDF_FAMILY = {"alphaDecay": "alpha", "betaMinusDecay": "beta-", "electronCaptureBetaPlusDecay": "EC/beta+"}
ENSDF_FAMILIES = set(ENSDF_FAMILY.values())

UNIT_S = {"ps": 1e-12, "ns": 1e-9, "us": 1e-6, "μs": 1e-6, "ms": 1e-3, "s": 1, "m": 60, "h": 3600,
          "d": 86400, "y": 365.2422 * 86400}


def local(u):
    return str(u).split("#", 1)[1]


def split(name):
    """'Tc-99m' -> ('Tc', 99, 1); 'Al-26m2' -> ('Al', 26, 2)."""
    m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m(\d?))?", name)
    return m.group(1), int(m.group(2)), (int(m.group(4) or 1) if m.group(3) else 0)


def load_nucleon(path):
    g = rdflib.Graph()
    g.parse(path)
    z_of = {str(g.value(e, RDFS.label)): int(g.value(e, N.atomicNumber))
            for e in g.subjects(RDF.type, N.Element) if (e, N.atomicNumber, None) in g}
    nuc = {}
    for n in g.subjects(N.massNumber, None):
        name = local(n)
        branches = []
        for b in sorted(g.objects(n, N.hasDecayBranch), key=lambda b: int(str(b).rsplit("_b", 1)[1])):
            mode = local(g.value(b, N.decayMode))
            d = g.value(b, N.decayProduct)
            r = g.value(b, N.branchingRatio)
            branches.append({"branch": local(b), "mode": mode, "family": FAMILY.get(mode, "other"),
                             "daughter": local(d) if d is not None else None,
                             "ratio": float(r) if r is not None else None})
        hl = g.value(n, N.halfLife)
        nuc[name] = {"branches": branches, "halfLife": float(hl) if hl is not None else None,
                     "dicom": (n, N.dicomCode, None) in g,
                     "excitation_keV": float(g.value(n, N.excitationEnergy)) if (n, N.excitationEnergy, None) in g else None}
    return nuc, z_of


def za(name, z_of):
    sym, a, _ = split(name)
    return z_of[sym], a


def load_icrp():
    import radioactivedecay as rd
    data = rd.DEFAULTDATA
    out = {}
    for n in data.nuclides:
        n = str(n)
        nuc = rd.Nuclide(n)
        prog = nuc.progeny()
        if not prog:
            continue
        sym, rest = n.split("-")
        m = re.fullmatch(r"(\d+)([mnp]?)", rest)
        state = {"": "", "m": "m", "n": "m2", "p": "m3"}[m.group(2)]
        key = f"{sym}-{m.group(1)}{state}"
        out[key] = {
            "daughters": [(str(p), ICRP_FAMILY.get(str(mode), "other"), float(bf))
                          for p, mode, bf in zip(prog, nuc.decay_modes(), nuc.branching_fractions())],
            "halfLife": float(nuc.half_life("s")),
            "rdname": n,
        }
    return out, data.dataset_name


def load_ensdf():
    import paceENSDF
    root = os.path.join(os.path.dirname(paceENSDF.__file__), "ENSDF_JSON")
    by_nucleus = defaultdict(lambda: defaultdict(list))  # (sym, A) -> level energy -> [(Z', A', family)]
    for f in glob.glob(os.path.join(root, "*.json")):
        d = json.load(open(f, encoding="utf-8"))
        fam = ENSDF_FAMILY.get(d["decayMode"])
        if fam is None:
            continue
        sym = re.match(r"[A-Za-z]+", d["parentID"]).group(0)
        lev = d["levelEnergyParentDecay"]
        m = re.match(r"[\d.]+", lev)
        e = float(m.group(0)) if m else 0.0
        if lev != "0" and e == 0.0:
            e = 1e-3  # unknown offset level ("0.0+X") is still an excited state
        by_nucleus[(sym, int(d["parentAtomicMass"]))][e].append(
            (int(d["daughterAtomicNumber"]), int(d["daughterAtomicMass"]), fam, os.path.basename(f)))
    return by_nucleus


def ensdf_for(name, nuc_names, ensdf):
    """ENSDF datasets for a NUCLEON nuclide, and how the level was matched."""
    sym, a, state = split(name)
    levels = ensdf.get((sym, a))
    if not levels:
        return None, None
    if state == 0:
        return (levels.get(0.0), "ground state") if 0.0 in levels else (None, None)
    excited = sorted(e for e in levels if e > 0)
    if not excited:
        return None, None
    # match the isomer to the ENSDF parent level with the closest excitation energy
    ex = nuc_names[name].get("excitation_keV")
    if ex is not None:
        best = min(excited, key=lambda e: abs(e - ex))
        if abs(best - ex) <= max(2.0, 0.01 * ex):
            return levels[best], f"isomer level {best:g} keV (NUCLEON {ex:g} keV)"
        return None, None
    isomers = sorted(s for s in (split(x) for x in nuc_names) if s[0] == sym and s[1] == a and s[2] > 0)
    if len(isomers) == len(excited):
        return levels[excited[state - 1]], f"isomer level {excited[state - 1]:g} keV (matched by order)"
    return None, None


def remap_icrp_isomers(icrp, nuc):
    """ICRP-107 labels a nucleus's isomers m, n in its own order, which need not
    match NUCLEON's m, m2, ... (e.g. ICRP Hf-178m = NUCLEON Hf-178m2, 31 y).
    Re-key each ICRP isomer to the NUCLEON isomer of the same nucleus with the
    closest half-life (log scale). Returns (remapped dict, list of renamings)."""
    import math
    out, renamed = {}, []
    for k, v in icrp.items():
        sym, a, state = split(k)
        if state == 0:
            out[k] = v
            continue
        cands = [n for n in nuc if split(n)[:2] == (sym, a) and split(n)[2] > 0 and nuc[n]["halfLife"]]
        if not cands:
            out[k] = v
            continue
        best = min(cands, key=lambda n: abs(math.log10(nuc[n]["halfLife"] / v["halfLife"])))
        if best != k:
            renamed.append((k, best))
        out[best] = v
    return out, renamed


def main(ttl="NUCLEON.ttl", outdir="validation"):
    nuc, z_of = load_nucleon(ttl)
    icrp, icrp_name = load_icrp()
    icrp, renamed = remap_icrp_isomers(icrp, nuc)
    ensdf = load_ensdf()

    rows, ref_only, hl_rows, br_rows = [], [], [], []
    for name in sorted(nuc, key=lambda x: (z_of[split(x)[0]], split(x)[1], split(x)[2])):
        info = nuc[name]
        if not info["branches"] and name not in icrp:
            continue
        # ---- choose reference -------------------------------------------------------
        if name in icrp:
            source, match_note = "ICRP-107", ""
            ref = [(za(re.sub(r"(\d+)[mnp]$", r"\1", p), z_of) if p != "SF" else None, fam,
                    p if p != "SF" else None, bf) for p, fam, bf in icrp[name]["daughters"]]
            ref_families = {r[1] for r in ref}
            covered = None  # all families
        else:
            ds, match_note = ensdf_for(name, nuc, ensdf)
            if ds is None:
                source, ref, ref_families = "none", [], set()
            else:
                source = "ENSDF"
                ref = [((z2, a2), fam, f"Z={z2} A={a2}", None) for z2, a2, fam, _ in ds]
                ref_families = {r[1] for r in ref}
            covered = ENSDF_FAMILIES
        ref_za = {r[0] for r in ref if r[0] is not None}

        # ---- each NUCLEON branch ---------------------------------------------------
        for b in info["branches"]:
            d = b["daughter"]
            ref_names = sorted({r[2] for r in ref if r[1] == b["family"] and r[2]})
            if d is None:
                verdict = "no daughter asserted (by design)" if b["family"] == "SF" or b["mode"].endswith("Fission") \
                    else "no daughter asserted"
            elif source == "none":
                verdict = "not checkable: parent in neither reference"
            elif covered is not None and b["family"] not in covered:
                verdict = f"not checkable: {b['mode']} not in ENSDF package"
            elif za(d, z_of) in {r[0] for r in ref if r[1] == b["family"] and r[0]}:
                isomeric = any(r[0] == za(d, z_of) and r[2] and re.search(r"\d[mnp]$", r[2]) for r in ref)
                verdict = "agree (reference daughter is an isomer; NUCLEON gives ground state)" if isomeric \
                    else "agree"
            elif za(d, z_of) in ref_za:
                verdict = "agree on daughter; reference assigns a different mode"
            elif b["family"] not in ref_families:
                verdict = "branch absent from reference"
            else:
                verdict = "DISAGREE"
            rows.append({"parent": name, "branch": b["branch"], "mode": b["mode"], "ratio_pct": b["ratio"],
                         "nucleon_daughter": d or "", "reference": source, "reference_level": match_note,
                         "reference_daughters_same_mode": "; ".join(ref_names), "verdict": verdict,
                         "dicom_coded": info["dicom"]})

        # ---- reference pairs NUCLEON lacks ----------------------------------------
        nuc_za = {za(b["daughter"], z_of) for b in info["branches"] if b["daughter"]}
        nuc_fams = {b["family"] for b in info["branches"]}
        for r in ref:
            if r[0] is None or r[0] in nuc_za:
                continue
            why = ("NUCLEON has no branch for this parent" if not info["branches"] else
                   "NUCLEON has this mode with a different daughter" if r[1] in nuc_fams else
                   "NUCLEON lacks this mode")
            ref_only.append({"parent": name, "reference": source, "mode_family": r[1],
                             "reference_daughter": r[2],
                             "reference_branching_pct": "" if r[3] is None else round(100 * r[3], 6),
                             "classification": why})

        # ---- half-life and branching ratio vs ICRP-107 ----------------------------
        if source == "ICRP-107":
            h1, h2 = info["halfLife"], icrp[name]["halfLife"]
            rel = abs(h1 - h2) / h2 if h1 and h2 else None
            hl_rows.append({"nuclide": name, "nucleon_s": h1, "icrp107_s": h2,
                            "rel_diff": None if rel is None else round(rel, 6),
                            "flag": "" if rel is not None and rel <= 0.05 else
                            ("no NUCLEON half-life" if h1 is None else "differs >5%")})
            fams = {b["family"] for b in info["branches"]} | {r[1] for r in ref}
            for fam in sorted(fams - {"other"}):
                fb = [b for b in info["branches"] if b["family"] == fam]
                if fam == "EC/beta+" and any(b["mode"] == "BetaPlusDecay" for b in fb):
                    # the source's EC+B+ total already includes separately listed B+ / EC parts
                    fb = [b for b in fb if b["mode"] == "BetaPlusDecay"]
                n_pct = sum(b["ratio"] or 0 for b in fb)
                n_known = any(b["ratio"] is not None for b in fb)
                i_pct = 100 * sum(r[3] for r in ref if r[1] == fam)
                diff = abs(n_pct - i_pct)
                br_rows.append({"nuclide": name, "mode_family": fam,
                                "nucleon_pct": round(n_pct, 6) if n_known else "",
                                "icrp107_pct": round(i_pct, 6),
                                "abs_diff_pct": round(diff, 6) if n_known else "",
                                "flag": "" if not n_known or diff <= 1 else "differs >1 percentage point"})

    os.makedirs(outdir, exist_ok=True)

    def write(fname, data):
        with open(os.path.join(outdir, fname), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys()))
            w.writeheader()
            w.writerows(data)

    write("daughter_comparison.csv", rows)
    write("daughter_comparison_reference_only.csv", ref_only)
    write("halflife_vs_icrp107.csv", hl_rows)
    write("branching_vs_icrp107.csv", br_rows)

    # ---- summary ---------------------------------------------------------------------
    linked = [r for r in rows if r["nucleon_daughter"]]
    checked = [r for r in linked if r["verdict"].startswith(("agree", "DISAGREE", "branch absent"))]
    by_src = Counter((r["reference"], r["verdict"]) for r in linked)
    print(f"ICRP dataset: {icrp_name}; ICRP radionuclides: {len(icrp)}; "
          f"in NUCLEON: {sum(1 for k in icrp if k in nuc)}")
    print(f"NUCLEON branches: {len(rows)}; with daughter: {len(linked)}; checked against a reference: {len(checked)}")
    for (src, v), c in sorted(by_src.items()):
        print(f"  {src:9} {c:6}  {v}")
    print("ICRP isomers re-keyed by half-life:", renamed)
    print("ICRP radionuclides absent from NUCLEON:", sorted(k for k in icrp if k not in nuc))
    print("reference-only pairs:", Counter((r["reference"], r["classification"]) for r in ref_only))
    print("half-life flags:", Counter(r["flag"] for r in hl_rows))
    print("branching flags:", Counter(r["flag"] for r in br_rows))
    return rows, ref_only, hl_rows, br_rows


if __name__ == "__main__":
    main(*sys.argv[1:])
